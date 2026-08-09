"""Project reference region markup onto a query via the C++ hot path.

Takes a parsed mmseqs hit plus the reference entry for the matched scaffold and
returns an AIRR-style record dict for the query.

Junction handling follows AIRR strictly: ``junction`` spans the conserved Cys104
through the [FW]118 that opens FR4; ``junction_aa`` starts with C and ends with
F/W for a canonical rearrangement. A junction is reported **even when not
canonical** (out-of-frame, missing the conserved residues). For an out-of-frame
junction (V and J in different frames) the amino-acid translation inserts 1-2 N
bases after the V germline end to restore the J frame; the codon that then
contains an inserted N is rendered as ``_``. The V/J split inside the junction is
located from the transferred ``v_sequence_end`` / ``j_sequence_start``.
"""

from __future__ import annotations

import math

from .. import _markup
from ..refbuild.constant import isotype_class
from ..refbuild.translate import translate, aa_coords_from_nt, detect_coding_frame
from .cigar import build_cigar, segment_cigars
from .reference import RefEntry, REGIONS

__all__ = ["transfer_hit", "AIRR_COLUMNS"]

# Output column order (AIRR-compatible subset + locus). ``mmseqs2_score``/``_evalue``/
# ``_identity`` carry the alignment quality of the chosen scaffold hit. The score is the
# mmseqs2 bit score over the *whole* V+J scaffold (not a per-segment AIRR ``v_score``), so
# it is named after its source; callers use it to rank references and filter weak hits.
AIRR_COLUMNS = (
    ["sequence_id", "sequence", "locus", "v_call", "d_call", "d2_call", "j_call",
     # Constant region. `c_call` names the CH1 exon(s) of the winning `J + C` scaffold; `c_class` is
     # the ISOTYPE CLASS (IGHG / IGHM / IGHA / ...). Report the class, never the subclass: IGHG1-4 are
     # ~95 % identical over CH1, so the top gene is a coin-flip between them -- it ties on 26.7 % of
     # real reads -- while the top class is unique on every one. Both empty on a V-J scaffold hit.
     "c_call", "c_class",
     "mmseqs2_score", "mmseqs2_evalue", "mmseqs2_identity",
     # Alignment GEOMETRY on the V-J scaffold, 1-based inclusive. A bit score is a scalar and
     # confounds read length with identity; these let a consumer ask the sharper question
     # "is the alignment ANCHORED?" -- i.e. does it run to the end of the read, or stop exactly
     # at a scaffold boundary (so the unaligned tail is a linker / adapter / C-region absent
     # from the scaffold)? A chance alignment stops mid-read AND mid-scaffold.
     "mmseqs2_qstart", "mmseqs2_qend", "mmseqs2_qlen",
     "mmseqs2_tstart", "mmseqs2_tend", "mmseqs2_tlen",
     # The scaffold's own V/J boundary. A read that runs off the V into the non-templated CDR3
     # STOPS at t_vend -- that is an explained clip, not a partial alignment. Without these two
     # numbers a consumer cannot tell "ran out of germline" from "stopped for no reason".
     "mmseqs2_t_vend", "mmseqs2_t_jstart",
     # End of the V-J part of the scaffold: its full length for a V-J scaffold, the J length for a
     # `J + C` scaffold. `tstart >= t_vjend` means the alignment lies WHOLLY inside the constant
     # region -- real receptor mRNA, but with no V(D)J and therefore no clonotype.
     "mmseqs2_t_vjend",
     "rev_comp", "productive",
     # AIRR alignment fields. `sequence_alignment`/`germline_alignment` are the mmseqs aligned
     # query and scaffold strings (coding strand) -- the scaffold IS `V + N*pad + J [+ C]`, so its
     # non-templated stretch reads as N, exactly what AIRR's germline_alignment expects. The flags
     # and identity are computed during translation and were previously discarded into `productive`.
     "stop_codon", "vj_in_frame", "v_identity",
     "sequence_alignment", "germline_alignment",
     # Per-segment CIGARs. Each is the sub-walk of the one query->scaffold alignment whose target
     # falls in that segment's germline range, with the rest of the read soft-clipped (see cigar.py).
     "v_cigar", "j_cigar", "c_cigar",
     # SHM in each segment's OWN germline frame: `G45A,C112T` -- germline base, 1-based germline
     # position, read base. Same alignment walk as the cigars, so it costs nothing extra, and it is
     # what a lineage/SHM tool consumes. `germline_alignment` already CONTAINED this (verified: the
     # germline it reports matches the shipped allele on 28,365 of 28,365 bulk IG reads), but
     # recovering it needs arda's scaffold geometry -- a consumer that just diffs the two alignment
     # strings attributes 20.1 % of the mismatches it finds to a germline that is not there, because
     # the scaffold's N-pad and C region are in the same strings. Scoped to V and J STRUCTURALLY:
     # the pad is not a segment, so a junction position has no germline coordinate to be filed under
     # and cannot enter the list by any code path (see cigar.py's module docstring).
     # Germline coordinates (1-based, in the V or J allele). The scaffold's V part is the V germline
     # verbatim (target pos == V-germline pos) and its J part is the full J allele, so these fall
     # straight out of the target span with a constant offset -- no separate lookup.
     "v_germline_start", "v_germline_end", "j_germline_start", "j_germline_end",
     "v_sequence_start", "v_sequence_end",
     "d_sequence_start", "d_sequence_end", "d2_sequence_start", "d2_sequence_end",
     # D germline coordinates (1-based, in the D allele) and CIGAR, emitted only when the D
     # call is a single allele -- a byte-identical-gene ambiguity list has no one germline to
     # anchor to. The D alignment is gapless, so the cigar is a single M run.
     "d_germline_start", "d_germline_end", "d_cigar",
     "d2_germline_start", "d2_germline_end", "d2_cigar",
     # AIRR `d_support`: the E-value of the D alignment. This is the quantity the call is
     # gated on (see `_d_min_score`), so shipping it lets a caller re-threshold downstream.
     "d_support", "d2_support",
     "j_sequence_start", "np1", "np2", "np3", "junction", "junction_aa"]
    + [c for r in REGIONS for c in (f"{r}_start", f"{r}_end", r, f"{r}_aa")]
    # ⛔ APPENDED AT THE END, deliberately. These are NON-AIRR-schema extras, and adding them
    # mid-list (where they first landed, after `c_cigar`) shifted every column from
    # `v_germline_start` onward -- silently breaking any consumer that reads the shipped set BY
    # POSITION. It is the same rule `airr_header(extra_columns)` already states for
    # `junction_quality`: new columns go last, so the shipped prefix never moves.
    + ["v_mutations", "j_mutations"]
    # ⛔ The JUNCTION BOUNDARY in GERMLINE coordinates, emitted so a downstream consumer can do
    # allele re-assignment and SHM calling WITHOUT arda's reference. `v_mutations` positions are
    # 1-based in the called V allele and `j_mutations` positions are 1-based in the called J allele,
    # but both lists span the junction: the V germline's 3' tail and the J germline's 5' head are
    # inside it, so chew-back and non-templated N/P bases appear as substitutions against a germline
    # that does not template them. Splitting framework from junction needs exactly these two
    # numbers, and until now they existed only in `cdr3_anchors.tsv` inside the reference.
    #
    # `v_anchor_nt` is the 0-based offset of the Cys104 codon in the called V allele's germline;
    # `j_anchor_nt` is the 0-based offset of the [FW]118 codon in the called J allele's. So a
    # `v_mutations` entry at 1-based position p is junction-internal iff p > v_anchor_nt, and a
    # `j_mutations` entry at p is junction-internal iff p <= j_anchor_nt + 3.
    #
    # ⚠ Emitted RAW, exactly as the reference stores them, rather than pre-classified: frequency
    # and position together decide what a variant is, and that is the consumer's call. Measured on
    # SRR5233636 (a TRA amplicon, 500 k reads, where TCRs cannot hypermutate so every entry is
    # spurious): 1.046 V and 1.658 J entries per read, 86.2 % of J entries at J position <= 10, and
    # the high-frequency ones sit INSIDE the junction (TRAV8-6*01 positions 281/282 at 0.88 against
    # its anchor at 270; TRAJ8*01 position 1 at 0.67 against its anchor at 26) -- i.e. an allele
    # difference in the templated V/J tail, not somatic mutation and not N/P diversity.
    + ["v_anchor_nt", "j_anchor_nt"]
)

_VSIDE = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

# D-segment mapping. D germlines are short (~8-31 nt) and trimmed on both ends, so
# they are mapped by gapless local alignment (``_markup.d_local_align``, net score
# match=+1/mismatch=-1) of every locus D allele against the V..J interior of the
# junction — not via the mmseqs scaffold DB.
#
# A raw score floor cannot work here, and the old `_D_MIN_SCORE = 6` was measurably wrong:
# six exact bases is strong evidence in a 10 nt mouse TRB interior scanned against 2 D
# alleles, and near-worthless in a 34 nt IGH interior scanned against 48. On OLGA-generated
# junctions whose D was excised entirely, a score of 6 still called a D in 65 % of IGH
# interiors (12 % TRB, 11 % TRD). The optimal raw floor is 9 for human IGH, 7 for human
# TRB and TRD, 6 for mouse TRB — four knobs, and more per organism.
#
# So gate on an E-value instead: the expected number of chance matches at least this good,
# given the interior length `m` and the total D-database length `n`. For match=+1/
# mismatch=-1 over four equiprobable bases, lambda solves `x/4 + 3/(4x) = 1` with
# `x = e^lambda`, giving `x = 3`; K is the usual ungapped constant.
#
#     E = K * m * n * 3^-S      ->      S_min = ceil( ln(K*m*n/E) / ln 3 )
#
# One knob for every locus and organism. Calibrated on OLGA ground truth (human IGH/TRB/TRD
# via vdjrearm, mouse TRB) with the D-excised interior as the null. At E <= 0.2:
#
#     locus         call rate   D-gene accuracy   false call on D-excised null
#     human IGH        81 %          95 %                    1 %
#     human TRB        47 %          88 %                    6 %
#     human TRD        71 %          99 %                    5 %
#     mouse TRB        67 %          93 %                    4 %
#
# The same gate serves the second (D-D) segment: a non-overlapping second D at E <= 0.2 is
# a false positive in 0-1 % of true single-D junctions, and is found in 61 % of injected
# IGH tandem D-D. The TR loci are far lower -- 15 % (human TRB) and 13 % (mouse TRB) on
# tandems injected in genomic order (TRBD1 5' of TRBD2, the only producible orientation).
# That is a length limit, not an orientation one: sensitivity is ~0 when either D survives
# under 6 nt and 57-72 % once both survive 7+, and trimming usually leaves less.
#
# The absolute E is not exactly calibrated -- insertions are Markov, not uniform ACGT, and
# K is a literature constant -- but it is monotone in S and it is what collapses four
# per-locus thresholds into one. Treat it as an m*n-corrected score, and the 0.2 as an
# empirical operating point (5 % false calls), not as a p-value.
#
# AMINO-ACID INPUT uses the same machinery over a different alphabet: `d_local_align` is a
# plain character comparison, so a D germline translated in its three frames aligns against
# an aa interior unchanged. Only lambda moves. In general lambda = ln((1-p)/p) where p is the
# chance two residues match -- p = 1/4 recovers ln 3 above. For amino acids p is NOT 1/20:
# N-region inserts and D germlines are both G/S/Y-rich, and the measured p over real middles
# x real D frames is 0.0613, giving lambda = 2.7285 (uniform would say 2.9444).
#
# The aa gate is stricter (E <= 0.05) because a 20-letter alphabet with a ~20-40 aa database
# leaves the E-value badly under-calibrated at small n -- at E <= 0.2 the measured false-call
# rate on a composition-preserving shuffled null reaches 13 % for IGH. At E <= 0.05 with the
# floor below (null = shuffled middle; positives = D surviving >= 6 nt):
#
#     locus         call rate   D-gene accuracy   false call on null   ambiguous
#     human IGH        69 %          99 %                2.0 %            5 %
#     human TRB         8 %         100 %                0.3 %            0 %
#     human TRD        12 %         100 %                0.0 %            0 %
#     mouse TRB        11 %         100 %                0.7 %           58 %
#
# Only IGH carries enough surviving D to see in protein. The TR loci fire on ~1 record in 10
# and are right when they do; mouse TRBD1/TRBD2 translate to near-identical poly-glycine, so
# most of those calls are honestly reported as a two-gene ambiguity list rather than a guess.
#
# Out of model, on the real GenBank fixtures (protein reconstructed from each read, compared
# with that same read's nucleotide D call): human IGH calls a D on 36 % of records where the
# nucleotide caller manages 68 %, and agrees with it on 98 % of those, 5 % ambiguous. Somatic
# hypermutation and codon degeneracy cost about half the recall and almost none of the
# precision -- which is the trade an aa-only consumer should expect.
_D_LAMBDA = math.log(3.0)                  # nt: p = 1/4
_D_AA_LAMBDA = 2.7285                      # aa: p = 0.0613, measured (see above)
_D_K = 0.33
_D_MAX_EVALUE = 0.2
_D_AA_MAX_EVALUE = 0.05
_D_SCORE_FLOOR = 4      # three exact bases is never evidence, whatever the arithmetic says
_D_AA_SCORE_FLOOR = 4   # governs only the TR loci, whose aa D database is 22-38 residues


def _d_gate(seqtype: str) -> tuple[float, float, int]:
    """``(lambda, default max E-value, score floor)`` for this alphabet."""
    if seqtype == "aa":
        return _D_AA_LAMBDA, _D_AA_MAX_EVALUE, _D_AA_SCORE_FLOOR
    return _D_LAMBDA, _D_MAX_EVALUE, _D_SCORE_FLOOR


def _d_db_nt(d_germlines) -> int:
    """Total residues in a locus D germline set (the `n` of the E-value).

    For an aa reference the set holds three translated frames per allele, so `n` counts all
    three: three frames really are three chances to match.
    """
    return sum(len(seq) for _, seq in d_germlines)


def _allowed_d(d_germlines, j_call: str):
    """Drop D genes lying 3' of the called J -- deletional joining cannot reach them.

    The TRB locus runs TRBD1 - TRBJ1 cluster - TRBC1 - TRBD2 - TRBJ2 cluster - TRBC2, and
    V(D)J recombination deletes the DNA between the joined segments. TRBD2 therefore cannot
    join any TRBJ1: it sits downstream of the entire J1 cluster. This is genomic order, not
    a usage preference, and it holds in every species with this architecture (human, mouse,
    rat, rhesus). IGH and TRD place every D 5' of every J, so nothing is excluded there.

    Without this, TRBD2 (16 nt) simply outscores TRBD1 (12 nt) on noise: 17 % of real human
    TRB J1-cluster records were assigned an impossible TRBD2, at E-values (median 0.096)
    sitting in the same band as chance hits, versus 0.014 for the genuinely producible
    TRBJ2 x TRBD2. An ambiguous J spanning both clusters excludes nothing.
    """
    # ⛔ ORPHONS FIRST, and unconditionally. IMGT ships `/OR` D genes -- `IGHD.../OR15-...` sit on
    # CHROMOSOME 15, outside the IGH locus, and cannot rearrange at all. They are not a usage
    # preference to down-weight; they are not producible. Measured on a real bulk IGH library:
    # **11 of 11 tandem D-D calls named `IGHD2/OR15-2a*01,IGHD2/OR15-2b*01` as their second D**, so
    # the entire tandem-D signal in that library was this one vocabulary artifact. Excluding them
    # leaves 639 of 650 single-D calls untouched, moves 9 to a real rearrangeable gene, loses 2,
    # and drops tandem D-D from 11 to 2.
    d_germlines = [(a, s) for a, s in d_germlines if "/OR" not in a]
    genes = {a.split("*")[0] for a in j_call.split(",") if a.strip()}
    if not genes or not all(g.startswith("TRBJ1-") for g in genes):
        return d_germlines
    return [(a, s) for a, s in d_germlines if not a.startswith("TRBD2")]


#: Genomic 5'->3' rank of the D genes, for the loci whose architecture pins it independently of
#: species. TRB runs TRBD1 - J1 cluster - TRBD2 - J2 cluster; TRD runs TRDD1 - TRDD2 - TRDD3.
#: Both hold in human, mouse, rat and rhesus -- the same architecture argument `_allowed_d` makes.
#:
#: ⛔ IGH IS DELIBERATELY ABSENT. In *human* IMGT the second number of `IGHD<family>-<position>`
#: is the genomic position (IGHD1-1 .. IGHD7-27), but in *mouse* it is a family-member index with
#: no locus meaning -- and the two vocabularies collide on real gene names (`IGHD1-1`, `IGHD2-15`,
#: `IGHD5-5`, `IGHD5-12`, `IGHD6-6` exist in both). `_map_d` is handed sequences, not an organism,
#: so a name-parsed IGH rank would silently mis-order mouse. IGH tandem D-D is already down from
#: 11 calls to 2 on real bulk IGH after the `/OR` orphon exclusion above; the remaining 2 are not
#: worth a species-dependent rule. Genes absent here impose NO constraint.
_D_GENOMIC_ORDER = {
    "TRBD1": 0, "TRBD2": 1,
    "TRDD1": 0, "TRDD2": 1, "TRDD3": 2,
}


def _dd_orientation_ok(alleles5, alleles3) -> bool:
    """Can a deletional D-D join put ``alleles5`` 5' of ``alleles3`` on the read?

    D-D fusion is a rearrangement like any other: the upstream D joins to the downstream D and
    everything between them is deleted, so the fused product carries them in GENOMIC order. A
    read whose 5' D lies 3' of its 3' D in the germline locus therefore names a product that
    deletional joining cannot make. Measured on a real TRB amplicon: **10 of 15 tandem calls
    were TRBD2 -> TRBD1**, i.e. the impossible direction, which is what a chance second hit
    looks like when the two germlines are 16 nt and 12 nt of shared-composition sequence.

    Strict ``<``: equal rank means the SAME gene twice, which needs two germline copies and is
    likewise not producible.

    Either side may be an allele ambiguity list (byte-identical germlines across genes). One
    producible assignment is enough -- the call does not claim which member it was. A gene
    outside :data:`_D_GENOMIC_ORDER` (all of IGH) constrains nothing and passes.
    """
    genes5 = {a.split("*")[0] for a in alleles5}
    genes3 = {a.split("*")[0] for a in alleles3}
    for g5 in genes5:
        for g3 in genes3:
            r5, r3 = _D_GENOMIC_ORDER.get(g5), _D_GENOMIC_ORDER.get(g3)
            if r5 is None or r3 is None or r5 < r3:
                return True
    return False


def _d_min_score(interior_len: int, db_nt: int, max_evalue: float | None = None,
                 seqtype: str = "nt") -> int:
    """Smallest net score whose E-value clears ``max_evalue`` for this interior."""
    lam, default_e, floor = _d_gate(seqtype)   # read at call time so the gate stays tunable
    if max_evalue is None:
        max_evalue = default_e
    if interior_len <= 0 or db_nt <= 0 or max_evalue <= 0:
        return 1 << 30
    raw = math.log(_D_K * interior_len * db_nt / max_evalue) / lam
    return max(floor, math.ceil(raw))


def _d_evalue(score: int, interior_len: int, db_nt: int, seqtype: str = "nt") -> float:
    """Expected number of chance matches scoring at least ``score``."""
    lam, _, _ = _d_gate(seqtype)
    if score <= 0 or interior_len <= 0 or db_nt <= 0:
        return float("inf")
    return _D_K * interior_len * db_nt * math.exp(-lam * score)

# A tandem D-D fusion is sought in EVERY locus that has D germlines -- IGH, TRB and TRD.
# There used to be a ``_DD_LOCI = {"IGH", "TRD"}`` here, which made a TRBD1->TRBD2 tandem
# unrepresentable. TRB D-D fusions are real, and excluding TRB was doubly costly: the
# PRJNA371303 amplicons are TRA + TRB, and TRA has no D, so TRB is the *only* locus in the
# benchmark with a matched amplicon D truth. The one locus we could measure was the one we
# did not search. The rate of the second call (real or spurious) is what the shuffle
# control above measures; do not re-introduce a hard-coded locus allow-list to suppress it.


#: Prototype for :func:`_empty_record`. Copying a built dict is **8.8x** faster than rebuilding it
#: from a comprehension (0.127 s -> 0.014 s per 54,178 records, ~2 % of an amplicon run's wall):
#: `dict.copy` presizes and memcpys the table instead of hashing 52 keys and growing it. Every
#: value is an immutable `""`, so the copy is safe to hand out.
_EMPTY_RECORD_TEMPLATE = dict.fromkeys(AIRR_COLUMNS, "")


def _empty_record(query_id: str, query_seq: str) -> dict:
    rec = _EMPTY_RECORD_TEMPLATE.copy()
    rec["sequence_id"] = query_id
    rec["sequence"] = query_seq
    return rec


def _aln_identity_py(qaln: str, taln: str, tstart: int, t_lo: int, t_hi: int):
    """Fractional identity over the germline positions in target range [t_lo, t_hi].

    Walks the aligned strings tracking the target (scaffold) position; counts each
    target-consuming column in range (a germline position covered), and how many are an
    exact base match. Returns matches / covered as a 0-1 fraction, or "" if none covered.
    """
    t, covered, ident = tstart, 0, 0
    for qa, ta in zip(qaln, taln):
        if ta != "-":                                # target-consuming column
            if t_lo <= t <= t_hi:
                covered += 1
                if qa != "-" and qa.upper() == ta.upper():
                    ident += 1
            t += 1
    return round(ident / covered, 4) if covered else ""


try:
    from .._markup import aln_identity as _aln_identity_cpp
except ImportError:  # pragma: no cover - source checkout without the built extension
    _aln_identity_cpp = None


def _aln_identity(qaln: str, taln: str, tstart: int, t_lo: int, t_hi: int):
    """Fractional identity over the germline positions in target range [t_lo, t_hi].

    Per-column loop, so it lives in the extension; `_aln_identity_py` is the reference it is
    tested against and the fallback when the extension is not built. The C++ side returns -1.0
    for "nothing covered" because 0.0 is a real identity of zero, a different statement.
    """
    if _aln_identity_cpp is not None:
        v = _aln_identity_cpp(qaln, taln, tstart, t_lo, t_hi)
        return round(v, 4) if v >= 0.0 else ""
    return _aln_identity_py(qaln, taln, tstart, t_lo, t_hi)


def _project_point(hit: dict, ref_pos: int) -> int:
    """Project a single scaffold (target) nt position onto query coords (or 0)."""
    if ref_pos <= 0:
        return 0
    qs = _markup.transfer_regions(
        hit["qaln"], hit["taln"], int(hit["qstart"]), int(hit["tstart"]),
        [ref_pos], [ref_pos])[0][0]
    return qs if qs > 0 else 0


def _junction_nt(query_seq, cs, f4, coding_start, v_end_q):
    """Build (junction_nt, junction_aa, cdr3_aa, phase) for nucleotide input.

    ``cs`` = CDR3 start, ``f4`` = FR4 start (query, 1-based). The nucleotide
    junction is the real query slice (no synthetic bases). For translation, if the
    rearrangement is out of frame, 1-2 N are inserted after the V germline end so
    the J side reads correctly; the codon containing inserted N becomes ``_``.
    """
    js, je = cs - 3, f4 + 2            # Cys104 codon start .. [FW]118 codon end
    if js < 1 or je > len(query_seq):
        return "", "", "", None
    junction_nt = query_seq[js - 1 : je]
    phase = (f4 - coding_start) % 3    # 0 => J in V frame (in-frame rearrangement)
    k = (3 - phase) % 3
    if k:
        # Insert N after the V germline end, but keep the bridge strictly inside
        # CDR3 so the conserved Cys (first codon) and [FW] (last codon) are
        # preserved; fall back to just before the [FW] codon if V end is unknown.
        ins_at = (v_end_q - js + 1) if v_end_q else (len(junction_nt) - 3)
        lo, hi = 3, len(junction_nt) - 3            # after Cys codon .. before [FW] codon
        ins_at = max(0, hi) if hi < lo else min(max(ins_at, lo), hi)
        corrected = junction_nt[:ins_at] + "N" * k + junction_nt[ins_at:]
        junction_aa = list(translate(corrected, 0))
        for ci in range(ins_at // 3, (ins_at + k - 1) // 3 + 1):   # codons holding N
            if ci < len(junction_aa):
                junction_aa[ci] = "_"
        junction_aa = "".join(junction_aa)
    else:
        junction_aa = translate(junction_nt, 0)
    cdr3_aa = junction_aa[1:-1]
    return junction_nt, junction_aa, cdr3_aa, phase


def _best_d(interior, d_germlines, min_score, exclude=None):
    """Best-aligning D against ``interior`` as ``(score, length, alleles, s, e)``.

    ``s``/``e`` are 0-based inclusive offsets within ``interior``. ``exclude`` is an
    optional ``(s, e)`` span the match must not overlap (used to find a second,
    non-overlapping D). Returns ``None`` if nothing scores at least ``min_score``.

    ``alleles`` is the sorted tuple of EVERY allele producing the winning alignment --
    same score, same span -- not one arbitrary member of it. Seven pairs of human IGH D
    germlines are byte-identical across *different genes* (``IGHD4-11*01``/``IGHD4-4*01``,
    ``IGHD5-18*01``/``IGHD5-5*01``, and the five ``IGHD*/OR15-*a``/``*b`` pairs). A read
    matching one matches the other, at every score, always. The previous tiebreak picked
    the lexicographically smaller allele NAME, so ``IGHD4-11*01`` beat ``IGHD4-4*01``
    because ``"1" < "4"`` in ASCII -- an arbitrary gene call, made 100 % of the time for
    those pairs rather than occasionally, and a mechanical contributor to the 46-69 % IGH
    D concordance. Report the ambiguity, as ``v_call``/``j_call`` already do.

    Alleles that tie on ``(score, length)`` but align at a DIFFERENT span are not
    coordinate-compatible and cannot share one ``d_sequence_start``/``_end``; the winning
    span is chosen deterministically (highest score, then longest, then 5'-most) and only
    the alleles sharing *that* span are reported.
    """
    best = None                              # (rank, [alleles], s, e, ds, de)
    for allele, dseq in d_germlines:
        score, s, e, ds, de = _markup.d_local_align(interior, dseq)
        if score < min_score or s < 0:
            continue
        if exclude is not None:
            xs, xe = exclude
            if not (e < xs or s > xe):       # overlaps the excluded span
                continue
        # Rank only on the query-side ALIGNMENT, never on the allele name: equal rank means
        # an identical query span, hence genuine ambiguity rather than a winner to be broken.
        # The D-germline offsets (ds, de) are the rank-winner's; they belong to a single
        # allele, so _map_d only emits d_germline_*/d_cigar when the call is unambiguous.
        rank = (score, e - s + 1, -s)
        if best is None or rank > best[0]:
            best = (rank, [allele], s, e, ds, de)
        elif rank == best[0]:
            best[1].append(allele)
    if best is None:
        return None
    (score, length, _), alleles, s, e, ds, de = best
    # `set`: an aa reference lists one entry per (allele, reading frame), so two frames of the
    # same allele can tie at one span. That is one allele, not an ambiguity between two.
    return score, length, tuple(sorted(set(alleles))), s, e, ds, de


def _common_prefix_py(a: str, b: str) -> int:
    n = 0
    while n < len(a) and n < len(b) and a[n] == b[n]:
        n += 1
    return n


def _common_suffix_py(a: str, b: str) -> int:
    n = 0
    while n < len(a) and n < len(b) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


# Per-character Python loops called 138,065 times per 100k-read amplicon run (once per V allele in
# `v_anchor_prefix`, twice per read in `_anchored_vj_bounds`). The `_py` versions above stay as the
# reference implementation and the fallback; `tests/unit/` asserts the two agree.
_common_prefix = getattr(_markup, "common_prefix", None) or _common_prefix_py
_common_suffix = getattr(_markup, "common_suffix", None) or _common_suffix_py


#: nt of the Cys104 codon that a called V's own junction germline must explain before a junction
#: is emitted. 2, not 3: a synonymous TGT/TGC substitution is the one SHM event that hits the
#: conserved codon, and it leaves the first two bases intact.
MIN_V_ANCHOR_PREFIX = 2


def v_anchor_prefix(junction: str, v_call: str, anchors) -> int:
    """Longest prefix of ``junction`` explained by any called V's own junction germline.

    A junction's first codon *is* the V's conserved Cys104, so this is normally the full
    templated stretch. It is 0 when the rearrangement trimmed V back past that codon: the
    projection still lands somewhere and emits a junction opening on bases the V germline never
    templated. Measured on a TRA amplicon (results/round18): 1,396 of 46,785 reads disagree with
    IgBLAST, **every one of them a pure 5' over-extension** with the 3' end correct, and
    ``prefix < 2`` separates 1,360 of them from 44,322 correct junctions.
    """
    best = 0
    for allele in (v_call or "").split(","):
        a = anchors.get(("V", allele.strip())) if anchors else None
        g = getattr(a, "germline_nt", "") if a else ""
        if a and a.status == "ok" and g:
            best = max(best, _common_prefix(junction, g))
    return best


def _anchored_vj_bounds(query_seq, cs, f4, v_call, j_call, anchors, seqtype="nt"):
    """``(v_end_q, j_start_q)`` from the germline junction anchors, or ``(0, 0)``.

    The scaffold projection cannot locate these. A scaffold is ``V + 9 nt N-pad + J``,
    so the *scaffold* has only 9 nt between the V germline end and the J germline start.
    A real read has a 20-40 nt N-D-N region there, and mmseqs -- unable to align anything
    to a run of N -- parks those bases against the flanking V and J instead. Both
    projected boundaries then march inward: on the committed human fixtures the projected
    ``v_sequence_end`` sits +6 nt too far and ``j_sequence_start`` -13 nt too early
    (IGH, medians), leaving an 11 nt "interior" where the truth is 37 nt. D mapping was
    being handed a window too small to hold an IGH D segment at all.

    The per-allele germlines in ``cdr3_anchors.tsv`` give the answer directly: the V
    templates the junction's 5' end and the J its 3' end, so the boundaries are the
    longest common prefix / suffix. Against IgBLAST on the same fixtures this lands
    within 2 nt for 85 % of IGH and 93 % of TRB ``v_sequence_end`` (projection: 43 %,
    62 %), and 99 % of TRB ``j_sequence_start`` (projection: 49 %).

    Somatic hypermutation truncates the exact match early, which *widens* the interior --
    the safe direction: a D is never clipped, only surrounded by a little more sequence.

    For ``seqtype="aa"`` the same argument runs one residue at a time against the anchors'
    ``templated_aa`` instead of ``germline_nt``, and the junction runs Cys104..[FW]118 in
    residues rather than codons.
    """
    if not anchors or cs is None or f4 is None:
        return 0, 0, 0
    if seqtype == "aa":
        js0, je, field = cs - 1, f4, "templated_aa"   # junction = Cys104 .. [FW]118 residues
    else:
        js0, je, field = cs - 3, f4 + 2, "germline_nt"  # .. the same, as codons
    if js0 < 1 or je > len(query_seq):
        return 0, 0, 0
    junction = query_seq[js0 - 1 : je]

    def best(calls, segment, fn):
        # An ambiguous call lists several alleles; their junction germlines can differ
        # (TRBV20-1*01 templates CSAR, *03 only CSA). Take whichever explains the most.
        best_n = 0
        for allele in (calls or "").split(","):
            a = anchors.get((segment, allele.strip()))
            g = getattr(a, field, "") if a else ""
            if a and a.status == "ok" and g:
                best_n = max(best_n, fn(junction, g))
        return best_n

    p = best(v_call, "V", _common_prefix)
    s = best(j_call, "J", _common_suffix)
    # `p` IS `v_anchor_prefix(junction, v_call, anchors)` -- the same scan over the same slice
    # (`_junction_nt` cuts `query_seq[cs-3-1 : f4+2]`, identical to `junction` above). Returning it
    # saves the caller a second pass over every called V allele per read.
    if not p or not s:
        return 0, 0, p
    return js0 + p - 1, je - s + 1, p


def _map_d(rec, query_seq, v_end_q, j_start_q, d_germlines, j_call: str = "",
           seqtype: str = "nt", d_max_evalue: float | None = None):
    """Map D segment(s) into the V..J interior and populate d_call/np regions.

    Coordinates emitted are AIRR (1-based closed, query space). A second, non-overlapping
    D is sought in every locus that has D germlines (IGH, TRB, TRD -- the caller passes
    ``None`` for VJ loci, which returns immediately). The two are ordered 5'->3' as
    ``d_call`` / ``d2_call`` -- POSITIONALLY, so ``d_call`` is the 5' segment and need not
    be the higher-scoring one -- with ``np1``/``np2``/``np3`` between V, the D(s), and J.

    ``j_call`` restricts the candidate D set to those genomically 5' of the J (see
    ``_allowed_d``); pass it whenever the J is known. A tandem pair is additionally required
    to run 5'->3' in genomic order (see ``_dd_orientation_ok``); a pair that does not is
    reported as the single, higher-scoring D.

    ``d_max_evalue`` overrides the shipped operating point (``_D_MAX_EVALUE`` = 0.2 for nt,
    0.05 for aa) for BOTH segments. Lower is stricter: measured against IgBLAST at gene level,
    ``0.01`` agrees .9985 on a TRB amplicon and 1.0000 on bulk IGH, against .9765 / .9417 for
    the shipped 0.2 -- at a lower call rate. ``None`` keeps the shipped value.

    ``d_call`` and ``d2_call`` are comma-separated allele ambiguity lists (see ``_best_d``).

    For ``seqtype="aa"`` the interior, the D set (three translated frames per allele) and the
    emitted coordinates are all residues, following this module's convention that an aa record
    reuses the nt column names in aa space. ``d_germline_*``/``d_cigar`` are NOT emitted there:
    the alignment offsets index a translated reading frame, not the D germline, and reporting
    them as germline coordinates would be a lie.
    """
    if not d_germlines or not v_end_q or not j_start_q:
        return
    i_lo, i_hi = v_end_q + 1, j_start_q - 1       # 1-based interior bounds (query)
    if i_hi < i_lo:
        return
    interior = query_seq[i_lo - 1 : i_hi]
    # `n` stays the full locus D set: _D_MAX_EVALUE is an operating point calibrated against
    # it, so rescoring against the smaller masked set would silently loosen the gate. Held
    # fixed, masking can only remove an impossible call, never admit a weak one -- measured
    # on the realworld TRB fixtures, rescoring instead admitted 58 new calls at median
    # E = 0.098, squarely in the chance band.
    db_nt = _d_db_nt(d_germlines)
    d_germlines = _allowed_d(d_germlines, j_call)
    if not d_germlines:
        return
    min_score = _d_min_score(len(interior), db_nt, max_evalue=d_max_evalue, seqtype=seqtype)
    d1 = _best_d(interior, d_germlines, min_score)
    if d1 is None:
        return

    segs = [d1]
    # The same E-value gate serves the second segment: on true single-D junctions a
    # non-overlapping second D at this stringency is a false positive 0-1 % of the time.
    d2 = _best_d(interior, d_germlines, min_score, exclude=(d1[3], d1[4]))
    if d2 is not None:
        segs.append(d2)
    segs.sort(key=lambda c: c[3])                 # order 5'->3' by interior start
    # ⛔ ORIENTATION, after sorting and never before: the rule is about the order the two D
    # segments occupy ON THE READ, not about which of them scored higher. A pair running against
    # genomic order is not a weaker tandem, it is not a tandem -- so drop the second call
    # entirely and report the higher-scoring segment alone (`d1`), rather than manufacturing a
    # producible partner by reaching further down the score list.
    if len(segs) == 2 and not _dd_orientation_ok(segs[0][2], segs[1][2]):
        segs = [d1]

    def q(off):                                   # interior 0-based offset -> query 1-based
        return i_lo + off

    qlen = len(query_seq)

    def germline(rec, pfx, alleles, ds, de, qs, qe):
        # D germline coords + AIRR cigar, only when the call is a single allele -- an ambiguity
        # list (byte-identical D genes) has no one germline to anchor to. Gapless (one M run), with
        # the query 5' offset as leading S and the D-germline 5' offset as leading N (AIRR spec).
        if seqtype == "nt" and len(alleles) == 1 and ds >= 0:
            rec[f"{pfx}_germline_start"], rec[f"{pfx}_germline_end"] = ds + 1, de + 1
            rec[f"{pfx}_cigar"] = build_cigar(qs - 1, ds, ["M"] * (de - ds + 1), qlen - qe)

    sc1, _, a1, s1, e1, ds1, de1 = segs[0]
    rec["d_call"] = ",".join(a1)
    rec["d_sequence_start"], rec["d_sequence_end"] = q(s1), q(e1)
    rec["d_support"] = f"{_d_evalue(sc1, len(interior), db_nt, seqtype):.3g}"
    germline(rec, "d", a1, ds1, de1, q(s1), q(e1))
    if len(segs) == 2:
        sc2, _, a2, s2, e2, ds2, de2 = segs[1]
        rec["d2_call"] = ",".join(a2)
        rec["d2_sequence_start"], rec["d2_sequence_end"] = q(s2), q(e2)
        rec["d2_support"] = f"{_d_evalue(sc2, len(interior), db_nt, seqtype):.3g}"
        germline(rec, "d2", a2, ds2, de2, q(s2), q(e2))
        rec["np1"] = query_seq[v_end_q : q(s1) - 1]
        rec["np2"] = query_seq[q(e1) : q(s2) - 1]
        rec["np3"] = query_seq[q(e2) : j_start_q - 1]
    else:
        rec["np1"] = query_seq[v_end_q : q(s1) - 1]
        rec["np2"] = query_seq[q(e1) : j_start_q - 1]


def transfer_hit(
    query_id: str,
    query_seq: str,
    hit: dict,
    ref: RefEntry,
    seqtype: str = "nt",
    rev_comp: bool = False,
    d_germlines: list[tuple[str, str]] | None = None,
    submitted_seq: str | None = None,
    anchors: dict | None = None,
    d_max_evalue: float | None = None,
) -> dict:
    """Build an AIRR record by projecting ``ref`` region coords onto the query.

    ``query_seq`` is the coding-strand sequence all markup/coords/CIGARs are computed on.
    ``submitted_seq`` is the read AS SUBMITTED, stored verbatim in the AIRR ``sequence`` field;
    for a reverse-strand hit it is the reverse complement of ``query_seq`` and ``rev_comp`` is set,
    per AIRR ("if rev_comp is True, all output data are based on the reverse complement of
    ``sequence``"). Defaults to ``query_seq`` (forward reads, where the two are identical).

    ``d_max_evalue`` overrides the D-call E-value gate; see :func:`_map_d`.
    """
    # ⛔ ONE walk of the alignment, not four. Besides the seven regions, this function needs three
    # single scaffold positions projected onto the query: the V germline end, the J germline start,
    # and the V coding-frame anchor. Each used to go through `_project_point`, i.e. its own
    # `transfer_regions` crossing -- a fresh 6-argument binding call, two fresh `std::string` copies
    # of the SAME alignment, and a fresh forward walk, measured at 443 ns each against 822 ns for
    # the real multi-region call. Projecting a point is exactly the degenerate region [p, p], so
    # they ride along as extra intervals and are read back by index.
    #
    # `_project_point` is kept: `_extra_points` can only fold in positions known BEFORE the walk,
    # and the aa path still projects one that is not.
    ts = int(hit["tstart"])
    t0 = ref.starts[0]
    # The V coding-frame anchor: first scaffold position >= tstart that is in the V reading frame.
    # Depends only on `t0` and `tstart`, so it is knowable here. -1 (no V, e.g. a J+C scaffold)
    # projects to 0, which is exactly what `_project_point` returned for it.
    frame_pos = ts + ((t0 - ts) % 3) if t0 > 0 else 0
    extra = [ref.v_sequence_end or 0, ref.j_sequence_start or 0, frame_pos]
    n_reg = len(ref.starts)
    all_coords = _markup.transfer_regions(
        hit["qaln"], hit["taln"], int(hit["qstart"]), ts,
        list(ref.starts) + extra, list(ref.ends) + extra)
    coords = all_coords[:n_reg]
    # `transfer_regions` returns (-1, -1) for an uncovered position; `_project_point` returned 0.
    v_end_pt, j_start_pt, frame_pt = (max(0, all_coords[n_reg + i][0]) for i in range(3))

    rec = _empty_record(query_id, query_seq if submitted_seq is None else submitted_seq)
    rec.update(locus=ref.locus, v_call=ref.v_call, j_call=ref.j_call,
               c_call=ref.c_call, c_class=isotype_class(ref.c_call),
               rev_comp="T" if rev_comp else "F", productive="")

    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return ""

    # Alignment quality of the scaffold hit (mmseqs2 bit score / E-value / % identity).
    rec["mmseqs2_score"] = _num(hit.get("bits"))
    rec["mmseqs2_evalue"] = _num(hit.get("evalue"))
    rec["mmseqs2_identity"] = _num(hit.get("pident"))
    for _c in ("qstart", "qend", "qlen", "tstart", "tend", "tlen"):
        rec[f"mmseqs2_{_c}"] = _num(hit.get(_c))
    rec["mmseqs2_t_vend"] = ref.v_sequence_end or ""
    rec["mmseqs2_t_jstart"] = ref.j_sequence_start or ""
    rec["mmseqs2_t_vjend"] = ref.vj_end or ""

    # AIRR alignment strings + germline coordinates. The scaffold is `V + N*pad + J [+ C]`, so the
    # target (scaffold) span maps to germline coords directly: V germline pos == target pos; J
    # germline pos == target pos - t_jstart + 1. C is not part of a germline V/J alignment.
    rec["sequence_alignment"] = hit.get("qaln") or ""
    rec["germline_alignment"] = hit.get("taln") or ""
    ts, te = int(hit["tstart"]), int(hit["tend"])
    t_vend = int(ref.v_sequence_end) if ref.v_sequence_end else 0
    t_jstart = int(ref.j_sequence_start) if ref.j_sequence_start else 0
    t_vjend = int(ref.vj_end) if ref.vj_end else 0
    if t_vend and ts <= t_vend:                      # alignment covers some V germline
        rec["v_germline_start"], rec["v_germline_end"] = ts, min(te, t_vend)
        rec["v_identity"] = _aln_identity(hit["qaln"], hit["taln"], ts, ts, t_vend)
    if t_jstart and t_vjend and te >= t_jstart:      # alignment reaches the J germline
        rec["j_germline_start"] = max(ts, t_jstart) - t_jstart + 1
        rec["j_germline_end"] = min(te, t_vjend) - t_jstart + 1
    # Per-segment CIGARs (v/j/c) in a single walk of the same aligned strings.
    rec.update(segment_cigars(hit["qaln"], hit["taln"], int(hit["qstart"]), ts,
                              len(query_seq), t_vend, t_jstart, t_vjend))
    # The junction boundary in GERMLINE coordinates -- see AIRR_COLUMNS for why it is emitted and
    # why it is raw. Taken from the FIRST call of a tie list, which is the one the positions in
    # `v_mutations` / `j_mutations` are expressed against.
    if anchors:
        for col, seg, call in (("v_anchor_nt", "V", ref.v_call), ("j_anchor_nt", "J", ref.j_call)):
            allele = (call or "").split(",")[0].strip()
            a = anchors.get((seg, allele)) if allele else None
            if a is not None and a.status == "ok":
                rec[col] = a.anchor_nt

    v_prefix: int | None = None      # set by `_anchored_vj_bounds`; reused by the junction gate
    region_q: dict[str, tuple[int, int]] = {}
    for name, (qs, qe) in zip(REGIONS, coords):
        if qs < 0:
            continue
        region_q[name] = (qs, qe)
        rec[f"{name}_start"], rec[f"{name}_end"], rec[name] = qs, qe, query_seq[qs - 1 : qe]

    # CDR3 end is J-anchored (somatic length is query-specific): from the V-anchored
    # start up to just before the [FW] that opens FR4.
    if "cdr3" in region_q and "fwr4" in region_q:
        cs, ce = region_q["cdr3"][0], region_q["fwr4"][0] - 1
        if ce >= cs:
            region_q["cdr3"] = (cs, ce)
            rec["cdr3_start"], rec["cdr3_end"], rec["cdr3"] = cs, ce, query_seq[cs - 1 : ce]

    # Transfer the V germline end and J germline start (extended scaffold markup), then
    # refine them against the per-allele junction germlines -- the projection systematically
    # collapses the V..J interior (see `_anchored_vj_bounds`).
    v_end_q, j_start_q = v_end_pt, j_start_pt      # from the single walk above
    if "cdr3" in region_q and "fwr4" in region_q:
        av, aj, v_prefix = _anchored_vj_bounds(
            query_seq, region_q["cdr3"][0], region_q["fwr4"][0],
            ref.v_call, ref.j_call, anchors, seqtype)
        if av and aj and aj > av:
            v_end_q, j_start_q = av, aj
    if "fwr1" in region_q:
        rec["v_sequence_start"] = region_q["fwr1"][0]
    if v_end_q:
        rec["v_sequence_end"] = v_end_q
    if j_start_q:
        rec["j_sequence_start"] = j_start_q
    if t_jstart and t_vjend and not rec["j_germline_start"] and not j_start_q:
        # Neither the alignment nor the junction anchor found any J in this read, so `ref.j_call`
        # names the J of whichever V*J scaffold the read landed on -- not a J the read carries
        # evidence for. On bulk RNA-seq that is 1,823 of 2,737 mapped reads (results/round18 §4).
        # Gated on the scaffold declaring where its J starts: a reference that cannot say (the aa
        # markup does not carry `j_sequence_start`) makes an absent `j_germline_start` mean
        # nothing, and blanking on it deleted every protein-input `j_call`.
        rec["j_call"] = ""

    if seqtype == "nt":
        # V coding frame from the alignment phase (works even without FR1). A `J + C` scaffold has no
        # V, so ``ref.starts[0]`` is -1 and there is no V frame to project -- guard the arithmetic
        # rather than feed -1 into it. FR4 still reads in its own (J) frame, so it is translated below,
        # outside this branch: on a V-less hit it is the only markup there is.
        coding_start = None
        if t0 > 0:
            pj = frame_pt                          # from the single walk above
            coding_start = pj or region_q.get("fwr1", (None,))[0]
            if pj == 0 and coding_start is not None:
                v_end = region_q.get("fwr3", (0, 0))[1] or len(query_seq)
                coding_start += detect_coding_frame(query_seq[coding_start - 1 : v_end])
        if coding_start is not None:
            protein = translate(query_seq[coding_start - 1:], 0)
            for name, (qs, qe) in region_q.items():       # V-side aa from V frame
                if name == "fwr4":
                    continue                              # J frame, not V frame -- see below
                a_s, a_e = aa_coords_from_nt(qs, qe, coding_start)
                rec[f"{name}_aa"] = protein[max(1, a_s) - 1 : a_e]
        # FR4 reads in its own (J) frame regardless of productivity, and regardless of whether a V
        # frame exists at all.
        if "fwr4" in region_q:
            f4s, f4e = region_q["fwr4"]
            rec["fwr4_aa"] = translate(query_seq[f4s - 1 : f4e], 0)
        if coding_start is not None:
            # Junction (+ CDR3 aa) with out-of-frame N-bridging. Needs the V's conserved Cys104.
            phase = None
            if "cdr3" in region_q and "fwr4" in region_q:
                jnt, jaa, c3aa, phase = _junction_nt(
                    query_seq, region_q["cdr3"][0], region_q["fwr4"][0],
                    coding_start, v_end_q)
                # Refuse a junction whose opening codon the called V's germline does not template:
                # the read's V was trimmed past Cys104, so this window is 5'-over-extended.
                # Declining is what IgBLAST and MiXCR both do here (see `v_anchor_prefix`).
                # `_anchored_vj_bounds` already scanned this exact slice, so reuse its answer
                # rather than walking every called V allele a second time.
                prefix = (v_prefix if v_prefix is not None
                          else v_anchor_prefix(jnt, ref.v_call, anchors))
                if jaa and anchors and prefix < MIN_V_ANCHOR_PREFIX:
                    jaa, phase = "", None
                if jaa:
                    rec["junction"], rec["junction_aa"], rec["cdr3_aa"] = jnt, jaa, c3aa
            vclean = all("*" not in rec.get(f"{r}_aa", "") for r in _VSIDE)
            jclean = "*" not in rec.get("junction_aa", "") and "_" not in rec.get("junction_aa", "")
            # `phase is None` means the read never reached BOTH cdr3 and fwr4, so no junction was
            # observed -- and productivity is a property of the V-J junction. Such a read is
            # `unevaluable`, not `non-productive`, exactly as a V-less read is (see the `else`
            # below, which has always got this right). Reporting "F" made a bare V fragment look
            # like a confirmed non-productive rearrangement: on real bulk output that is 72 % of
            # mapped reads, since most reads land wholly inside V and never reach CDR3. AIRR's
            # `productive` is "predicted to be productive" -- with no junction there is nothing to
            # predict from, and empty is the honest answer.
            if phase is None:
                rec["productive"] = ""
                rec["vj_in_frame"] = ""
            else:
                rec["productive"] = "T" if (phase == 0 and vclean and jclean) else "F"
                rec["vj_in_frame"] = "T" if phase == 0 else "F"
            # `stop_codon` is NOT gated the same way: a stop in the V-side regions is directly
            # observed whether or not the junction was reached, so it stays evaluable here.
            j_stop = "*" in rec.get("junction_aa", "")
            rec["stop_codon"] = "F" if (vclean and not j_stop) else "T"
        # else: `productive` stays "" -- a V-less read is not "non-productive", it is unevaluable.
        # D-segment mapping (VDJ loci only; gated by presence of D germlines).
        _map_d(rec, query_seq, v_end_q, j_start_q, d_germlines, ref.j_call,
               d_max_evalue=d_max_evalue)
    else:  # aa input: regions are already amino acids; no frame bridging needed.
        for name, (qs, qe) in region_q.items():
            rec[f"{name}_aa"] = query_seq[qs - 1 : qe]
        if "cdr3" in region_q and "fwr4" in region_q:
            cs, f4 = region_q["cdr3"][0], region_q["fwr4"][0]
            if cs >= 2 and f4 <= len(query_seq):
                rec["junction"] = query_seq[cs - 2 : f4]
                rec["junction_aa"] = rec["junction"]
                rec["cdr3_aa"] = rec["junction"][1:-1]
        if "fwr1" in region_q and "fwr4" in region_q:
            span = query_seq[region_q["fwr1"][0] - 1 : region_q["fwr4"][1]]
            rec["productive"] = "T" if "*" not in span else "F"
        # D (and tandem D-D) against the three translated frames of each D germline. Only IGH
        # keeps enough surviving D to see in protein; the TR loci mostly stay silent, which is
        # the honest answer. Same call, same columns, aa coordinates.
        _map_d(rec, query_seq, v_end_q, j_start_q, d_germlines, ref.j_call, seqtype="aa",
               d_max_evalue=d_max_evalue)
    return rec
