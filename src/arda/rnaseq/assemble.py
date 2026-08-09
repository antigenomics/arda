"""Stage 3 — contig assembly (anchored greedy overlap-extension).

Reconstruct the clonotypes that Stage 1 maps but cannot *call*: a long CDR3 (V(DD)J
ultralong, ~20-40 aa) does not fit in one 100-150 bp read, so no read spans the junction
and :mod:`~arda.rnaseq.correct`'s complete-junction filter drops every read of it. Assembly-based
extractors recover these by assembly; this module does the same on the reads Stage 1 already mapped.

Why overlap-extension and not a de Bruijn graph: every clonotype sharing a germline V/J
contributes identical k-mers, so a dBG collapses distinct clones exactly across the CDR3
(the region of interest). This is the reason the Pevzner-lab Ig assemblers use a read graph,
not a k-mer graph (Safonova 2015, 10.1093/bioinformatics/btv238). We exploit arda's own
anchors instead: Stage 1 gives every V-side read a ``cdr3_start`` offset, so reads of one
clone are already coordinate-aligned at the CDR3 -- we seed from those and extend 3' through
the CDR3 into J, where the sequence is clone-specific, and stop before running deep into the
(shared) constant region. Seeding never extends 5' into the germline V, which is what keeps
distinct clones apart and bounds the germline-k-mer blow-up.

The assembled contig physically unites a clone's V-side reads with its J/C-side reads under
one junction, so :mod:`~arda.rnaseq.correct` also gets the clone's isotype (from the J/C
mates' ``c_class``) for free -- the long clones were previously invisible to both.

Output is a per-member-read AIRR fragment (``sequence_id`` -> the contig's complete junction),
meant to be concatenated with the Stage-1 mapped AIRR and fed to ``correct`` in one pass: a
read that was incomplete in Stage 1 is dropped there and kept here, so fragments count once.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ..annotate.airr_out import read_airr as _read_airr

from ._res import Stage
from ..annotate.contig import reannotate_contigs

__all__ = ["assemble_contigs", "AssembleReport"]



# D columns carried from the contig annotation onto each member read. Coordinates are omitted
# on purpose: they index the contig, and these rows carry only the junction.
_D_COLUMNS = ("d_call", "d2_call", "d_support", "d2_support", "np1", "np2", "np3")
_OUT_COLUMNS = ("sequence_id", "locus", "v_call", "j_call", "c_call", "c_class",
                "junction", "junction_aa", *_D_COLUMNS)

_COMP = str.maketrans("ACGTN", "TGCAN")
_LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
# A complete junction spans both conserved anchors, in frame, no stop -- same rule as
# ``correct._COMPLETE`` (kept local to avoid importing a private symbol).
_CANON = re.compile(r"^C[ACDEFGHIKLMNPQRSTVWY]*[FW]$")

#: How much of the contig's junction a member read must actually cover before the contig's junction
#: is attributed to it. Membership only requires a ``min_overlap`` match, and after the extension
#: passes accumulate germline at the contig ends that match can be entirely germline V or J -- so a
#: read of a DIFFERENT clone of the same V gene qualifies. Ten nucleotides is a little over three
#: codons of clone-specific sequence: enough that the read cannot be explained by germline alone,
#: and short enough that a read clipping the junction's edge still counts.
#: ⚠ Failing this does NOT lose the read -- it stays in the Stage-1 frame for coverage assignment.
#: It only withholds the FORCED attribution to this clonotype.
_MIN_JUNCTION_COVER = 10


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _gene(x: str | None) -> str:
    return (x or "").split(",")[0].split("(")[0].split("*")[0].strip()


def _locus_of(row: dict) -> str:
    loc = (row.get("locus") or "")[:3]
    if loc in _LOCI:
        return loc
    for c in ("v_call", "j_call", "c_call"):
        g = _gene(row.get(c))
        if g[:3] in _LOCI:
            return g[:3]
    return ""


def _complete(ja: str | None, jn: str | None) -> bool:
    return bool(ja) and bool(jn) and len(jn) % 3 == 0 and "*" not in ja and "_" not in ja \
        and _CANON.match(ja) is not None


@dataclass
class AssembleReport:
    reads_in: int = 0
    seeds: int = 0
    contigs: int = 0
    contigs_complete: int = 0          # contigs whose re-annotated junction is complete
    reads_rescued: int = 0             # incomplete member reads attributed to a complete contig
    #: Members skipped because their span did not cover the contig's junction by
    #: ``_MIN_JUNCTION_COVER`` -- recruited on germline-only overlap, so there is no evidence they
    #: belong to this clonotype. They keep their reads; only the attribution is withheld.
    members_without_junction: int = 0
    # See `_res.Stage`: peak is the WHOLE-PROCESS high-water mark as of this stage's end
    # (monotone -- getrusage offers no per-stage reset), gain is this stage's contribution.
    wall_seconds: float = 0.0
    peak_rss_mb: float = 0.0
    rss_gain_mb: float = 0.0

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _greedy_contigs(
    oriented: list[str],
    seed_idx: list[int],
    cdr3_start: list[int | None],
    *,
    k: int,
    min_overlap: int,
    min_id: float,
    max_ext_past_cdr3: int,
    scan_cap: int,
    min_v: int = 70,
) -> list[tuple[str, list[int], list[tuple[int, int]]]]:
    """Seed from ``seed_idx`` (longest CDR3 tail first), extend 3' into J then 5' into V.

    Returns ``[(contig_seq, member_read_indices, member_spans)]``, where each span is the member's
    ``[start, end)`` in FINAL contig coordinates. The spans exist so attribution can be checked:
    membership is granted on a germline-shared overlap, so a member need not have covered the
    junction at all -- see :func:`assemble_contigs`. ``used`` partitions reads: each read
    joins at most one contig (dominant clones seed first), so member counts don't overlap.
    The 5' pass extends ~``min_v`` nt into the V so re-annotation can call the V gene (mmseqs
    needs enough germline to anchor it); that region is shared germline, so any V-read of the
    gene extends it correctly -- no chimera risk, and ``scan_cap`` bounds the germline-k-mer cost.
    """
    # ⛔ Cap the postings AT INSERT, not only when reading them. Both consumers already take
    # `[:scan_cap]` of the posting list, so keeping only the first `scan_cap` entries yields the
    # IDENTICAL candidate set -- this is a pure memory fix, not a behaviour change. Unbounded, the
    # index held every k-mer position of every mapped read of the locus (~90 postings per 100 nt
    # read) in a Python dict of int lists, which is ~460 MB per 100k reads; `--assemble` is ON by
    # default, so every run paid it. `_assign_coverage` bounds its equivalent index the same way.
    index: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(oriented):
        for p in range(0, len(s) - k + 1):
            lst = index[s[p:p + k]]
            if len(lst) < scan_cap:
                lst.append(i)

    used = [False] * len(oriented)
    # longest CDR3 tail first: that read spans the most of the junction, the best seed.
    seed_idx = sorted(seed_idx, key=lambda i: len(oriented[i]) - (cdr3_start[i] or 0), reverse=True)
    contigs: list[tuple[str, list[int], list[tuple[int, int]]]] = []
    for si in seed_idx:
        if used[si]:
            continue
        contig = oriented[si]
        anchor = cdr3_start[si] or 0            # CDR3 start within the contig
        members = [si]
        # Each member's [start, end) in CONTIG coordinates. The 5' pass prepends, so every existing
        # span shifts right by the length of what was prepended.
        spans: list[tuple[int, int]] = [(0, len(contig))]
        used[si] = True
        # 3' pass: extend through the CDR3 into J (clone-specific until the germline J).
        while len(contig) - anchor < max_ext_past_cdr3:
            tail = contig[-k:]
            best_j, best_ext = None, ""
            for j in index.get(tail, ())[:scan_cap]:
                if used[j]:
                    continue
                s = oriented[j]
                p = s.find(tail)                # where the shared k-mer sits in read j
                if p < 0:
                    continue
                start = len(contig) - k - p     # where read j would begin in contig coords
                if start < 0:                   # read reaches 5' of the contig; not a 3' extension
                    continue
                ov = len(contig) - start        # overlap length
                if ov < min_overlap:
                    continue
                a, b = contig[start:], s[:ov]
                if sum(1 for x, y in zip(a, b) if x != y) > (1 - min_id) * ov:
                    continue
                ext = s[ov:]
                if len(ext) > len(best_ext):
                    best_j, best_ext = j, ext
            if best_j is None or not best_ext:
                break
            contig += best_ext
            used[best_j] = True
            members.append(best_j)
            spans.append((len(contig) - len(oriented[best_j]), len(contig)))
        # 5' pass: extend into the V until there is enough germline for re-annotation to call V.
        while anchor < min_v:
            head = contig[:k]
            best_j, best_ext = None, ""
            for j in index.get(head, ())[:scan_cap]:
                if used[j]:
                    continue
                s = oriented[j]
                p = s.find(head)                # head sits at position p in read j
                if p <= 0:                      # j must extend 5' of the contig
                    continue
                ov = min(len(s) - p, len(contig))
                if ov < min_overlap:
                    continue
                a, b = contig[:ov], s[p:p + ov]
                if sum(1 for x, y in zip(a, b) if x != y) > (1 - min_id) * ov:
                    continue
                ext = s[:p]
                if len(ext) > len(best_ext):
                    best_j, best_ext = j, ext
            if best_j is None or not best_ext:
                break
            contig = best_ext + contig
            anchor += len(best_ext)
            spans = [(a + len(best_ext), b + len(best_ext)) for a, b in spans]
            used[best_j] = True
            members.append(best_j)
            spans.append((0, len(oriented[best_j])))
        if len(members) >= 2:
            contigs.append((contig, members, spans))
        else:
            # ⛔ RELEASE a rejected contig's reads. `used` is set as reads are recruited, but a
            # contig dropped here never gave them back, so a seed that failed to extend was
            # permanently consumed -- it could no longer join a LATER seed's contig even as an
            # ordinary extension member. Seeds are tried longest-CDR3-tail first, so the reads this
            # stranded were exactly the short-tailed ones that most need a contig to reach a
            # junction. A rejected contig here always has just its seed (any successful extension
            # would have made it >= 2), so this releases one read.
            for mi in members:
                used[mi] = False
    return contigs


def assemble_contigs(
    airr_tsv: str | Path,
    output: str | Path,
    *,
    organism: str = "human",
    k: int = 21,
    min_overlap: int = 21,
    min_id: float = 0.90,
    max_ext_past_cdr3: int = 130,
    scan_cap: int = 400,
    threads: int = 0,
    map_d: bool = True,
    d_max_evalue: float | None = None,
    report_path: str | Path | None = None,
) -> AssembleReport:
    """Assemble long-CDR3 contigs from Stage-1 mapped reads and attribute their junctions.

    Reads the Stage-1 mapped-reads AIRR (needs ``sequence``, ``rev_comp``, ``locus``,
    ``cdr3_start`` and the ``v/j/c_call`` columns), assembles per-locus contigs, re-annotates
    them (:func:`~arda.annotate.contig.reannotate_contigs`), and writes an AIRR TSV with one
    row per **incomplete** member read carrying the contig's complete junction (and the read's
    own ``c_class`` so isotype survives). Concatenate this with the mapped AIRR and run
    ``correct`` once: the read's incomplete Stage-1 row is dropped and this complete row kept,
    so each fragment is counted exactly once.

    The contig's D call travels with the junction (``d_call``, ``d2_call``, ``d_support``,
    ``d2_support``, ``np1``-``np3``). An ultralong CDR3 is the one place a tandem D-D is both
    most likely and least visible to a single read, so the contig is where it must be called.

    Args:
        airr_tsv: Stage-1 mapped-reads AIRR TSV.
        output: assembled-reads AIRR TSV (header only if nothing assembles).
        map_d: map D segments on the assembled contig (default ``True``).
        d_max_evalue: E-value gate on the D call(s); ``None`` keeps the shipped 0.2.
        max_ext_past_cdr3: stop extending a contig once it reaches this many nt past the CDR3
            start -- enough to cross the junction into J without running into the shared C region.
        scan_cap: per-step cap on candidate reads examined for a (germline-frequent) k-mer.

    Returns:
        An :class:`AssembleReport`.
    """
    output = Path(output)
    raw = _read_airr(airr_tsv)
    stage = Stage()
    report = AssembleReport(reads_in=raw.height)
    n = raw.height
    # Only the columns this function reads. A Stage-1 AIRR has 83, and `to_list()` on the
    # rest builds Python str objects for `sequence_alignment`, `germline_alignment` and
    # every region sequence -- measured 2.42 KB/row against 0.41 KB/row for the columns
    # actually used, i.e. ~2 KB wasted per mapped read (~7 GB at SRR5233639's full depth).
    # `col()` below already tolerates an absent column, so restricting the dict is safe.
    _USED = ("c_call", "c_class", "cdr3_start", "j_call", "junction", "junction_aa",
             "locus", "rev_comp", "sequence", "sequence_id", "v_call")
    cols = {c: raw[c].to_list() for c in raw.columns if c in _USED}

    def col(name):
        return cols.get(name, [None] * n)

    seq_c, rc_c = col("sequence"), col("rev_comp")
    ja_c, jn_c, cs_c = col("junction_aa"), col("junction"), col("cdr3_start")
    sid_c, cclass_c = col("sequence_id"), col("c_class")
    loc_c, v_c, j_c, c_c = col("locus"), col("v_call"), col("j_call"), col("c_call")

    # Bucket read indices by locus; orient to the receptor-coding strand.
    by_locus: dict[str, list[int]] = defaultdict(list)
    oriented: list[str] = [""] * n
    cdr3_start: list[int | None] = [None] * n
    complete: list[bool] = [False] * n
    for i in range(n):
        loc = _locus_of({"locus": loc_c[i], "v_call": v_c[i], "j_call": j_c[i], "c_call": c_c[i]})
        if not loc:
            continue
        s = seq_c[i] or ""
        if str(rc_c[i]).upper() in ("T", "TRUE", "1"):
            s = _rc(s)
        oriented[i] = s
        cs = cs_c[i]
        cdr3_start[i] = int(cs) if cs not in (None, "", "None") else None
        complete[i] = _complete(ja_c[i], jn_c[i])
        by_locus[loc].append(i)

    # Assemble per locus. Seeds = incomplete reads that reach the CDR3 (the ones `correct` drops).
    contig_records: list[tuple[str, str]] = []          # (contig_id, seq) for reannotate
    contig_members: list[list[int]] = []
    contig_spans: list[list[tuple[int, int]]] = []      # each member's [start, end) in the contig
    for loc, idxs in by_locus.items():
        local_seq = [oriented[i] for i in idxs]
        local_cs = [cdr3_start[i] for i in idxs]
        seeds_local = [p for p, i in enumerate(idxs)
                       if cdr3_start[i] is not None and not complete[i] and len(oriented[i]) >= k]
        if not seeds_local:
            continue
        report.seeds += len(seeds_local)
        for cseq, members_local, spans_local in _greedy_contigs(
                local_seq, seeds_local, local_cs,
                k=k, min_overlap=min_overlap, min_id=min_id,
                max_ext_past_cdr3=max_ext_past_cdr3, scan_cap=scan_cap):
            gid = len(contig_records)
            contig_records.append((f"contig_{loc}_{gid}", cseq))
            contig_members.append([idxs[p] for p in members_local])
            contig_spans.append(spans_local)
    report.contigs = len(contig_records)
    if not contig_records:
        _write_empty(output)
        stage.finish(report)
        if report_path:
            Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
        return report

    ann = reannotate_contigs(contig_records, organism, threads=threads, map_d=map_d,
                             d_max_evalue=d_max_evalue)
    ann_by_id = {a.get("sequence_id"): a for a in ann}

    rows: list[dict] = []
    for (cid, cseq), members, spans in zip(contig_records, contig_members, contig_spans):
        a = ann_by_id.get(cid)
        if a is None or not _complete(a.get("junction_aa"), a.get("junction")):
            continue
        report.contigs_complete += 1
        jn, ja = a["junction"], a["junction_aa"]
        vc, jc = a.get("v_call") or "", a.get("j_call") or ""
        lc = (a.get("locus") or "")[:3] or _locus_of(a)
        # The contig is the only place a long CDR3's D is visible at all -- no single read spans
        # it. Carry the contig's D across to every member read. Calls, E-values and the np
        # stretches are junction-scoped facts and stay true on a row that carries only the
        # junction; d_sequence_start/end are CONTIG coordinates and would be meaningless here,
        # so they are not propagated.
        d = {c: a.get(c) or "" for c in _D_COLUMNS}
        # ⛔ ATTRIBUTION NEEDS EVIDENCE. Membership is granted on a >= `min_overlap` match, and after
        # the extension passes have accumulated germline at the contig ends that overlap can be pure
        # germline -- the 5' pass says so in its own docstring ("that region is shared germline, so
        # any V-read of the gene extends it correctly"). Stamping the contig's junction onto such a
        # member credits a read carrying ZERO clone-specific evidence to this clonotype's
        # duplicate_count, and a read of a different clone of the same V gene is exactly the kind of
        # thing that overlap admits. So a member is only attributed if it actually COVERED the
        # junction. Locate the junction in the contig and require a real overlap with it.
        #
        # ⚠ Dropping a rescued row does not lose the read: it stays in the Stage-1 frame and
        # `correct._assign_coverage` can still place it on evidence. What it loses is the FORCED
        # attribution.
        j0 = cseq.find(jn)
        j1 = j0 + len(jn) if j0 >= 0 else -1
        for mi, (ms, me) in zip(members, spans):
            if complete[mi]:          # already counted via the mapped AIRR; don't double-attribute
                continue
            if j0 >= 0 and min(me, j1) - max(ms, j0) < _MIN_JUNCTION_COVER:
                report.members_without_junction += 1
                continue
            rows.append({
                "sequence_id": sid_c[mi], "locus": lc, "v_call": vc, "j_call": jc,
                "c_call": "", "c_class": cclass_c[mi] or "", "junction": jn, "junction_aa": ja,
                **d,
            })
    report.reads_rescued = len(rows)

    if rows:
        pl.DataFrame(rows).select(_OUT_COLUMNS).write_csv(output, separator="\t")
    else:
        _write_empty(output)
    stage.finish(report)
    if report_path:
        Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
    return report


def _write_empty(output: Path) -> None:
    pl.DataFrame(schema={c: pl.Utf8 for c in _OUT_COLUMNS}).write_csv(output, separator="\t")
