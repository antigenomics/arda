"""Place the junction by arithmetic instead of by alignment.

AIRR ``junction`` runs from the Cys104 codon to the [FW]118 codon that opens FR4, both **included**
(IMGT ``cdr3`` excludes them and is two residues shorter). Today arda learns those two positions by
aligning the read against a ``V + pad + J`` scaffold, and the 15,414-scaffold reference exists for
that and essentially nothing else -- the V and J *gene* calls already come out of the 924-target
segment pass at .9997 / .9998 agreement.

But the segment pass already returns, per read per side, ``(target, tstart, qstart)``; and
``cdr3_anchors.tsv`` already records ``anchor_nt``, the 0-based offset of the anchor codon inside
that segment's germline. The segment targets are built so target coordinates *are* germline
coordinates -- ``refbuild.segments`` writes the V target as ``scaffold[:v_sequence_end]`` and the J
target as ``scaffold[j_sequence_start - 1 : vj_end]``, so each germline starts at offset 0. So the
anchor's position in the read is two integer adds:

.. code-block:: text

    offset = (anchor_nt + 1) - tstart          # germline distance from the hit start to the anchor
    pos    = qstart + offset                   # forward hit
    pos    = (qlen - qstart + 1) + offset      # reverse-complement hit, in revcomp coordinates

Measured against IgBLAST on a TRA amplicon, a TRB amplicon and a bulk library:
**54,740 / 54,756 = .99971 byte-exact junctions**, at no cost beyond the segment pass that already
runs. See ``arda-benchmark/results/round14/README.md`` §2.

⛔ **This is a fast path, not a replacement.** It yields ``junction`` and its coordinates. It does
*not* yield ``v_identity``, ``sequence_alignment``, ``germline_alignment``, the per-segment CIGARs
and mutation lists, or the ``mmseqs2_*`` block, all of which ``annotate.transfer`` derives from the
alignment's
``qaln``/``taln``. Anything needing those still needs the scaffold alignment.

⛔ **And it refuses rather than degrades.** A junction that is well-formed but wrong is the worst
output this codebase can produce -- the reference-geometry bug shipped junctions that started ``C``,
ended ``[FW]``, passed ``--complete-only`` and were short by exactly the allele's truncation. Every
condition this projection cannot verify sends the read back to the aligner instead of guessing. Same
rule as :func:`arda.igblast.auxiliary_data`, which raises rather than returning ``None``.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Projection", "project_junction", "REFUSALS"]

#: Every reason a read is handed back to the scaffold aligner. Counted in the run report so the
#: fast-path yield is auditable rather than inferred -- 87.0 % of hit TRA-amplicon reads carry both
#: anchors against 7.2 % of bulk reads, and a silent fast path would hide that difference entirely.
REFUSALS = ("no_anchor", "unvalidated_locus", "indel_unchecked", "indel_split",
            "strand_mismatch", "order", "off_read", "bad_codon")

#: Loci the projection declines regardless of how well the arithmetic works on them.
#:
#: ⛔ TRD is here because it is UNVALIDATED, not because it is known bad. The pre-registered bar for
#: shipping a locus was "byte-exact >= 0.99 at n >= 2,000, or the locus goes on the refusal list",
#: and TRD came back at **n = 0**: across two TR amplicons the segment pass never handed a single
#: TRD read both anchors, so all 767 TRD junctions in the IgBLAST truth fell through to the aligner
#: untouched. Zero coverage is not a pass. The one measurement that does exist is the design pass's
#: 43/51 = 0.843, against .999 pooled -- far too small to be a rate and far too poor to ignore.
#:
#: TRD is also the locus most likely to break a projection: TRAV/DV alleles rearrange to either TRAJ
#: or TRDJ and **the J decides the locus**, so a V-derived coordinate can be read in the wrong
#: frame of reference. Declining costs nothing today (the path never fires) and stops a future
#: reference change from silently routing TRD reads through untested arithmetic.
UNVALIDATED_LOCI = frozenset({"TRD"})


@dataclass(frozen=True)
class Projection:
    """A junction placed from segment coordinates alone."""

    junction: str
    #: 1-based, inclusive, on ``strand_seq`` -- i.e. on the reverse complement when ``rev_comp``.
    start: int
    end: int
    rev_comp: bool

    @property
    def cdr3(self) -> str:
        """IMGT ``cdr3``: the junction with both anchor codons removed."""
        return self.junction[3:-3]


def _anchor(anchors: dict, side: str, call: str):
    """Resolve a possibly comma-joined segment call to a single usable anchor.

    A segment target may name several alleles arda cannot separate (23 of 775 human ``V|`` targets
    do, ``IGHV3-23*01,IGHV3-23D*01`` among them). They are duplicate SEQUENCES, so their anchors
    agree by construction -- verified zero disagreement across every multi-allele header -- and
    taking the first resolvable one is exact, not a heuristic.
    """
    for allele in call.split(","):
        a = anchors.get((side, allele.strip()))
        if a is not None and a.status == "ok" and a.anchor_nt >= 0:
            return a
    return None


def project_junction(strand_seq: str, qlen: int, *, v_row: dict, j_row: dict,
                     v_call: str, j_call: str, anchors: dict,
                     split_checked: bool) -> tuple[Projection | None, str]:
    """Project the junction from two segment hits. Returns ``(projection, refusal_reason)``.

    Args:
        strand_seq: the read on the strand the hits are measured in. For a reverse-complement hit
            the caller passes ``revcomp(read)``, because that is the frame the arithmetic below is
            expressed in and translating coordinates twice is how sign errors get in.
        qlen: length of the ORIGINAL read, which is what ``segmap`` used to build a minus-strand
            ``qstart``. Equal to ``len(strand_seq)``; passed explicitly so a caller that trims cannot
            silently disagree with the mapper.
        v_row, j_row: segment rows -- ``qstart``, ``qend``, ``tstart``, ``split``.
        split_checked: did indel detection actually run? ``segment_rows`` only populates ``split``
            when ``max_indel > 0``, so a ``False`` here means every ``split`` is 0 for want of
            checking, not for want of indels. **No default** -- the caller must state it.

    Exactly one of the return slots is set: a ``Projection`` and ``""``, or ``None`` and a member of
    :data:`REFUSALS`.
    """
    va, ja = _anchor(anchors, "V", v_call), _anchor(anchors, "J", j_call)
    if va is None or ja is None:
        return None, "no_anchor"

    # Locus from the J, never the V -- TRAV/DV alleles pair with either TRAJ or TRDJ and the J
    # decides which, so a V-derived locus would mislabel exactly the reads this check protects.
    if ja.locus in UNVALIDATED_LOCI:
        return None, "unvalidated_locus"

    # An indel between the read and the germline shifts every downstream base, so a projection that
    # assumes a constant offset is wrong by exactly the indel length. `segmap`'s two-diagonal
    # signature is what detects that, and it is the load-bearing refusal here.
    #
    # ⛔ `split` IS ONLY POPULATED WHEN INDEL DETECTION RAN. `segment_rows` passes `max_indel = 0`
    # unless `--indel-rescue` is on, and then every `split` is 0 -- indistinguishable from "checked
    # and clean". A caller that forgets would get the projection with its indel protection SILENTLY
    # INERT, which is this codebase's most repeated failure shape (mmseqs `createdb` on its first
    # byte, `fetch_database` across a filesystem boundary, a reference swap into the wrong cache).
    # So the caller must SAY whether the check ran; there is no default.
    #
    # It matters, measured: on IGH_repertoire (91.77 % median V identity) the gate costs 3.97 % of
    # fast-path yield and takes byte-exact accuracy from **.99634 to .99915** -- 332 wrong junctions
    # down to 74. On IGH_naive (99.41 %) it costs 2.92 % and buys +0.00002. Like `--indel-rescue`
    # itself, its value tracks SHM load, so the right default is per-library, not global.
    if not split_checked:
        return None, "indel_unchecked"
    if v_row.get("split") or j_row.get("split"):
        return None, "indel_split"

    # `segmap` signals a minus-strand hit with qstart > qend, the convention `_align_implied` also
    # relies on. A read whose V and J hits disagree on strand is not a rearrangement this can place.
    v_rc, j_rc = v_row["qstart"] > v_row["qend"], j_row["qstart"] > j_row["qend"]
    if v_rc != j_rc:
        return None, "strand_mismatch"

    def pos(row: dict, anchor_nt: int) -> int:
        # `tstart` is 1-BASED on the forward target (segmap.cpp:376); `anchor_nt` is 0-BASED
        # (cdr3fix.Anchor). Hence the +1. Getting this wrong shifts every junction by one base and
        # still produces something that starts with a plausible codon.
        offset = (anchor_nt + 1) - row["tstart"]
        q = row["qstart"]
        # For a minus-strand hit `qstart` is already in FORWARD coordinates (segmap.cpp:383), so it
        # must be reflected back into the revcomp frame `strand_seq` is in before the offset applies.
        return (qlen - q + 1 if v_rc else q) + offset

    start, fw = pos(v_row, va.anchor_nt), pos(j_row, ja.anchor_nt)
    end = fw + 2                                   # the [FW]118 codon is the junction's last 3 nt

    if end <= start:
        return None, "order"
    if start < 1 or end > len(strand_seq):
        return None, "off_read"

    junction = strand_seq[start - 1:end]
    # A junction whose length is not a multiple of 3 cannot be the thing AIRR says it is, whatever
    # the arithmetic produced. Cheap, and it catches a coordinate error that survives every check
    # above -- which is exactly the class of bug that shipped junctions "short by the truncation".
    if len(junction) % 3:
        return None, "bad_codon"
    return Projection(junction=junction, start=start, end=end, rev_comp=v_rc), ""
