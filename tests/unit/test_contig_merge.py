"""The two contig-cigar paths agree — merge reproduces re-annotate.

A Stage-3 contig can get its AIRR cigars two ways (see :mod:`arda.annotate.contig`):
re-annotate the whole contig, or stitch its reads' alignments (the C++
``_markup.merge_alignment``). Both end in ``transfer_hit``, so if the merge reconstructs the
same query->scaffold alignment the re-annotation produced, the records are byte-identical.

That is exactly what :func:`test_merge_reconstructs_the_reannotated_alignment` checks, on all 29
real human GenBank receptors: take each contig's re-annotated alignment, cut it into overlapping
read-sized windows (a stand-in for the assembly layout), and assert the merge stitches them back
into the original alignment. The scale/timing comparison lives in the arda-benchmark Phase-D harness.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from arda import _markup, paths
from arda.annotate.cigar import check_cigar
from arda.annotate.contig import Contig, ReadPlacement, merge_contig
from arda.annotate.mapper import annotate_records
from arda.annotate.reference import load_reference

from tests.conftest import requires_mmseqs, requires_human_db

_DATA = Path(__file__).parent.parent / "data"


def _load(fname: str) -> list[tuple[str, str]]:
    recs, name, seq = [], None, []
    for line in (_DATA / fname).read_text().splitlines():
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name:
        recs.append((name, "".join(seq)))
    return recs


def _slice_into_reads(qaln: str, taln: str, qstart: int, tstart: int,
                      n_windows: int = 4) -> list[ReadPlacement]:
    """Cut a query->scaffold alignment into overlapping read windows (an assembly layout stand-in).

    Windows begin and end on aligned (M) columns so each read's offset is exact, and overlap so
    every column is covered by at least one read -- the contiguous-consensus layout the merge assumes.
    """
    cols = list(zip(qaln, taln))
    # 1-based contig / scaffold position at the START of each column (cumulative non-gap counts).
    qcum = [0] * (len(cols) + 1)
    tcum = [0] * (len(cols) + 1)
    for i, (a, b) in enumerate(cols):
        qcum[i + 1] = qcum[i] + (a != "-")
        tcum[i + 1] = tcum[i] + (b != "-")
    m_cols = [i for i, (a, b) in enumerate(cols) if a != "-" and b != "-"]
    total = len(m_cols)
    if total < 2 * n_windows:                       # too few aligned columns to tile: one window
        return [ReadPlacement(qaln, taln, 1, tstart, qstart - 1)]

    win = max(3, (total * 2) // n_windows)          # ~50 % overlap between neighbours
    step = max(1, (total - win) // (n_windows - 1))
    starts = list(range(0, max(1, total - win + 1), step))
    if starts[-1] != total - win:
        starts.append(total - win)

    reads = []
    for s in starts:
        p = m_cols[s]                               # first column of the window (an M)
        q = m_cols[min(s + win, total - 1)]         # last column of the window (an M)
        a, b = p, q + 1
        reads.append(ReadPlacement(
            qaln=qaln[a:b], taln=taln[a:b],
            qstart=1,                               # read starts at its first aligned base
            tstart=tstart + tcum[a],                # 1-based scaffold pos at column a
            offset=qstart - 1 + qcum[a],            # 0-based contig pos of the read's first base
        ))
    return reads


pytestmark = [requires_mmseqs, requires_human_db]


@pytest.fixture(scope="module")
def human_annot():
    recs = _load("genbank_receptors.fa")
    assert len(recs) == 29
    return {r["sequence_id"]: r for r in annotate_records(recs, "human", "nt", threads=8)}


def test_merge_reconstructs_the_reannotated_alignment(human_annot):
    """The core agreement: for every real contig, cutting its re-annotated alignment into
    overlapping reads and merging them yields the SAME (qaln, taln, qstart, qend, tstart, tend).
    Identical hit => identical transfer_hit output => the two paths produce the same AIRR record."""
    checked = 0
    for acc, r in human_annot.items():
        qaln, taln = r["sequence_alignment"], r["germline_alignment"]
        if not qaln:
            continue
        qs, qe = int(float(r["mmseqs2_qstart"])), int(float(r["mmseqs2_qend"]))
        ts, te = int(float(r["mmseqs2_tstart"])), int(float(r["mmseqs2_tend"]))
        reads = _slice_into_reads(qaln, taln, qs, ts)
        assert len(reads) >= 1
        merged = _markup.merge_alignment(
            [x.qaln for x in reads], [x.taln for x in reads],
            [x.qstart for x in reads], [x.tstart for x in reads],
            [x.offset for x in reads], r["sequence"])
        assert merged == (qaln, taln, qs, qe, ts, te), f"{acc}: merge != re-annotate alignment"
        checked += 1
    assert checked >= 25


def test_merge_contig_end_to_end_matches_across_read_splits():
    """The full merge_contig wrapper (C++ merge -> transfer_hit) on a real scaffold as the contig:
    valid cigars, correct calls, and -- the stitching invariant -- one read and two overlapping
    reads give a byte-identical record. Uses a reference scaffold so the target is known offline."""
    from arda.refbuild.imgt import read_fasta

    ref = load_reference("human", "nt")
    seqs = dict((h.split()[0], s) for h, s in read_fasta(ref.target_fasta))
    sid = next(s for s, e in ref.entries.items()
               if e.v_call and e.j_call and e.vj_end and e.v_sequence_end
               and e.j_sequence_start and s in seqs and len(seqs[s]) > 120)
    seq = seqs[sid]

    whole = Contig(sid + "_contig", seq, sid, [ReadPlacement(seq, seq, 1, 1, 0)])
    rec = merge_contig(whole, ref)
    assert rec["v_call"] == ref.entries[sid].v_call
    assert rec["j_call"] == ref.entries[sid].j_call
    assert rec["sequence_alignment"] == seq
    assert rec["v_cigar"] and check_cigar(rec["v_cigar"], len(seq))
    assert rec["j_cigar"] and check_cigar(rec["j_cigar"], len(seq))

    mid = len(seq) // 2
    r1 = ReadPlacement(seq[: mid + 20], seq[: mid + 20], 1, 1, 0)
    r2 = ReadPlacement(seq[mid - 20:], seq[mid - 20:], 1, (mid - 20) + 1, mid - 20)
    rec2 = merge_contig(Contig(sid + "_contig", seq, sid, [r1, r2]), ref)
    assert rec2 == rec, "two overlapping reads must stitch to the same record as one"


def test_merge_alignment_rejects_a_coverage_gap():
    """A scaffold column inside the aligned span that no read covers is not one contig -> throws,
    rather than silently emitting a bogus alignment."""
    # two reads on scaffold cols 1..3 and 7..9, nothing on 4..6
    with pytest.raises(Exception):
        _markup.merge_alignment(["ACG", "TAC"], ["ACG", "TAC"], [1, 1], [1, 7], [0, 6], "ACGNNNTAC")
