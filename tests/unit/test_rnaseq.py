"""Unit tests for the RNA-seq mode (``arda.rnaseq``).

``read_pairs`` is pure and always runs. ``map_rnaseq`` needs mmseqs + the human DB
(skips otherwise). ``correct_airr`` needs the optional ``seqtree`` dep.
"""

from __future__ import annotations

import gzip
import json
import random

import polars as pl
import pytest
from typer.testing import CliRunner

from arda import __version__, paths
from arda.cli import app

from arda.rnaseq.map import read_pairs, map_rnaseq, merge_pair
from arda.refbuild.translate import reverse_complement
from tests.conftest import requires_mmseqs, requires_human_db, requires_imgt


def test_merge_pair_overlapping():
    # Fragment ACGT...; R1 = first 40, R2 = rc of last 40, overlapping in the middle.
    frag = "ACGTACGTACGTTTGGCCAATTGGCCAAGGTTCCAAGGTTCCAAGATCGATCGATCG"
    r1 = frag[:40]
    r2 = reverse_complement(frag[-40:])  # sequenced from the other end
    assert merge_pair(r1, r2) == frag


def test_merge_pair_no_overlap_returns_none():
    r1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    r2 = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"  # rc is CCCC..., no shared anchor
    assert merge_pair(r1, r2) is None


def test_merge_pair_quality_wins_the_overlap_mismatch():
    """In the overlap the mates disagree at one base; the higher-Phred base wins. Without qualities
    R2 wins the whole overlap (the historical behaviour) -- the regression this fixes."""
    s1 = "AAAACCCCGGGGTTTTACGT"                  # R1
    rcs2 = "AAAACCCCGGGGTTTAACGT"                # rc(R2): identical but s1[15] T -> A
    s2 = reverse_complement(rcs2)
    assert merge_pair(s1, s2) == rcs2            # no quality: R2 wins the overlap (== rc(R2))

    hi = "I" * 20                                # Phred 40
    q2_lowat4 = "I" * 4 + "#" + "I" * 15         # rc(R2) quality at col 15 == q2[::-1][15] == q2[4]
    assert merge_pair(s1, s2, q1=hi, q2=q2_lowat4) == s1        # R1 higher there -> R1 base ('T')

    q1_lowat15 = "I" * 15 + "#" + "I" * 4        # R1 low at the mismatch
    assert merge_pair(s1, s2, q1=q1_lowat15, q2=hi) == rcs2     # R2 higher there -> R2 base ('A')


def test_read_pairs_reconstruct_merges_overlap(tmp_path):
    frag = "ACGTACGTACGTTTGGCCAATTGGCCAAGGTTCCAAGGTTCCAAGATCGATCGATCG"
    r1p = _write_fastq(tmp_path / "r1.fastq", [("f", frag[:40])])
    r2p = _write_fastq(tmp_path / "r2.fastq", [("f", reverse_complement(frag[-40:]))])
    assert list(read_pairs(r1p, r2p, reconstruct=True)) == [("f", frag)]  # merged, bare id


def _write_fastq(path, records):
    with open(path, "w") as fh:
        for sid, seq in records:
            fh.write(f"@{sid}\n{seq}\n+\n{'I' * len(seq)}\n")
    return path


def _truncated_gzip_fastq(path, n_reads=200, chop=100):
    """A valid-prefix, truncated gzip: write ``n_reads``, then chop the last ``chop``
    compressed bytes so the stream ends before its gzip end-of-stream marker."""
    with gzip.open(path, "wt") as fh:
        for i in range(n_reads):
            seq = "ACGT" * 30
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    path.write_bytes(path.read_bytes()[:-chop])
    return path


def test_read_sequences_truncated_gzip_raises_valueerror(tmp_path):
    """A truncated .gz must surface as a clear ValueError, not a bare EOFError.

    A bare EOFError from the gzip layer reaches Typer/Click, which mistakes it for a Ctrl-D
    and prints "Aborted." with no cause -- a truncated FASTQ then looks like a user interrupt
    instead of corrupt input (release blocker, arda 2.3.2)."""
    fq = _truncated_gzip_fastq(tmp_path / "trunc.fq.gz")
    with pytest.raises(ValueError, match="truncated or corrupt"):
        list(read_pairs(fq))  # single-end -> read_sequences


@requires_mmseqs
@requires_human_db
def test_map_rnaseq_truncated_gzip_raises_not_aborts(tmp_path):
    """map_rnaseq must reject a truncated gzip with a ValueError, not die as "Aborted.".

    The reader runs in a daemon thread; its exception is captured and re-raised on the
    consumer side (map.py). This locks the end-to-end contract for the reported blocker."""
    fq = _truncated_gzip_fastq(tmp_path / "trunc.fq.gz")
    with pytest.raises(ValueError, match="truncated or corrupt"):
        map_rnaseq(fq, tmp_path / "out.tsv", threads=1)


@requires_mmseqs
@requires_human_db
def test_annotate_file_truncated_gzip_raises_not_silent(tmp_path):
    """annotate_file shares map's daemon-reader design and must have the same guard.

    Without it, a truncated gzip dies unheard in the reader thread while the sentinel is
    still posted, so annotate_file returns partial output with exit 0 -- silent truncation."""
    from arda.annotate.mapper import annotate_file
    fq = _truncated_gzip_fastq(tmp_path / "trunc.fq.gz")
    with pytest.raises(ValueError, match="truncated or corrupt"):
        annotate_file(fq, tmp_path / "out.tsv", threads=1)


def test_read_pairs_single(tmp_path):
    fq = _write_fastq(tmp_path / "s.fastq", [("a", "ACGT"), ("b", "TTTT")])
    assert list(read_pairs(fq)) == [("a", "ACGT"), ("b", "TTTT")]


def test_read_pairs_paired_tags_mates(tmp_path):
    # Both mates share the SRA id; read_pairs must tag them /1 and /2 and interleave.
    r1 = _write_fastq(tmp_path / "r1.fastq", [("read1", "AAAA"), ("read2", "CCCC")])
    r2 = _write_fastq(tmp_path / "r2.fastq", [("read1", "GGGG"), ("read2", "TTTT")])
    assert list(read_pairs(r1, r2)) == [
        ("read1/1", "AAAA"), ("read1/2", "GGGG"),
        ("read2/1", "CCCC"), ("read2/2", "TTTT"),
    ]


@requires_mmseqs
@requires_human_db
def test_map_rnaseq_filters_and_keys_by_read_id(tmp_path, human_scaffolds):
    # 4 receptor-window reads (real IGH/TRB transcript slices) + 16 random reads.
    rng = random.Random(1)
    pool = [s for i, s in human_scaffolds if i.startswith(("IGH_", "TRB_")) and len(s) > 150]
    recs = []
    for k in range(4):
        s = rng.choice(pool)
        p = rng.randrange(len(s) - 150)
        recs.append((f"rec{k}", s[p:p + 150]))
    for k in range(16):
        recs.append((f"rand{k}", "".join(rng.choice("ACGT") for _ in range(150))))
    rng.shuffle(recs)
    fq = _write_fastq(tmp_path / "in.fastq", recs)

    out = tmp_path / "out.airr.tsv"
    rep = map_rnaseq(fq, out, threads=2, report_path=tmp_path / "rep.json")

    df = pl.read_csv(out, separator="\t", infer_schema_length=0)
    ids = set(df["sequence_id"].to_list())
    assert rep.total_reads == 20
    # Only the 4 receptor windows should survive (recall-first, but random reads must
    # not map): mapped output is exactly the receptor reads, keyed by their read id.
    assert ids == {"rec0", "rec1", "rec2", "rec3"}
    assert rep.mapped_reads == 4
    assert all(v for v in df["v_call"].to_list())  # every kept row has a V call
    # the cutoff is recorded (a mapped_reads count without it is uninterpretable)
    # and every surviving read clears it
    assert rep.min_score == 75.0
    assert all(float(s) >= 75.0 for s in df["mmseqs2_score"].to_list())


def test_map_rnaseq_min_score_gates_low_scoring_reads(tmp_path, human_scaffolds):
    """``--min-score`` must actually drop reads, and ``0`` must disable the filter.

    150 bp receptor windows score far above 75 bits, so an absurd cutoff is the only way
    to prove the gate fires; the paired assertion (``min_score=0`` keeps all 4) proves the
    reads are otherwise mappable and that the gate is what removed them.
    """
    rng = random.Random(1)
    pool = [s for i, s in human_scaffolds if i.startswith(("IGH_", "TRB_")) and len(s) > 150]
    recs = []
    for k in range(4):
        s = rng.choice(pool)
        p = rng.randrange(len(s) - 150)
        recs.append((f"rec{k}", s[p:p + 150]))
    fq = _write_fastq(tmp_path / "in.fastq", recs)

    keep_all = map_rnaseq(fq, tmp_path / "all.tsv", threads=2, min_score=0)
    gated = map_rnaseq(fq, tmp_path / "gated.tsv", threads=2, min_score=10_000)

    assert keep_all.mapped_reads == 4 and keep_all.min_score == 0
    assert gated.mapped_reads == 0 and gated.min_score == 10_000
    assert gated.per_locus == {}  # per-locus counts reflect the filter, not the raw hits


def test_map_rnaseq_max_seqs_changes_calls_not_the_read_set(tmp_path, human_scaffolds):
    """More MMseqs2 candidates may only IMPROVE the best hit, never lose a read.

    A bit score is a maximum over candidates, so widening the candidate list cannot lower it
    and cannot drop a read. This is what makes --max-seqs safe to raise: the filter (and the
    --min-score calibration) are untouched, only which V/J scaffold wins. See mapper._MAX_SEQS.
    """
    rng = random.Random(5)
    pool = [s for i, s in human_scaffolds if i.startswith(("IGH_", "TRB_")) and len(s) > 200]
    recs = []
    for k in range(6):
        s = rng.choice(pool)
        p = rng.randrange(len(s) - 150)
        recs.append((f"rec{k}", s[p:p + 150]))
    fq = _write_fastq(tmp_path / "in.fastq", recs)

    lo = map_rnaseq(fq, tmp_path / "lo.tsv", threads=2, min_score=0, max_seqs=1)
    hi = map_rnaseq(fq, tmp_path / "hi.tsv", threads=2, min_score=0, max_seqs=300)

    a = pl.read_csv(tmp_path / "lo.tsv", separator="\t", infer_schema_length=0)
    b = pl.read_csv(tmp_path / "hi.tsv", separator="\t", infer_schema_length=0)
    assert set(a["sequence_id"]) == set(b["sequence_id"])       # identical read set
    assert lo.mapped_reads == hi.mapped_reads

    bl = {i: float(s) for i, s in zip(a["sequence_id"], a["mmseqs2_score"])}
    bh = {i: float(s) for i, s in zip(b["sequence_id"], b["mmseqs2_score"])}
    assert all(bh[i] >= bl[i] for i in bl), "widening the candidate list lowered a bit score"


def test_prep_rejects_unknown_strand():
    """An unrecognised --strand must fail loudly, not silently search forward-only.

    On stranded paired libraries R2 is antisense, so a typo here quietly discards the R2-only
    fragments -- ~40% of the recoverable repertoire.
    """
    from arda.annotate import mapper

    with pytest.raises(ValueError, match="strand must be"):
        mapper._prep("human", "nt", 1, None, "Both")
    with pytest.raises(ValueError, match="strand must be"):
        mapper._prep("human", "nt", 1, None, "reverse")


def test_correct_parent_child_collapse(tmp_path):
    seqtree = pytest.importorskip("seqtree")  # noqa: F841 — optional dep gate
    from arda.rnaseq.correct import correct_airr

    # 30 nt (in frame) and a canonical C...F junction_aa — `correct` keeps only complete
    # junctions by default, so the fixture has to be a biologically valid one.
    P = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"
    C1 = "TGTGACAGCAGCTTAGACGGGACAGGGTTC"   # 1 substitution from P, same V/J
    U = "TGTACCCCGGGGTTTTAAAACCCCGGGTTC"    # unrelated
    rows = []
    # error_rate 0.01 => a 1-sub child collapses once the parent is >=100x deeper: 500*0.01=5 >= 2.
    for tag, junc, n in [("p", P, 500), ("c", C1, 2), ("u", U, 5)]:
        for k in range(n):
            rows.append({"sequence_id": f"{tag}{k}", "junction": junc,   # unique ids (real reads never collide)
                         "junction_aa": "CASSLDGTF", "v_call": "TRBV20-1*01",
                         "j_call": "TRBJ2-1*01", "locus": "TRB"})
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    rep = correct_airr(airr, tmp_path / "clones.tsv", read_map=tmp_path / "map.tsv")
    out = pl.read_csv(tmp_path / "clones.tsv", separator="\t", infer_schema_length=0)
    counts = {r["junction"]: int(r["duplicate_count"]) for r in out.iter_rows(named=True)}

    assert rep.clonotypes_in == 3 and rep.clonotypes_out == 2
    assert counts[P] == 502 and counts[U] == 5   # C1 absorbed into P, read count conserved
    assert C1 not in counts
    rm = pl.read_csv(tmp_path / "map.tsv", separator="\t", infer_schema_length=0)
    assert rm.height == 507 and set(rm["junction"].unique().to_list()) == {P, U}


def test_correct_drops_incomplete_junctions(tmp_path):
    """A read that stops short of [FW]118 yields a *prefix* of a junction, not a clonotype.

    Stage 1 reports a junction even when the read does not span it, so without this gate the
    clonotype table fills with truncated/out-of-frame/stop-codon artefacts — measured at 42 %
    of IGH clonotypes on 100 bp RNA-seq. `productive` does not catch truncation on its own.
    """
    pytest.importorskip("seqtree")
    from arda.rnaseq.correct import correct_airr

    good = ("TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTF")   # canonical, in frame
    trunc = ("TGTGCCAGCAGCTTAGACGGGACAGGGTAC", "CASSLDGTY")  # no [FW] anchor: read ran out
    stop = ("TGTGCCAGCTGATTAGACGGGACAGGGTTC", "CAS*LDGTF")   # stop codon
    frameshift = ("TGTGCCAGCAGCTTAGACGGGACAGGGTT", "CASSLD_TF")  # 29 nt, out of frame + N

    rows = [{"sequence_id": f"{aa}{k}", "junction": nt, "junction_aa": aa,
             "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01", "locus": "TRB"}
            for nt, aa in (good, trunc, stop, frameshift) for k in range(3)]
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    rep = correct_airr(airr, tmp_path / "clones.tsv")
    out = pl.read_csv(tmp_path / "clones.tsv", separator="\t", infer_schema_length=0)
    assert out["junction"].to_list() == [good[0]]        # only the complete junction survives
    assert rep.reads_with_junction == 12 and rep.reads_incomplete == 9

    # ...and --all-junctions restores the old (wrong) behaviour, proving the gate is what filtered
    raw = correct_airr(airr, tmp_path / "raw.tsv", complete_only=False)
    assert raw.clonotypes_in == 4 and raw.reads_incomplete == 0


def test_correct_keys_on_locus_v_j_not_junction_alone(tmp_path):
    """A clonotype is (locus, v_call, j_call, junction). The SAME nucleotide junction from a
    different V/J is a different clonotype -- grouping on the junction alone merged them and kept
    an arbitrary member's calls."""
    pytest.importorskip("seqtree")
    from arda.rnaseq.correct import correct_airr

    junc, aa = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTF"
    rows = []
    for vj in [("TRBV20-1*01", "TRBJ2-1*01"), ("TRBV28*01", "TRBJ2-7*01")]:   # same junction, diff V/J
        for k in range(4):
            rows.append({"sequence_id": f"{vj[0]}{k}", "junction": junc, "junction_aa": aa,
                         "v_call": vj[0], "j_call": vj[1], "locus": "TRB"})
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    correct_airr(airr, tmp_path / "clones.tsv")
    out = pl.read_csv(tmp_path / "clones.tsv", separator="\t", infer_schema_length=0)
    assert out.height == 2, "same junction with different V/J must stay two clonotypes"
    assert set(out["v_call"].to_list()) == {"TRBV20-1*01", "TRBV28*01"}


def test_correct_read_and_consensus_counts_from_mates(tmp_path):
    """AIRR counts: ``duplicate_count`` = READS (every mate row; the standard read-counting
    convention), ``consensus_count`` = distinct fragment consensuses (the two mates of one molecule
    are one). 5 molecules x 2 spanning mates -> duplicate_count 10, consensus_count 5."""
    pytest.importorskip("seqtree")
    from arda.rnaseq.correct import correct_airr

    junc, aa = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTF"
    rows = []
    for frag in range(5):                                    # 5 molecules, both mates each span
        for mate in ("1", "2"):
            rows.append({"sequence_id": f"read{frag}/{mate}", "junction": junc, "junction_aa": aa,
                         "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01", "locus": "TRB"})
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    correct_airr(airr, tmp_path / "clones.tsv")
    out = pl.read_csv(tmp_path / "clones.tsv", separator="\t", infer_schema_length=0).to_dicts()[0]
    assert int(out["duplicate_count"]) == 10, "10 mate rows are 10 reads"
    assert int(out["consensus_count"]) == 5, "5 molecules are 5 consensuses"


def test_correct_isotype_from_constant_mate(tmp_path):
    """The isotype lives on the constant-region MATE (``c_class``, no junction), which the complete-only
    filter drops; the junction read carries none. correct links them by fragment id and reports the
    dominant RESOLVED class -- the ambiguous ``IGHC`` wins only if nothing resolves."""
    pytest.importorskip("seqtree")
    from arda.rnaseq.correct import correct_airr

    junc, aa = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTF"
    isos = ["IGHG", "IGHG", "IGHG", "IGHC", "IGHC", "IGHC", "IGHC", ""]   # 3 resolved vs 4 ambiguous
    rows = []
    for i, iso in enumerate(isos):                       # /1 = junction read (no c_class), /2 = C mate
        rows.append({"sequence_id": f"f{i}/1", "junction": junc, "junction_aa": aa,
                     "v_call": "IGHV3-23*01", "j_call": "IGHJ4*02", "locus": "IGH", "c_class": ""})
        rows.append({"sequence_id": f"f{i}/2", "junction": "", "junction_aa": "",
                     "v_call": "", "j_call": "", "locus": "", "c_class": iso})
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    correct_airr(airr, tmp_path / "clones.tsv")
    out = pl.read_csv(tmp_path / "clones.tsv", separator="\t", infer_schema_length=0).to_dicts()[0]
    assert out["c_call"] == "IGHG", "3 resolved IGHG (from constant mates) must beat 4 ambiguous IGHC"
    assert int(out["duplicate_count"]) == 8, "8 fragments"


def test_correct_keeps_inframe_indel_variant(tmp_path):
    """A 3 bp (in-frame SHM) indel costs indel_rate**3 and must NOT collapse even onto a very deep
    parent -- three coincident indel errors are ~1e-6, so it is a real clonotype, not an error."""
    pytest.importorskip("seqtree")
    from arda.rnaseq.correct import correct_airr

    P = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"                 # CASSLDGTF (30 nt, in frame)
    V3 = "TGTGCCAGCAGCGACGGGACAGGGTTC"                   # P minus codon 'TTA' -> CASSDGTF (27 nt)
    rows = []
    for tag, junc, aa, n in [("p", P, "CASSLDGTF", 500), ("v", V3, "CASSDGTF", 5)]:
        for k in range(n):
            rows.append({"sequence_id": f"{tag}{k}", "junction": junc, "junction_aa": aa,
                         "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01", "locus": "TRB"})
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    rep = correct_airr(airr, tmp_path / "clones.tsv", max_indel=3)   # search 3 bp indels: found but not collapsed
    junctions = set(pl.read_csv(tmp_path / "clones.tsv", separator="\t",
                                infer_schema_length=0)["junction"].to_list())
    assert rep.clonotypes_out == 2 and V3 in junctions, "in-frame 3 bp indel must survive"


def test_correct_pileup_binom_collapses_low_freq_variant(tmp_path):
    """The binom pileup path piles reads up per position and collapses a low-frequency 1-sub variant
    onto a deep parent (child allele depth is consistent with sequencing error of the parent depth)."""
    pytest.importorskip("seqtree")
    from arda.rnaseq.correct import correct_airr

    P = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"
    C1 = "TGTGACAGCAGCTTAGACGGGACAGGGTTC"                # 1 substitution from P
    rows = []
    for tag, junc, n in [("p", P, 500), ("c", C1, 2)]:
        for m in range(n):                              # full-span reads: sequence == junction
            rows.append({"sequence_id": f"{tag}{m}", "sequence": junc, "rev_comp": "F",
                         "junction": junc, "junction_aa": "CASSLDGTF", "v_call": "TRBV20-1*01",
                         "j_call": "TRBJ2-1*01", "locus": "TRB"})
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    rep = correct_airr(airr, tmp_path / "clones.tsv", error_method="binom")
    assert rep.clonotypes_out == 1, "C1 (2 reads) is a sequencing error of P (500) at its position"


def test_correct_rejects_invalid_error_rates(tmp_path):
    """error_rate/indel_rate must be in (0, 1): p_err < 1 keeps counts strictly increasing along
    parent pointers (no cycles); 0 collapses everything, >=1 collapses nothing."""
    from arda.rnaseq.correct import correct_airr
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame([{"sequence_id": "a", "junction": "TGTGCCAGCAGCTTAGACGGGACAGGGTTC",
                   "junction_aa": "CASSLDGTF", "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01",
                   "locus": "TRB"}]).write_csv(airr, separator="\t")
    for bad in (0.0, 1.0, 1.5, -0.1):
        with pytest.raises(ValueError, match="error_rate"):
            correct_airr(airr, tmp_path / "out.tsv", error_rate=bad)
        with pytest.raises(ValueError, match="indel_rate"):
            correct_airr(airr, tmp_path / "out.tsv", indel_rate=bad)


def test_kmer_is_plumbed_to_mmseqs_and_defaults_to_12(monkeypatch):
    """``-k`` sizes MMseqs2's 4**k prefilter table, so it must actually reach MMseqs2.

    Measured on the V+J | J+C reference, 100 k reads: k=15 (MMseqs2's default) ~8.4 GB, k=13 697 MB,
    k=12 298 MB, k=11 202 MB -- while recall AND precision are invariant over k=11..14 (recall 1.0000,
    precision .9463-.9469). 12 is also the fastest measured. A silent regression here multiplies
    memory with no other symptom, so pin it.
    """
    from arda import mmseqs
    from arda.annotate import mapper

    seen: list[list[str]] = []
    monkeypatch.setattr(mmseqs, "run", lambda args: seen.append(args))

    mmseqs.search("q", "t", "r", "tmp", kmer=12)
    assert "-k" in seen[-1] and seen[-1][seen[-1].index("-k") + 1] == "12"

    mmseqs.search("q", "t", "r", "tmp")            # kmer=None -> let MMseqs2 choose
    assert "-k" not in seen[-1]

    # nucleotide default is 12; the amino-acid prefilter is a different index and must stay untouched
    assert mapper._KMER["nt"] == 12
    assert mapper._KMER["aa"] is None


def test_top_hit_reduces_the_result_db_before_convertalis(monkeypatch):
    """`filterdb --extract-lines 1` must sit between `search` and `convertalis`.

    Without it, `--max-seqs 300` makes convertalis write every alignment's cigar/qaln/taln: 804 k rows
    and 194 MB of TSV for the 4 k hits that survive. Parsing that was arda's single largest memory
    consumer -- 877 MB peak, against 284 MB for the mmseqs subprocess itself. The output is
    bit-identical either way, so only a test can stop this from regressing silently.
    """
    from arda import mmseqs

    seen: list[list[str]] = []
    monkeypatch.setattr(mmseqs, "run", lambda args: seen.append(args))
    mmseqs.top_hit("resDB", "bestDB")
    assert seen[-1] == ["filterdb", "resDB", "bestDB", "--extract-lines", "1"]


@requires_mmseqs
@requires_human_db
def test_default_kmer_keeps_every_receptor_read(tmp_path, human_scaffolds):
    """Dropping k must not cost recall: a shorter seed is strictly more sensitive."""
    rng = random.Random(3)
    pool = [s for i, s in human_scaffolds if i.startswith(("IGH_", "TRB_")) and len(s) > 150]
    recs = [(f"rec{k}", (lambda s, p: s[p:p + 150])(rng.choice(pool), rng.randrange(50)))
            for k in range(5)]
    fq = _write_fastq(tmp_path / "in.fastq", recs)

    default = map_rnaseq(fq, tmp_path / "k13.tsv", threads=2, min_score=0)              # k=13
    mmseqs_default = map_rnaseq(fq, tmp_path / "k15.tsv", threads=2, min_score=0, kmer=None)

    a = pl.read_csv(tmp_path / "k13.tsv", separator="\t", infer_schema_length=0)
    b = pl.read_csv(tmp_path / "k15.tsv", separator="\t", infer_schema_length=0)
    assert set(b["sequence_id"]) <= set(a["sequence_id"]), "k=13 lost a read that k=15 found"
    assert default.mapped_reads >= mmseqs_default.mapped_reads == 5


@requires_mmseqs
@requires_human_db
def test_rnaseq_run_maps_then_corrects_and_merges_the_report(tmp_path, human_scaffolds):
    """``arda rnaseq run`` is ``map`` piped into ``correct``: three named outputs + a merged report.

    The mapping and correction themselves are covered above; this pins the one-shot glue -- the
    ``<prefix>.airr/.clones/.arda.json`` naming from ``--out-prefix``, and that the report carries
    both stages plus the arda version the module actually ran.
    """
    pytest.importorskip("seqtree")
    rng = random.Random(2)
    pool = [s for i, s in human_scaffolds if i.startswith(("IGH_", "TRB_")) and len(s) > 150]
    recs = [(f"rec{k}", (lambda s, p: s[p:p + 150])(rng.choice(pool), rng.randrange(50)))
            for k in range(5)]
    fq = _write_fastq(tmp_path / "in.fastq", recs)

    res = CliRunner().invoke(
        app, ["rnaseq", "run", "--r1", str(fq), "-p", "SAMPLE", "-d", str(tmp_path), "--threads", "2"])
    assert res.exit_code == 0, res.output

    airr, clones, rep = (tmp_path / f"SAMPLE.{ext}" for ext in ("airr.tsv", "clones.tsv", "arda.json"))
    assert airr.exists() and clones.exists() and rep.exists()

    # the clonotype table carries the correct-stage schema even when biology yields few rows.
    # The four D columns are appended by `correct --map-d` (on by default): D is a function of
    # the corrected junction, so it is called once per clonotype rather than voted over reads.
    cols = pl.read_csv(clones, separator="\t", infer_schema_length=0).columns
    assert cols == ["junction", "junction_aa", "v_call", "j_call", "c_call", "locus",
                    "duplicate_count", "consensus_count",
                    "d_call", "d2_call", "d_support", "d2_support"]

    # the merged report carries both stages, at the version the module used (defaults: k=12, min 75)
    r = json.loads(rep.read_text())
    assert r["arda_version"] == __version__
    assert r["map"]["mapped_reads"] == 5 and r["map"]["min_score"] == 75.0
    assert "clonotypes_out" in r["correct"]


# --- constant region: `J + C` scaffolds, the P1 rule, and isotype from a gapped mate -------------

def test_isotype_class_reports_class_never_subclass():
    from arda.refbuild.constant import isotype_class
    assert isotype_class("IGHG1") == "IGHG"
    assert isotype_class("IGHG1,IGHG3") == "IGHG"       # IGHG1-4 are ~95 % identical over CH1
    assert isotype_class("IGKC") == "IGKC"
    assert isotype_class("IGLC1,IGLC2") == "IGLC"       # ambiguous within IGL -> locus constant
    assert isotype_class("IGHG1,IGHM") == "IGHC"        # straddles classes -> locus constant
    assert isotype_class("") == ""                      # no C hit at all is NOT the same as ambiguous


def test_constant_only_predicate_uses_the_vj_end_of_the_scaffold():
    from arda.rnaseq.map import _constant_only
    # a V-J scaffold: vj_end == tlen, so no alignment can start at or past it
    assert not _constant_only({"mmseqs2_t_vjend": 350, "mmseqs2_tstart": 349})
    # a J+C scaffold with a 52 nt J: tstart 53 never touched the J
    assert _constant_only({"mmseqs2_t_vjend": 52, "mmseqs2_tstart": 53})
    assert not _constant_only({"mmseqs2_t_vjend": 52, "mmseqs2_tstart": 40})
    # a reference built before c_call existed carries no vj_end: never drop
    assert not _constant_only({"mmseqs2_t_vjend": "", "mmseqs2_tstart": 900})


# --- Stage 3: contig assembly (the long-CDR3 clones no single read spans) -------------------------

def test_greedy_contigs_reconstructs_a_split_cdr3():
    """A long CDR3 lands across the read ends: V-side reads enter it but truncate, J-side reads
    cover its 3' + J, no single read spans it. Anchored greedy overlap-extension (3' into J, then
    5' into V) must stitch the tiling reads back into the full V(D)J sequence."""
    from arda.rnaseq.assemble import _greedy_contigs

    rng = random.Random(0)
    true = "".join(rng.choice("ACGT") for _ in range(180))     # V(0..70) | CDR3(70..150) | J(150..180)
    cdr3_pos = 70
    reads = [true[o:o + 100] for o in (0, 20, 40, 60, 80, 100)]  # 100 bp reads, stride 20 -> 80 bp overlap
    # cdr3_start within each read (None once the read starts past the CDR3 = a pure J/downstream read)
    cs = [cdr3_pos - o if 0 <= cdr3_pos - o < 100 else None for o in (0, 20, 40, 60, 80, 100)]
    seeds = [i for i, c in enumerate(cs) if c is not None]

    contigs = _greedy_contigs(reads, seeds, cs, k=21, min_overlap=21, min_id=0.9,
                              max_ext_past_cdr3=130, scan_cap=400, min_v=70)
    assert contigs, "the tiling reads must assemble into at least one contig"
    assert any(true in seq for seq, _ in contigs), "the full V(D)J sequence was not reconstructed"


def test_assemble_rescues_incomplete_reads_via_contig(tmp_path, monkeypatch):
    """``assemble_contigs`` attributes a contig's complete junction to its INCOMPLETE member reads
    (a read that already spanned its junction is left for the mapped AIRR, not double-counted), and
    carries each member's own ``c_class`` so isotype survives. The junction call itself
    (``reannotate_contigs``) is mmseqs-backed and mocked here to keep the test pure."""
    import arda.rnaseq.assemble as asm

    rng = random.Random(1)
    true = "".join(rng.choice("ACGT") for _ in range(180))
    offs = (0, 20, 40, 60, 80, 100)
    rows = []
    for i, o in enumerate(offs):
        cs = 70 - o
        rows.append({
            "sequence_id": f"f{i}/1", "sequence": true[o:o + 100], "rev_comp": "F", "locus": "IGH",
            "v_call": "IGHV3-23*01", "j_call": "IGHJ4*02", "c_call": "",
            "c_class": "IGHG" if o >= 80 else "",                 # the 3' reads carry the isotype
            "cdr3_start": str(cs) if 0 <= cs < 100 else "",
            "junction": "", "junction_aa": "",                    # every read is INCOMPLETE
        })
    airr = tmp_path / "mapped.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    JN, JA = "TGT" + "GCTAGA" * 12 + "TGG", "C" + "AR" * 12 + "W"   # canonical, in frame, no stop
    monkeypatch.setattr(asm, "reannotate_contigs", lambda records, organism, threads=0, map_d=True: [
        {"sequence_id": cid, "junction": JN, "junction_aa": JA, "v_call": "IGHV3-23*01",
         "j_call": "IGHJ4*02", "locus": "IGH", "d_call": "IGHD3-10*01", "d_support": "0.01",
         "np1": "GG", "np2": "TT"} for cid, _ in records])

    out = tmp_path / "assembled.airr.tsv"
    rep = asm.assemble_contigs(airr, out)
    assert rep.contigs >= 1 and rep.contigs_complete >= 1
    df = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert df.height >= 2, "incomplete member reads should be rescued with the contig junction"
    assert set(df["junction_aa"].to_list()) == {JA}
    assert "IGHG" in df["c_class"].to_list(), "the 3' member's isotype must survive onto the contig"
    # The contig is the only thing that spans a long CDR3, so its D must travel to the reads.
    assert set(df["d_call"].to_list()) == {"IGHD3-10*01"}
    assert set(df["np1"].to_list()) == {"GG"}
    assert "d_sequence_start" not in df.columns, "contig coords are meaningless on a member read"


def test_constant_rule_is_per_fragment_and_donates_isotype_from_a_gapped_mate():
    """Insert size exceeds 2x read length for 36 % of pairs, so the commonest informative layout is
    R1 across V/CDR3 and R2 deep in C with no J. Read-level filtering throws that R2 -- and the only
    isotype evidence the fragment has -- away."""
    from arda.rnaseq.map import _apply_constant_rule
    vdj = {"sequence_id": "f1/1", "mmseqs2_t_vjend": 350, "mmseqs2_tstart": 10,
           "c_call": "", "c_class": ""}
    cmate = {"sequence_id": "f1/2", "mmseqs2_t_vjend": 52, "mmseqs2_tstart": 80,
             "c_call": "IGHG1,IGHG3", "c_class": "IGHG"}
    # a fragment lying entirely inside C: real receptor mRNA, but no rearrangement
    conly_a = {"sequence_id": "f2/1", "mmseqs2_t_vjend": 52, "mmseqs2_tstart": 60,
               "c_call": "IGHM", "c_class": "IGHM"}
    conly_b = {"sequence_id": "f2/2", "mmseqs2_t_vjend": 52, "mmseqs2_tstart": 90,
               "c_call": "IGHM", "c_class": "IGHM"}

    kept, dropped, donated = _apply_constant_rule([vdj, cmate, conly_a, conly_b])

    assert [r["sequence_id"] for r in kept] == ["f1/1"]   # the C-only mate is not itself an AIRR row
    assert dropped == 1 and donated == 1
    assert kept[0]["c_call"] == "IGHG1,IGHG3"             # ...but it donated the isotype
    assert kept[0]["c_class"] == "IGHG"                   # class, never subclass


def test_read_pairs_rejects_shuffled_or_truncated_mates(tmp_path):
    """A shuffled or truncated R2 silently pairs mate 1 of one fragment with mate 2 of another.
    That exact corruption produced a published false discovery in this project."""
    import pytest
    from arda.rnaseq.map import read_pairs
    r1 = tmp_path / "r1.fq"; r2 = tmp_path / "r2.fq"; r2s = tmp_path / "r2s.fq"; r2t = tmp_path / "r2t.fq"
    rec = lambda n, s: f"@{n}\n{s}\n+\n{'I' * len(s)}\n"  # noqa: E731
    r1.write_text(rec("a", "ACGT") + rec("b", "TTTT"))
    r2.write_text(rec("a", "ACGT") + rec("b", "TTTT"))
    r2s.write_text(rec("b", "TTTT") + rec("a", "ACGT"))   # shuffled
    r2t.write_text(rec("a", "ACGT"))                       # truncated

    assert len(list(read_pairs(r1, r2))) == 4              # the good case still works
    with pytest.raises(ValueError, match="mate mismatch"):
        list(read_pairs(r1, r2s))
    with pytest.raises(ValueError, match="truncated"):
        list(read_pairs(r1, r2t))


@requires_imgt
@pytest.mark.skipif(not (paths.database_dir() / "c_genes" / "human.fasta").exists(),
                    reason="C-gene bundle not present")
def test_jc_scaffolds_translate_through_the_splice_junction():
    """`CExon1` begins MID-CODON -- the codon straddles the J-C splice -- so `J + CExon1`
    reconstructs the mRNA and reads through. If the exon boundary or the frame were wrong, it
    would not. This is the acceptance test for the whole constant-region reference."""
    from arda.refbuild.constant import build_jc_scaffolds, isotype_class
    from arda.refbuild.translate import translate

    jc = build_jc_scaffolds("human", "Homo_sapiens")
    assert len(jc) > 300
    assert {s.locus for s in jc} == {"IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"}

    for j, c, expect in [("IGHJ4*02", "IGHG1*01", "WGQGTLVTVSSASTKGPSVFP"),
                         ("IGHJ4*02", "IGHM*01", "WGQGTLVTVSSGSASAPTLFP"),
                         ("IGKJ1*01", "IGKC*01", "WTFGQGTKVEIKRTVAAPSVFI")]:
        hits = [s for s in jc if j in s.j_call.split(",") and c in s.c_call.split(",")]
        assert hits, f"no scaffold {j}+{c}"
        s = hits[0]
        assert any(expect in translate(s.sequence[f:]) for f in range(3)), f"{j}+{c} frameshifted"
        assert s.sequence[:s.j_len] and s.j_len < len(s.sequence)   # J part, then C part

    # the whole point: this must stay a rounding error next to the 17,244 V-J scaffolds
    assert len(jc) < 600, "J+C scaffolds should be additive, not multiplicative"
    assert isotype_class(next(s for s in jc if "IGHG1*01" in s.c_call).c_call) == "IGHG"


@requires_imgt
@pytest.mark.skipif(not (paths.database_dir() / "c_genes").exists(),
                    reason="C-gene bundle not present")
def test_every_shipped_c_gene_maps_to_a_known_locus():
    """A C gene whose name does not resolve to a locus is dropped *silently* and its whole locus loses
    its constant region. Mouse shipped `TCRG-C1..C4` (an alternate nomenclature) where IMGT names them
    `TRGC1..4`, so mouse TRG had no J+C scaffolds and nothing said so. Assert the mapping is total."""
    from arda.refbuild.constant import _locus_of, build_jc_scaffolds
    from arda.refbuild.imgt import read_fasta
    from arda.refbuild.loci import LOCI

    known = {loc.name for loc in LOCI}
    for fa in sorted((paths.database_dir() / "c_genes").glob("*.fasta")):
        for header, _ in read_fasta(fa):
            gene = header.split()[0]
            assert _locus_of(gene) in known, f"{fa.name}: {gene} -> {_locus_of(gene)!r} is not a locus"

    # mouse must reach all seven loci, human too -- the regression that motivated this test
    for organism, species_dir, loci in (("human", "Homo_sapiens", 7), ("mouse", "Mus_musculus", 7)):
        jc = build_jc_scaffolds(organism, species_dir)
        assert len({s.locus for s in jc}) == loci, f"{organism}: {sorted({s.locus for s in jc})}"


# `J + CH1` must translate contiguously across the splice -- CH1 begins mid-codon, so this is the one
# property that proves the exon's 5' boundary is right. Eight functional C genes fail it: their CH1
# anchor in the source is 1-2 nt short (or long), which is provable from a paralog -- mouse TRBC1 and
# TRBC2 splice onto the SAME J cluster and so share a splice phase, yet TRBC2 reads through and TRBC1
# does not, and TRBC1's sequence is TRBC2's minus one leading base. The missing bases are not
# reconstructed here: inventing sequence to satisfy a test is worse than a documented defect.
#
# Consequence: for these genes the scaffold is `J + C[k:]`, so a read crossing the splice pays a k-nt
# gap. It still aligns and still carries the right `c_call`; only its bit score dips. Human is clean.
# If a data fix lands, this test FAILS -- delete the entry, that is the point of pinning it.
KNOWN_BAD_CH1_START = {
    ("mouse", "TRBC1*01"), ("rat", "TRBC1*01"), ("rat", "TRBC2*01"),
    ("rabbit", "IGHA5*01"), ("rabbit", "IGLC4*01"), ("rabbit", "IGLC5*01"), ("rabbit", "IGLC6*01"),
    ("rhesus_monkey", "IGHM*01"),
}
_SPECIES_DIRS = {"human": "Homo_sapiens", "mouse": "Mus_musculus", "rat": "Rattus_norvegicus",
                 "rabbit": "Oryctolagus_cuniculus", "rhesus_monkey": "Macaca_mulatta"}


@requires_imgt
@pytest.mark.skipif(not (paths.database_dir() / "c_genes").exists(),
                    reason="C-gene bundle not present")
def test_jc_scaffolds_read_through_the_splice_for_every_functional_c_gene():
    import re
    from arda.refbuild.constant import build_jc_scaffolds
    from arda.refbuild.imgt import read_fasta
    from arda.refbuild.translate import translate

    # IMGT J-TRP (W118) opens FR4 in IGH; J-PHE (F118) in every other locus. Both match [FW]G.G.
    motif = re.compile(r"[FW]G.G")
    func = {}
    for fa in (paths.database_dir() / "c_genes").glob("*.fasta"):
        for hdr, _ in read_fasta(fa):
            parts = hdr.split()
            func[(fa.stem, parts[0])] = parts[1].split("=")[1] if len(parts) > 1 else "?"

    # per functional C gene: how many of its scaffolds read through, out of how many
    tally: dict[tuple[str, str], list[int]] = {}
    for organism, species_dir in _SPECIES_DIRS.items():
        for sc in build_jc_scaffolds(organism, species_dir):
            for f in range(3):
                j_aa = translate(sc.sequence[f:sc.j_len])
                if not (motif.search(j_aa) and "*" not in j_aa):
                    continue                      # not the J's own reading frame
                reads_through = "*" not in translate(sc.sequence[f:])[max(0, len(j_aa) - 6):len(j_aa) + 15]
                for gene in sc.c_call.split(","):
                    if func.get((organism, gene)) != "F":
                        continue                  # a pseudogene's CH1 need not be an open frame
                    t = tally.setdefault((organism, gene), [0, 0])
                    t[0] += 1
                    t[1] += reads_through
                break

    # A BOUNDARY defect makes a gene fail on EVERY J: the wrong exon start is a property of the gene.
    # A single J failing against an otherwise-clean gene is a property of that J allele (odd length or
    # an internal stop) and is caught by the aggregate below, not misattributed to the C gene.
    seen_bad = {g for g, (n, ok) in tally.items() if n and ok == 0}
    assert seen_bad == KNOWN_BAD_CH1_START, (
        f"CH1 5'-boundary defects changed.\n  newly broken: {sorted(seen_bad - KNOWN_BAD_CH1_START)}"
        f"\n  now fixed (remove from KNOWN_BAD_CH1_START): {sorted(KNOWN_BAD_CH1_START - seen_bad)}")

    # Away from the boundary defects, read-through is near-total. The residue (human TRDJ4*02+TRDC,
    # mouse TRAJ45*02+TRAC) is a property of those J alleles' splice phase, not of the C exon -- human
    # TRDC reads through on 4 of its 5 J alleles. It does not affect the NUCLEOTIDE scaffold, which is
    # the only thing arda aligns against; translation is a diagnostic for the C exon's 5' boundary.
    checked = sum(n for g, (n, _) in tally.items() if g not in KNOWN_BAD_CH1_START)
    passed = sum(ok for g, (_, ok) in tally.items() if g not in KNOWN_BAD_CH1_START)
    assert checked > 600, checked
    assert passed / checked > 0.99, f"only {passed}/{checked} functional scaffolds read through"

    # no gene-wide boundary defect may appear in human, the reference species
    assert not any(org == "human" for org, _ in seen_bad), sorted(seen_bad)


def test_clonotype_row_order_is_deterministic_under_tied_abundance(tmp_path):
    """Tied clonotypes must not reorder between runs.

    Ranking on `(duplicate_count, consensus_count)` alone left ties in *read* order, and read
    order comes out of a threaded mmseqs search -- so the same FASTQ produced the same rows in
    a different sequence each run, and `examples/rnaseq/clones.tsv` was not byte-reproducible.
    Shuffling the input reads must not move a single output row.
    """
    import random

    from arda.rnaseq.correct import correct_airr

    # Three distinct clonotypes with identical read support: every pairwise tie is live.
    juncs = ["TGTGCCAGCAGCTTAGACGGGACAGGGTTC",
             "TGTGCCAGCAGCTTAGACGGGACAGGTTTC",
             "TGTGCCAGCAGCTTAGACGGGACAGGCTTC"]
    rows = [{"sequence_id": f"r{i}_{k}", "junction": jn, "junction_aa": "CASSLDGTF",
             "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01", "locus": "TRB"}
            for i, jn in enumerate(juncs) for k in range(6)]

    outs = []
    for seed in (0, 1, 2):
        shuffled = rows[:]
        random.Random(seed).shuffle(shuffled)
        airr, clones = tmp_path / f"in{seed}.tsv", tmp_path / f"out{seed}.tsv"
        pl.DataFrame(shuffled).write_csv(airr, separator="\t")
        correct_airr(airr, clones, map_d=False)
        outs.append(clones.read_text())

    assert outs[0] == outs[1] == outs[2], "clonotype table must be byte-stable under read order"
    assert len(pl.read_csv(outs[0].encode(), separator="\t", infer_schema_length=0)) == 3
