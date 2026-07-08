"""Unit tests for the RNA-seq mode (``arda.rnaseq``).

``read_pairs`` is pure and always runs. ``map_rnaseq`` needs mmseqs + the human DB
(skips otherwise). ``correct_airr`` needs the optional ``seqtree`` dep.
"""

from __future__ import annotations

import random

import polars as pl
import pytest

from arda import paths

from arda.rnaseq.map import read_pairs, map_rnaseq, merge_pair
from arda.refbuild.translate import reverse_complement
from tests.conftest import requires_mmseqs, requires_human_db


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
    C1 = "TGTGACAGCAGCTTAGACGGGACAGGGTTC"   # 1 mismatch from P
    U = "TGTACCCCGGGGTTTTAAAACCCCGGGTTC"    # unrelated
    rows = []
    for junc, n in [(P, 100), (C1, 2), (U, 5)]:   # 100 >= 2*20 -> C1 child of P
        for k in range(n):
            rows.append({"sequence_id": f"{junc[:4]}{k}", "junction": junc,
                         "junction_aa": "CASSLDGTF", "v_call": "TRBV20-1*01",
                         "j_call": "TRBJ2-1*01", "locus": "TRB"})
    airr = tmp_path / "in.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    rep = correct_airr(airr, tmp_path / "clones.tsv", read_map=tmp_path / "map.tsv")
    out = pl.read_csv(tmp_path / "clones.tsv", separator="\t", infer_schema_length=0)
    counts = {r["junction"]: int(r["count"]) for r in out.iter_rows(named=True)}

    assert rep.clonotypes_in == 3 and rep.clonotypes_out == 2
    assert counts[P] == 102 and counts[U] == 5   # C1 absorbed into P, read count conserved
    assert C1 not in counts
    rm = pl.read_csv(tmp_path / "map.tsv", separator="\t", infer_schema_length=0)
    assert rm.height == 107 and set(rm["junction"].unique().to_list()) == {P, U}


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
