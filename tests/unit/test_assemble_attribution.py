"""A contig's junction may only be attributed to members that actually COVERED it.

⛔ Contig membership is granted on a ``min_overlap`` match, and after the extension passes have
accumulated germline at the contig ends that match can be **pure germline** -- the 5' pass says so
itself: "that region is shared germline, so any V-read of the gene extends it correctly". Stamping
the contig's junction onto such a member credits a read carrying zero clone-specific evidence to
this clonotype's ``duplicate_count``, and a read of a DIFFERENT clone of the same V gene is exactly
what a germline overlap admits.

⚠ Withholding the row does not lose the read: it stays in the Stage-1 frame and
``correct._assign_coverage`` can still place it on evidence. Only the FORCED attribution is
withheld, and the count is reported as ``members_without_junction``.
"""

from __future__ import annotations

import random

import polars as pl

# A canonical junction: starts C, ends W, in frame, no stop.
JN = "TGT" + "GCTAGA" * 12 + "TGG"
JA = "C" + "AR" * 12 + "W"


def _mock_reannotate(records, organism, threads=0, map_d=True, d_max_evalue=None):
    return [{"sequence_id": cid, "junction": JN, "junction_aa": JA, "v_call": "IGHV3-23*01",
             "j_call": "IGHJ4*02", "locus": "IGH"} for cid, _ in records]


def test_a_germline_only_member_is_not_credited_with_the_contig_junction(tmp_path, monkeypatch):
    """The contig carries the junction only in its middle; a member lying entirely in the 5'
    germline V must not be attributed to it."""
    import arda.rnaseq.assemble as asm

    rng = random.Random(3)
    v = "".join(rng.choice("ACGT") for _ in range(120))      # germline V, shared by every clone
    j = "".join(rng.choice("ACGT") for _ in range(60))
    true = v + JN + j                                        # the real molecule
    monkeypatch.setattr(asm, "reannotate_contigs", _mock_reannotate)

    rows, offs = [], (0, 30, 60, 90, 120, 150, 180)
    for i, o in enumerate(offs):
        cs = len(v) - o
        rows.append({
            "sequence_id": f"f{i}/1", "sequence": true[o:o + 100], "rev_comp": "F", "locus": "IGH",
            "v_call": "IGHV3-23*01", "j_call": "IGHJ4*02", "c_call": "", "c_class": "",
            "cdr3_start": str(cs) if 0 <= cs < 100 else "",
            "junction": "", "junction_aa": "",
        })
    airr = tmp_path / "mapped.airr.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")

    out = tmp_path / "assembled.airr.tsv"
    rep = asm.assemble_contigs(airr, out)
    assert rep.contigs_complete >= 1

    df = pl.read_csv(out, separator="\t", infer_schema_length=0)
    rescued = set(df["sequence_id"].to_list()) if df.height else set()
    # f0 spans true[0:100], entirely inside the 120 nt of germline V -- it never reaches the
    # junction, so it must not be credited.
    assert "f0/1" not in rescued, "a read lying entirely in germline V was credited with the junction"
    # ...and reads that really do cross the junction still are.
    assert rescued, "reads that covered the junction must still be rescued"
    assert set(df["junction_aa"].to_list()) == {JA}


def test_every_rescued_member_actually_overlaps_the_junction(tmp_path, monkeypatch):
    """The invariant, stated directly: no rescued row may come from a span that misses the
    junction by more than the allowed clip."""
    import arda.rnaseq.assemble as asm

    rng = random.Random(5)
    true = ("".join(rng.choice("ACGT") for _ in range(150)) + JN
            + "".join(rng.choice("ACGT") for _ in range(90)))
    monkeypatch.setattr(asm, "reannotate_contigs", _mock_reannotate)

    rows = []
    for i, o in enumerate(range(0, 260, 25)):
        cs = 150 - o
        rows.append({
            "sequence_id": f"r{i}/1", "sequence": true[o:o + 100], "rev_comp": "F", "locus": "IGH",
            "v_call": "IGHV3-23*01", "j_call": "IGHJ4*02", "c_call": "", "c_class": "",
            "cdr3_start": str(cs) if 0 <= cs < 100 else "",
            "junction": "", "junction_aa": "",
        })
    airr = tmp_path / "m.tsv"
    pl.DataFrame(rows).write_csv(airr, separator="\t")
    out = tmp_path / "a.tsv"
    rep = asm.assemble_contigs(airr, out)

    df = pl.read_csv(out, separator="\t", infer_schema_length=0)
    for sid in (df["sequence_id"].to_list() if df.height else []):
        i = int(sid.split("/")[0][1:])
        s, e = i * 25, i * 25 + 100                      # this read's span in `true`
        j0, j1 = 150, 150 + len(JN)                      # the junction's span in `true`
        assert min(e, j1) - max(s, j0) >= asm._MIN_JUNCTION_COVER, (
            f"{sid} was rescued but covers only {min(e, j1) - max(s, j0)} nt of the junction")
    assert rep.reads_rescued == (df.height if df.height else 0)
