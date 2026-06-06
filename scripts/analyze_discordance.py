"""Categorize arda-vs-IgBLAST discordances and confirm V/J output."""
from pathlib import Path
import polars as pl
from arda.paths import data_dir
from arda.annotate.mapper import annotate_records
from arda.annotate.io import read_sequences

REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3")
fasta = data_dir() / "realworld" / "igh_mrna.fasta"
queries = list(read_sequences(fasta))
igb = {r["sequence_id"]: r for r in pl.read_csv(
    fasta.with_suffix(".igblast.airr.tsv"), separator="\t", infer_schema_length=0
).iter_rows(named=True)}
arda = {r["sequence_id"]: r for r in annotate_records(queries, "human", "nt", threads=8)}

# --- V/J gene-level agreement (arda best-hit allele set vs igblast) ---
def gene(call): return (call or "").split("*")[0].split(",")[0]
vj_match = vj_tot = 0
for sid, ig in igb.items():
    ar = arda.get(sid)
    if not ar or not ar["v_call"]:
        continue
    vj_tot += 1
    arda_v_genes = {g.split("*")[0] for g in ar["v_call"].split(",")}
    vj_match += gene(ig["v_call"]) in arda_v_genes
print(f"V gene agreement (igblast V in arda best-hit set): {vj_match}/{vj_tot}")
print("sample v_call/j_call output:")
for sid in list(igb)[:3]:
    ar = arda.get(sid, {})
    print(f"  {sid}: arda v={ar.get('v_call','')[:30]} j={ar.get('j_call','')[:20]} | "
          f"igblast v={igb[sid]['v_call']} j={igb[sid]['j_call']}")

# --- Categorize discordances ---
cats = {"exact": 0, "boundary(substr)": 0, "frameshift": 0, "no_arda": 0, "other": 0}
examples = {"frameshift": [], "other": []}
for sid, ig in igb.items():
    ar = arda.get(sid)
    for r in REGIONS:
        i = (ig.get(f"{r}_aa") or "").strip()
        if not i:
            continue
        a = ((ar or {}).get(f"{r}_aa") or "").strip()
        if not a:
            cats["no_arda"] += 1
            continue
        if i == a:
            cats["exact"] += 1
        elif i in a or a in i:
            cats["boundary(substr)"] += 1
        elif a.count("*") >= 1 or (len(a) == len(i) and sum(x != y for x, y in zip(a, i)) > len(i) // 2):
            cats["frameshift"] += 1
            if len(examples["frameshift"]) < 4:
                examples["frameshift"].append((sid, r, i, a))
        else:
            cats["other"] += 1
            if len(examples["other"]) < 4:
                examples["other"].append((sid, r, i, a))
print("\ndiscordance categories (region-level):")
for k, v in cats.items():
    print(f"  {k:18s} {v}")
for cat, exs in examples.items():
    if exs:
        print(f"\n{cat} examples:")
        for sid, r, i, a in exs:
            print(f"  {sid} {r}: IG={i!r} ARDA={a!r}")
