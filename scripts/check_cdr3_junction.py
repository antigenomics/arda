"""Verify CDR3/junction AIRR semantics vs IgBLAST, and trimmed-input behavior."""
from pathlib import Path
import polars as pl
from arda.paths import data_dir
from arda.annotate.mapper import annotate_records
from arda.annotate.io import read_sequences

fasta = data_dir() / "realworld" / "igh_mrna.fasta"
queries = list(read_sequences(fasta))
igb = {r["sequence_id"]: r for r in pl.read_csv(
    fasta.with_suffix(".igblast.airr.tsv"), separator="\t", infer_schema_length=0
).iter_rows(named=True)}
arda = {r["sequence_id"]: r for r in annotate_records(queries, "human", "nt", threads=8)}

# --- AIRR junction invariants + exact match vs igblast ---
n = c_start = fw_end = cdr3_is_inner = junc_match = cdr3_match = 0
bad = []
for sid, ar in arda.items():
    ja, c3 = (ar.get("junction_aa") or ""), (ar.get("cdr3_aa") or "")
    if not ja:
        continue
    n += 1
    c_start += ja.startswith("C")
    fw_end += ja.endswith(("F", "W"))
    cdr3_is_inner += (c3 == ja[1:-1])
    ig = igb.get(sid, {})
    junc_match += (ja == (ig.get("junction_aa") or ""))
    cdr3_match += (c3 == (ig.get("cdr3_aa") or ""))
    if not (ja.startswith("C") and ja.endswith(("F", "W")) and c3 == ja[1:-1]) and len(bad) < 5:
        bad.append((sid, ja, c3, ig.get("junction_aa"), ig.get("cdr3_aa")))

print(f"records with junction: {n}")
print(f"  junction_aa starts with C : {c_start}/{n}")
print(f"  junction_aa ends F/W      : {fw_end}/{n}")
print(f"  cdr3_aa == junction[1:-1] : {cdr3_is_inner}/{n}")
print(f"  junction_aa == igblast    : {junc_match}/{n}")
print(f"  cdr3_aa     == igblast    : {cdr3_match}/{n}")
for sid, ja, c3, ij, ic in bad:
    print(f"  BAD {sid}: arda junc={ja!r} cdr3={c3!r} | ig junc={ij!r} cdr3={ic!r}")

# --- Trimmed inputs: V-only (FR1-FR3) and J-side (CDR3-FR4) fragments ---
print("\n--- trimmed inputs ---")
sid0 = next(s for s in arda if arda[s].get("fwr3_end"))
full = dict(queries)[sid0]
ar_full = arda[sid0]
f1s = int(ar_full["fwr1_start"]); f3e = int(ar_full["fwr3_end"])
c3s = int(ar_full["cdr3_start"]); f4e = int(ar_full["fwr4_end"])
v_only = full[f1s-1:f3e]                 # FR1..FR3, no CDR3/J
j_side = full[c3s-1:f4e]                  # CDR3..FR4, no V
for name, frag in [("V-only(FR1-3)", v_only), ("J-side(CDR3-FR4)", j_side)]:
    r = annotate_records([(name, frag)], "human", "nt", threads=4)[0]
    present = [reg for reg in ("fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4") if r.get(reg)]
    print(f"  {name} (len {len(frag)}): v={r['v_call'][:16]} regions={present} cdr3_aa={r.get('cdr3_aa')!r}")
