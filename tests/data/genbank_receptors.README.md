# genbank_receptors{,_mouse}.fa — provenance

Full-length antibody / T-cell-receptor mRNA records, one FASTA entry each, used by
`tests/unit/test_genbank_receptors.py` and `tests/unit/test_jc_fr4.py`:

- `genbank_receptors.fa` — **29 human** records (IGH 5, IGK 6, IGL 6, TRA 6, TRB 6).
- `genbank_receptors_mouse.fa` — **20 mouse** records (IGH 5, IGK 5, TRA 5, TRB 5). Mouse λ is
  gene-poor and was not queried. Accessions: DQ452620.1 EF154514.1 MG753988.1 MG753989.1
  PX682345.1 PX682346.1 PX682347.1 PX682348.1 PX682349.1 PX682350.1 PX682351.1 PX682352.1
  PX682353.1 PX682354.1 U07658.1 U07659.1 U07660.1 U07661.1 U07662.1 U95921.1

- **Source**: NCBI GenBank (nuccore), fetched 2026-07-08 via E-utilities `efetch` (`rettype=fasta`).
- **Query**, per locus: `"Homo sapiens"[Organism] AND "<chain>"[Title] AND "complete cds"[Title]
  AND biomol_mrna[PROP] AND 700:2500[SLEN]`, `sort=relevance`, top ~6 per locus — where `<chain>` ∈
  {immunoglobulin heavy chain, immunoglobulin kappa, immunoglobulin lambda, T cell receptor beta,
  T cell receptor alpha}. Human-only records kept.
- **Composition**: IGH 5, IGK 6, IGL 6, TRA 6, TRB 6.
- **Accessions**: AY359884.1 BC100294.1 BC110354.1 GQ179995.1 GU122926.1 GU122927.1 GU326331.1
  JN833726.1 JN833727.1 KT207830.1 KT207831.1 OM687940.1 OP311729.1 OP311730.1 OP311731.1
  OR378286.1 PQ177856.1 PQ879411.1 PQ879415.1 PQ879417.1 PQ879419.1 PQ879421.1 PQ879423.1
  PQ879425.1 PQ879427.1 PQ879429.1 PV604442.1 PV604443.1 PX873516.1
- **Re-fetch**: `efetch -db nuccore -id <acc> -format fasta` (Entrez Direct), or
  `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<acc>&rettype=fasta&retmode=text`.
- **Provenance**: experimental — deposited cDNA sequences, not derived or computed by this project.

## Two records that are edge cases, not errors

- **BC100294.1** — a real 5′-truncated TRA cDNA. arda's alignment starts *inside* the J (no V),
  runs into C, so it yields `j_call` + `c_call` but no junction. It is the natural V-less J→C case.
- **PQ879427.1** — a published spike-binding λ clone whose *own GenBank CDS translation* is
  frameshifted after CDR3 (`...YYCQTWGTGTQLGVRRR...`, never reaching FGGGTKLTVL). arda calls it
  non-productive with an out-of-frame `_` marker — correct, and confirmed against the record's own
  `/translation`.
- **OP311729.1 / OP311730.1 / OP311731.1** carry an XhoI (`CTCGAG`) vector cloning site a few nt
  into the J. The V-J scaffold ends at the J terminus and stops the alignment at the site,
  truncating FR4; the `J + C` scaffold anchors in C beyond it and recovers the full FR4. This is
  the FR4-recovery the constant-region reference exists to provide.
