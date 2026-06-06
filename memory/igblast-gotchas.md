# IgBLAST gotchas

## VJ loci still require -germline_db_D
Even for D-less loci (TRA, TRG, IGK, IGL), `igblastn` errors without
`-germline_db_D` (it falls back to a built-in `<org>_gl_D` not on our search
path). Fix: `airr_extract._dummy_d_db` builds a one-sequence placeholder D
database and passes it for VJ loci. See `refbuild/airr_extract.py`.

## IGDATA
`igblast.igdata_env()` sets `IGDATA` to `bin/` so IgBLAST finds `internal_data/`
and `optional_file/`. The downloaded release is laid out by
`scripts/fetch_igblast.py` so that `bin/` contains the executables **and** those
two trees.

## Germline DBs
Build with `makeblastdb -parse_seqids -dbtype nucl` from **ungapped** IMGT files
(via `edit_imgt_file.pl`). Gapped IMGT sequences (with `.`) must not go into
makeblastdb.

## AIRR output (`-outfmt 19`)
- Coordinates (`*_start`/`*_end`, `v_sequence_start`…) are **1-based, closed**.
- Region nt seqs are `fwr1..cdr3`; AA versions are `fwr1_aa..cdr3_aa` (we reuse
  these directly instead of re-translating regions).
- `productive` is `T`/`F`/empty; rev_comp `T`/`F`.

## Binaries / versions
IgBLAST 1.22.0 (NCBI LATEST as of build). Organisms with internal_data: human,
mouse, rat, rabbit, rhesus_monkey. MMseqs2 from bioconda (binary `mmseqs`).
