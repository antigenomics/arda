# Multi-node (SLURM) sharding

`src/arda/cluster.py` + CLI `split` / `merge` / `slurm`.

- **split(input, out_dir, shards)**: one streaming pass, round-robin (record k →
  shard k%shards) into `shard_<i>.fasta`. Round-robin balances load and is the
  reason shard sizes differ by ≤1. Writes FASTA (quality dropped; arda ignores it).
- **merge(dir|list, output)**: concatenate per-shard AIRR TSVs, header once.
- **render_submit_script(...)**: emits a `submit.sh` that runs `arda split`, then
  `sbatch --array=0-(N-1) --wrap 'arda annotate ... shard_${SLURM_ARRAY_TASK_ID}'`,
  then `sbatch --dependency=afterok:$ARRAY_JID --wrap 'arda merge'`. `ARDA_MMSEQS`
  is exported into the script if set in the caller's env.
- CLI `arda slurm -i ... -o ... --shards N [--submit]` writes `arda_slurm/submit.sh`
  (chmod +x) and optionally runs it.

Design choice: pre-split once (low I/O) rather than each array task reading the
whole file (round-robin in `read_sequences` would mean N× read I/O on a huge file).

Tested without a cluster: `tests/unit/test_cluster.py` (split partitions all records
disjointly + balanced; merge keeps one header; script chains the three steps with
the dependency). A live SLURM run is still TODO (user will provide a cluster).
Single-node path is unchanged (streaming `annotate_file`).
