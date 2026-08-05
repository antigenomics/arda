#!/bin/bash
# Pack the curated VDJ reference into the release asset that arda's auto-fetch downloads.
#
# `arda._database_fetch` fetches
#   https://github.com/antigenomics/arda/releases/download/v<version>/arda-reference-vdj.tar.gz
# on first use of a plain `pip install arda-mapper`, and extracts a tarball whose root holds `vdj/`
# (i.e. built with `tar -C database ... vdj`). It ships the allele FASTAs + region markup only; the
# version-specific precompiled MMseqs2 indexes (`vdj/<org>/mmseqs/`) are NOT shipped -- arda rebuilds
# them on demand for the locally installed mmseqs version. So exclude them (this is what the working
# v2.0.3 asset did: 3.1 MB, not ~25 MB).
#
# Usage:  bash scripts/pack_reference.sh [output.tar.gz]
set -euo pipefail
OUT="${1:-arda-reference-vdj.tar.gz}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export COPYFILE_DISABLE=1   # macOS bsdtar: no ._AppleDouble members

tar czf "$OUT" -C "$ROOT/database" --exclude='vdj/*/mmseqs' --exclude='vdj/*/segments.*' vdj
# Capture the listing once. Piping `tar tzf | grep -q` SIGPIPEs tar under `pipefail` (grep exits on the
# first match, closing the pipe) -> a false failure; here-strings have no producer to break.
members=$(tar tzf "$OUT")
echo "wrote $OUT ($(du -h "$OUT" | cut -f1)); $(grep -c . <<<"$members") members"
# Every file a pip user needs at runtime must be in here, per organism. `arda markup`,
# `arda.dpost` and D mapping all read these, and a missing one is a 404-shaped failure that
# only shows up on someone else's machine.
for org in human mouse rat rabbit rhesus_monkey; do
  for f in alleles.fasta alleles.aa.fasta markup.tsv markup.aa.tsv cdr3_anchors.tsv; do
    grep -q "^vdj/$org/$f\$" <<<"$members" || { echo "ERROR: vdj/$org/$f missing"; exit 1; }
  done
done
# D germlines exist only for organisms with a D locus; d_prior only where a model was published.
for org in human mouse; do
  grep -q "^vdj/$org/d_germlines.fasta\$" <<<"$members" || { echo "ERROR: vdj/$org/d_germlines.fasta missing"; exit 1; }
  grep -q "^vdj/$org/d_prior.tsv\$" <<<"$members" || { echo "ERROR: vdj/$org/d_prior.tsv missing"; exit 1; }
done
if grep -q 'mmseqs' <<<"$members"; then echo "ERROR: mmseqs indexes leaked into the asset"; exit 1; fi
echo "OK: asset root holds vdj/ (markup + anchors + D priors), no mmseqs indexes"
