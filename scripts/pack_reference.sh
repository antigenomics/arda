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
# Use the pinned manifest, not the working directory's generated files or OS metadata.
python3 - "$ROOT" "$OUT" <<'PYTHON'
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

root, output = map(Path, sys.argv[1:])
manifest = json.loads((root / 'src/arda/_reference_manifest.json').read_text())['files_sha256']
payloads = []
for name, expected in sorted(manifest.items()):
    data = (root / 'database' / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise SystemExit(f'ERROR: reference does not match pinned manifest: {name}')
    payloads.append((name, data))
with output.open('wb') as stream, gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as compressed:
    with tarfile.open(fileobj=compressed, mode='w|', format=tarfile.USTAR_FORMAT) as archive:
        for name, data in payloads:
            member = tarfile.TarInfo(name)
            member.size = len(data)
            member.mode = 0o644
            archive.addfile(member, io.BytesIO(data))
print(f'wrote {output} ({output.stat().st_size} bytes); {len(payloads)} pinned runtime files')
PYTHON
