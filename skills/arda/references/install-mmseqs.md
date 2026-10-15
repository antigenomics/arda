# MMseqs2 Setup & Auto-Fetch

Annotation requires the `mmseqs` binary. arda makes this hands-off.

**Contents:** discovery order · auto-fetch · env vars · shipped indexes ·
version mismatch.

## Discovery order

`arda.mmseqs.mmseqs_binary()` resolves, in order:

1. `$ARDA_MMSEQS` — explicit path to a binary (highest priority).
2. `<project>/bin/mmseqs` — a binary placed in the arda checkout's `bin/`.
3. `mmseqs` on `PATH` — e.g. the bioconda build in the `arda` conda env.
4. **Auto-fetch** — download a static MMseqs2 binary into `bin/mmseqs` (one-time).

So both conda users (env ships `mmseqs2`) and plain-pip users get a working
binary with no manual step.

## Auto-fetch

On first use with nothing found, arda downloads the platform-appropriate static
release from the MMseqs2 GitHub `latest` release into `bin/mmseqs`. Logic lives in
the packaged module `arda._mmseqs_fetch` (so it works for wheel installs), with
`scripts/fetch_mmseqs.py` as a CLI wrapper:

```bash
python scripts/fetch_mmseqs.py            # eager fetch into bin/
python scripts/fetch_mmseqs.py --force    # re-download
```

`setup.sh --no-conda` runs the eager fetch for you.

Implementation note: the downloader sends a browser `User-Agent` — GitHub release
downloads return HTTP 504 under the default `Python-urllib` UA on some networks.

## Environment variables

| Variable | Effect |
|----------|--------|
| `ARDA_MMSEQS` | Absolute path to a specific mmseqs binary (skips all other lookup). |
| `ARDA_MMSEQS_ASSET` | Override the release asset, e.g. `mmseqs-linux-sse41.tar.gz` on pre-AVX2 CPUs. |
| `ARDA_NO_AUTO_FETCH` | Disable auto-fetch; then install mmseqs yourself. |

Asset defaults: macOS → `mmseqs-osx-universal.tar.gz`; Linux x86-64 →
`mmseqs-linux-avx2.tar.gz`; Linux aarch64 → `mmseqs-linux-arm64.tar.gz`.

## Shipped indexes & version mismatch

`database/vdj/<organism>/mmseqs/<nt|aa>/` can hold immutable index generations.
Reuse requires the source FASTA SHA256, normalized MMseqs version, and SHA256 of
all sequence/header/index files to match the generation manifest. Legacy `VERSION`
and mtime markers alone cannot establish identity; old files remain untouched,
while a validated private generation is built under `data/mmseqs_db/`.
`arda build-index --organism all` builds validated generations in the source tree.

`arda.reference_index_identity(organism, seqtype="aa")` returns a JSON-compatible
mapping with `source_fasta` (resolved path), `source_fasta_sha256`, `database_path`
(exact resolved generation prefix), `files_sha256` (sidecar basename → SHA256),
`mmseqs_version`, `dbtype`, `schema_version`, `organism`, and `seqtype`.
It validates or builds an index without running an annotation search. Capture it
before and after a batch and require equality to pin that batch's reference.
