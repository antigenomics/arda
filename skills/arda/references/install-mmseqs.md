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

`database/vdj/<organism>/mmseqs/<nt|aa>/` holds precompiled MMseqs2 indexes plus a
`VERSION` marker. They are used out of the box **only when the local mmseqs
version matches**; otherwise arda builds a private DB into `data/mmseqs_db/` on
first run and reuses it. `arda build-index --organism all` rebuilds the shipped
indexes for your installed mmseqs version.
