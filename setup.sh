#!/usr/bin/env bash
# arda bootstrap (uv) — reproducible dev install, no conda.
#
# Steps:
#   1. uv venv .venv + editable install (compiles the _markup C++ extension).
#   2. Download the latest IgBLAST release into ./bin (gitignored; DB build only).
#   3. Fetch a static MMseqs2 binary into ./bin unless one is already on PATH.
#
# Conda is used ONLY by the Nextflow integration (integrations/nextflow/arda, the
# Gamaleya/ISP pipeline), which ships its own environment.yml + Dockerfile.
#
# Usage:
#   bash setup.sh [--build-db] [--tests]
set -euo pipefail

# Script dir, portable to both `bash setup.sh` and `zsh setup.sh` (BASH_SOURCE is
# bash-only; $0 is the script path when executed, in both shells).
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
DO_BUILD_DB=0
DO_TESTS=0
for arg in "$@"; do
  case "$arg" in
    --build-db) DO_BUILD_DB=1 ;;
    --tests)    DO_TESTS=1 ;;
    *) echo "Unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m[arda]\033[0m %s\n' "$*"; }

command -v uv >/dev/null 2>&1 || {
  echo "uv not found — install it: curl -LsSf https://astral.sh/uv/install.sh | sh" >&2
  exit 1
}

# --- 1. venv + editable install --------------------------------------------
# Build deps go in the venv so the scikit-build editable on-import rebuild can
# find pybind11 (hence --no-build-isolation).
log "creating .venv and installing arda (uv)"
# ⛔ A stale build/ dir is not harmless. scikit-build-core caches CMake's configuration, including
# the ABSOLUTE PATH of the interpreter it configured against; if that venv is gone (a previous
# checkout, a deleted conda env) every later on-import rebuild fails with "Could NOT find Python"
# and arda silently falls back to the pure-Python markup path. Start clean.
rm -rf "$ROOT/build"
uv venv "$ROOT/.venv"
# shellcheck disable=SC1091
source "$ROOT/.venv/bin/activate"
# Pinned to the same range as pyproject.toml's build-system: `PYBIND11_MODULE` changed to
# multi-phase init in 3.0.0, so an unbounded install builds against whatever PyPI serves today.
uv pip install 'pybind11>=3.0.2,<4' scikit-build-core ninja
# `.[test,dev]` -- not a bare install. `--tests` used to run `pytest` that was never installed,
# print "No module named pytest", and be swallowed by `|| true` while the script reported success.
uv pip install -e "$ROOT[test,dev]" --no-build-isolation

# --- 2. IgBLAST release (offline DB build only) ----------------------------
python "$ROOT/scripts/fetch_igblast.py" --dest "$ROOT/bin"

# --- 3. MMseqs2 static binary (runtime) ------------------------------------
# arda also auto-fetches lazily on first use, so this is just eager.
command -v mmseqs >/dev/null 2>&1 || \
  python "$ROOT/scripts/fetch_mmseqs.py" --dest "$ROOT/bin" || true

# --- 4. verification -------------------------------------------------------
# ⛔ `import arda` is NOT a check that the build worked. arda falls back to a pure-Python markup
# path when `_markup` is missing, so a failed C++ build looks like a successful install and shows
# up later as a silent ~2x slowdown. Assert the extension, and assert the CLI surface on the CLI --
# a deploy into the wrong environment prints a correct version and still lacks the commands.
log "verifying"
python -c "import arda; print('arda', arda.__version__)"
python -c "import arda._markup as m; print('_markup', m.__version__)"
python -c "import arda._segmap, arda._denoise" \
  && log "_segmap + _denoise built" \
  || { echo "FAIL: a C++ extension did not build" >&2; exit 1; }
python -c "from arda.mmseqs import version; print('mmseqs', version())"
# Resolve each command rather than grepping `--help`: typer renders that through rich, in a box,
# wrapped to $COLUMNS, so a grep over it tests the renderer as much as the CLI.
for cmd in rnaseq amplicon singlecell map correct assemble shm cluster annotate; do
  arda "$cmd" --help >/dev/null 2>&1 \
    || { echo "FAIL: 'arda ${cmd}' does not resolve" >&2; exit 1; }
done
log "CLI surface OK (modes: rnaseq/amplicon/singlecell; stages: map/correct/assemble/shm)"
"$ROOT/bin/igblastn" -version | head -1 || true

# --- optional follow-ups ---------------------------------------------------
if [[ "$DO_BUILD_DB" -eq 1 ]]; then
  log "building reference database"
  arda build-db --organism all
fi
if [[ "$DO_TESTS" -eq 1 ]]; then
  # ⛔ NOT `|| true`. A swallowed test failure under a script that then prints "done" is worse
  # than no test run at all -- that is exactly what this flag did before.
  log "running fast tests"
  python -m pytest "$ROOT/tests/unit" "$ROOT/tests/synthetic" -q
fi

log "done. activate with: source .venv/bin/activate"
