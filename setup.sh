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
uv venv "$ROOT/.venv"
# shellcheck disable=SC1091
source "$ROOT/.venv/bin/activate"
uv pip install pybind11 scikit-build-core ninja
uv pip install -e "$ROOT" --no-build-isolation

# --- 2. IgBLAST release (offline DB build only) ----------------------------
python "$ROOT/scripts/fetch_igblast.py" --dest "$ROOT/bin"

# --- 3. MMseqs2 static binary (runtime) ------------------------------------
# arda also auto-fetches lazily on first use, so this is just eager.
command -v mmseqs >/dev/null 2>&1 || \
  python "$ROOT/scripts/fetch_mmseqs.py" --dest "$ROOT/bin" || true

# --- 4. verification -------------------------------------------------------
log "verifying"
python -c "import arda; print('arda', arda.__version__)"
python -c "from arda.mmseqs import version; print('mmseqs', version())"
"$ROOT/bin/igblastn" -version | head -1 || true

# --- optional follow-ups ---------------------------------------------------
if [[ "$DO_BUILD_DB" -eq 1 ]]; then
  log "building reference database"
  arda build-db --organism all
fi
if [[ "$DO_TESTS" -eq 1 ]]; then
  log "running fast tests"
  python -m pytest "$ROOT/tests/unit" "$ROOT/tests/synthetic" -q || true
fi

log "done. activate with: source .venv/bin/activate"
