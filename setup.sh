#!/usr/bin/env bash
# arda bootstrap — source of truth for a reproducible install.
#
# Steps:
#   1. Create/update the `arda` conda environment (python + mmseqs2 + toolchain).
#   2. Download the latest IgBLAST release into ./bin (gitignored).
#   3. pip install -e . (compiles the _markup C++ extension).
#
# Flags:
#   --no-conda    Skip conda env creation (use the already-active environment).
#   --build-db    After install, build the reference DB (arda build-db).
#   --tests       After install, run the fast test suites.
#
# Usage:
#   bash setup.sh [--no-conda] [--build-db] [--tests]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="arda"
USE_CONDA=1
DO_BUILD_DB=0
DO_TESTS=0

for arg in "$@"; do
  case "$arg" in
    --no-conda) USE_CONDA=0 ;;
    --build-db) DO_BUILD_DB=1 ;;
    --tests)    DO_TESTS=1 ;;
    *) echo "Unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m[arda]\033[0m %s\n' "$*"; }

# --- 1. conda environment --------------------------------------------------
if [[ "$USE_CONDA" -eq 1 ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found on PATH; install miniconda/anaconda or pass --no-conda." >&2
    exit 1
  fi
  if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    log "conda env '$ENV_NAME' exists — updating from environment.yml"
    conda env update -n "$ENV_NAME" -f "$ROOT/environment.yml" --prune
  else
    log "creating conda env '$ENV_NAME' from environment.yml"
    conda env create -f "$ROOT/environment.yml"
  fi
  # Run subsequent python/pip inside the env.
  PY="conda run -n $ENV_NAME python"
else
  PY="python"
fi

# --- 2. IgBLAST release ----------------------------------------------------
"$PY" "$ROOT/scripts/fetch_igblast.py" --dest "$ROOT/bin"

# --- 3. editable install (builds the C++ extension) ------------------------
log "pip install -e . (builds _markup)"
$PY -m pip install -e "$ROOT"

# --- 3b. mmseqs2 (no-conda only) -------------------------------------------
# The conda env provides mmseqs2; without conda, fetch a static binary into bin/
# (arda also auto-fetches lazily on first use, so this is just eager).
if [[ "$USE_CONDA" -eq 0 ]]; then
  if ! command -v mmseqs >/dev/null 2>&1; then
    log "fetching static mmseqs2 binary into ./bin"
    $PY "$ROOT/scripts/fetch_mmseqs.py" --dest "$ROOT/bin" || true
  fi
fi

# --- 4. verification -------------------------------------------------------
log "verifying toolchain"
if [[ "$USE_CONDA" -eq 1 ]]; then
  conda run -n "$ENV_NAME" mmseqs version || true
fi
"$ROOT/bin/igblastn" -version | head -1 || true
$PY -c "import arda._markup as m; print('arda._markup', m.__version__)"
$PY -c "import arda; print('arda', arda.__version__)"

# --- optional follow-ups ---------------------------------------------------
if [[ "$DO_BUILD_DB" -eq 1 ]]; then
  log "building reference database"
  $PY -m arda.cli build-db --organism all
fi
if [[ "$DO_TESTS" -eq 1 ]]; then
  log "running fast tests"
  $PY -m pytest "$ROOT/tests/unit" "$ROOT/tests/synthetic" -q || true
fi

log "done."
