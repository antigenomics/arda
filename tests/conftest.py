"""Shared pytest fixtures and capability detection.

Unit tests run anywhere. Integration tests need both the MMseqs2 binary and a
built reference DB; they skip cleanly when either is unavailable.
"""

from __future__ import annotations

import pytest

from arda import mmseqs
from arda.paths import vdj_dir


def mmseqs_available() -> bool:
    try:
        mmseqs.mmseqs_binary()
        return True
    except Exception:
        return False


def human_db_available() -> bool:
    return (vdj_dir("human") / "alleles.fasta").exists()


requires_mmseqs = pytest.mark.skipif(
    not mmseqs_available(), reason="mmseqs binary not found (set $ARDA_MMSEQS or activate env)"
)
requires_human_db = pytest.mark.skipif(
    not human_db_available(), reason="human reference DB not built (run `arda build-db`)"
)


@pytest.fixture(scope="session")
def human_scaffolds():
    """A handful of (id, nt_seq) reference scaffolds for integration tests."""
    from pathlib import Path
    from arda.refbuild.imgt import read_fasta

    path = vdj_dir("human") / "alleles.fasta"
    if not path.exists():
        pytest.skip("human DB not built")
    return read_fasta(Path(path))
