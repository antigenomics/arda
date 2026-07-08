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


def imgt_reference_available() -> bool:
    """True when the IMGT V-QUEST reference has been downloaded.

    Reference-*build* tests (e.g. those that call ``build_jc_scaffolds``, which loads J alleles
    from IMGT) need it; a plain checkout / CI has the committed ``database/`` but not the IMGT
    download, so those tests must skip there rather than fail. The committed reference is enough
    for *annotation* tests -- this gate is only for tests that re-run the build.
    """
    from arda.refbuild import imgt
    root = imgt.reference_dir()
    return root.is_dir() and any(root.iterdir())


requires_mmseqs = pytest.mark.skipif(
    not mmseqs_available(), reason="mmseqs binary not found (set $ARDA_MMSEQS or activate env)"
)
requires_human_db = pytest.mark.skipif(
    not human_db_available(), reason="human reference DB not built (run `arda build-db`)"
)
requires_imgt = pytest.mark.skipif(
    not imgt_reference_available(),
    reason="IMGT V-QUEST reference not downloaded (build-pipeline test; run `arda build-db`)",
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
