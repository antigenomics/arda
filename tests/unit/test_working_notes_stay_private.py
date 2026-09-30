"""Working notes must never enter the public repository or a published artifact.

Several markdown files in this tree exist to hold internal detail -- benchmark repository paths,
private dataset names, measured dead ends, in-flight decisions. ``SOURCES.md`` is the clearest case:
its *purpose* is provenance that is not for publication. ``CLAUDE.md`` is how to work in the repo,
and ``CHANGELOG_local.md`` / ``ROADMAP_local.md`` are the detailed versions of the two documents
whose published forms are written for a USER.

They are kept out by two independent mechanisms, and this pins both:

* ``.gitignore`` keeps them untracked, so ``git add -A`` cannot sweep one in;
* ``pyproject.toml``'s sdist exclude keeps them out of the published package.

Neither is self-announcing when it breaks -- a newly tracked ``NOTES.md`` looks like any other file
in a diff -- which is why this is a test and not a convention. (Ported from vdjtools, which has the
same policy and learned it the same way.)
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

#: Filenames whose content is internal by design. Matched against the basename, case-insensitively.
PRIVATE = re.compile(
    r"^(SOURCES|CLAUDE|TODO|NOTES|ROADMAP_local|CHANGELOG_local|STATUS|AGENTS|NULLS|PLAN)\.md$"
    r"|^ISSUES.*\.md$",
    re.IGNORECASE,
)

#: Every working-note name, as the sdist exclude and the gitignore must both spell them.
NAMES = ("SOURCES.md", "CLAUDE.md", "CHANGELOG_local.md", "ROADMAP_local.md",
         "TODO.md", "NOTES.md", "STATUS.md", "ISSUES.md", "ISSUES_ext.md")


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True,
                          check=True).stdout


def test_no_working_note_is_tracked_by_git():
    tracked = [p for p in _git("ls-files").split("\n") if p]
    leaked = sorted(p for p in tracked if PRIVATE.match(Path(p).name))
    assert not leaked, (
        f"working notes are tracked and would be published: {leaked}. "
        f"Run `git rm --cached <file>` and confirm .gitignore covers it."
    )


@pytest.mark.parametrize("name", NAMES)
def test_each_working_note_name_is_gitignored(name):
    """Ignored by NAME, not merely absent -- an absent file that is not ignored is a future leak."""
    r = subprocess.run(["git", "check-ignore", "-q", name], cwd=ROOT)
    assert r.returncode == 0, (
        f"{name} is not gitignored. It may not exist today, but nothing stops it being created "
        f"and committed tomorrow."
    )


def test_the_sdist_exclude_names_every_working_note():
    """Asserted against the file's text, not a parsed key path.

    scikit-build-core and hatchling spell the exclude table differently, and a test that follows
    one of them silently stops checking if the backend changes.
    """
    text = (ROOT / "pyproject.toml").read_text()
    missing = [n for n in NAMES if f'"{n}"' not in text and f'"/{n}"' not in text]
    assert not missing, (
        f"not excluded from the sdist: {missing}. gitignore does not reach the build backend for "
        f"a file that is present but untracked."
    )


def test_the_published_documents_are_not_caught_by_the_matcher():
    """The guard on the guard: it must catch what it exists for and leave the public files alone."""
    assert PRIVATE.match("SOURCES.md") and PRIVATE.match("claude.md")
    assert PRIVATE.match("CHANGELOG_local.md") and PRIVATE.match("ROADMAP_local.md")
    assert not PRIVATE.match("README.md")
    assert not PRIVATE.match("CHANGELOG.md"), "the slim changelog is for a user, and ships"
    assert not PRIVATE.match("ROADMAP.md"), "so is the slim roadmap"
    assert not PRIVATE.match("SKILL.md")
