"""`segments.fasta` is GENERATED, not shipped — so a plain pip install has none.

The auto-fetched reference tarball does not carry it, so before this every `--two-pass` /
`--fast-segments` run on a pip install degraded to the one-pass search behind a single log line.
The flagship configuration was unreachable out of the box, and the failure was invisible: correct
output, exit 0, none of the speed.

⛔ The two questions "is this file STALE" and "does this file EXIST" need DIFFERENT done-predicates,
which is why there are two functions. `_has_jc_targets` is false for a missing file, so reusing the
stale-format predicate would make a missing file read as *already regenerated* and the lock would
skip the build — silently, in the same direction as the bug being fixed.
"""

from __future__ import annotations

from pathlib import Path

from arda.annotate.mapper import _has_jc_targets


def test_has_jc_targets_is_total_and_false_for_a_missing_file(tmp_path):
    """It used to raise FileNotFoundError, which would blow up inside the build lock."""
    assert _has_jc_targets(tmp_path / "nope.fasta") is False


def test_has_jc_targets_detects_the_pre_2_8_0_format(tmp_path):
    old = tmp_path / "old.fasta"
    old.write_text(">V|TRBV1*01\nACGT\n>JC|TRB_5\nACGT\n")
    assert _has_jc_targets(old) is True
    new = tmp_path / "new.fasta"
    new.write_text(">V|TRBV1*01\nACGT\n>J|TRBJ1-1*01\nACGT\n>C|TRBC1*01\nACGT\n")
    assert _has_jc_targets(new) is False


def test_generate_and_regenerate_are_separate_functions_with_different_predicates():
    """Pinned because collapsing them is the silent-failure direction."""
    import inspect

    from arda.annotate import mapper

    def body(fn):
        # Strip the docstring: `_generate_segments` NAMES `_has_jc_targets` in prose precisely to
        # explain why it must not call it, and an assertion that cannot tell the two apart would
        # fail on the explanation rather than on the code.
        src = inspect.getsource(fn)
        head, _, rest = src.partition('"""')
        _, _, tail = rest.partition('"""')
        return head + tail

    gen, regen = body(mapper._generate_segments), body(mapper._regenerate_segments)
    assert "fasta.exists()" in gen, "generation must gate on EXISTENCE"
    assert "_has_jc_targets" in regen, "regeneration must gate on the stale FORMAT"
    assert "_has_jc_targets" not in gen, (
        "generation must not reuse the stale-format predicate: a missing file has no JC| targets, "
        "so it would read as already-done and the build would be skipped")


def test_the_map_path_generates_a_missing_segments_reference():
    """`_cached_segment_db` must build the file, not return None and fall back."""
    import inspect

    from arda.annotate import mapper

    src = inspect.getsource(mapper._cached_segment_db)
    assert "_generate_segments" in src, "a missing segments.fasta must be generated, not tolerated"
    # And it must still degrade rather than raise if that fails.
    assert "return None" in src
