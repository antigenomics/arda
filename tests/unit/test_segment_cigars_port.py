"""The C++ `segment_cigars` must agree with the Python one on every input, not just the easy ones.

`segment_cigars` walks the alignment column by column and was 63 % of `transfer_hit` -- the term
that caps arda on receptor-rich libraries. It was ported to C++ for that reason, and a port of a
per-column loop is exactly the kind of change that is correct on the happy path and wrong on a
boundary: an insertion at a segment edge, a `J + C` scaffold where `t_vend == 0`, an alignment
that starts inside the N-pad and therefore belongs to no germline.

So this is differential, over randomised alignments that deliberately straddle those boundaries.
"""
from __future__ import annotations

import random

import pytest

from arda.annotate.cigar import _segment_cigars_py

pytest.importorskip("arda._markup", reason="native markup extension not built")

from arda._markup import segment_cigars as cpp  # noqa: E402


def _random_alignment(rnd: random.Random, n: int) -> tuple[str, str]:
    """An aligned pair with gaps on both sides, never a gap aligned to a gap."""
    q, t = [], []
    for _ in range(n):
        r = rnd.random()
        if r < 0.10:
            q.append("-"); t.append(rnd.choice("ACGT"))       # deletion
        elif r < 0.20:
            q.append(rnd.choice("ACGT")); t.append("-")       # insertion
        else:
            b = rnd.choice("ACGT")
            q.append(b); t.append(b if rnd.random() < 0.9 else rnd.choice("ACGT"))
    return "".join(q), "".join(t)


@pytest.mark.parametrize("seed", range(40))
def test_cpp_matches_python_on_random_alignments(seed):
    rnd = random.Random(seed)
    qaln, taln = _random_alignment(rnd, rnd.randint(1, 160))
    qstart = rnd.randint(1, 40)
    tstart = rnd.randint(1, 300)
    qlen = qstart + len(qaln.replace("-", "")) + rnd.randint(0, 50)
    # A V+pad+J scaffold: V ends, an N-pad, then J, then C.
    t_vend = rnd.choice([0, rnd.randint(50, 250)])
    t_jstart = t_vend + rnd.randint(1, 30) if t_vend else rnd.randint(1, 60)
    t_vjend = t_jstart + rnd.randint(10, 60)
    args = (qaln, taln, qstart, tstart, qlen, t_vend, t_jstart, t_vjend)
    assert cpp(*args) == _segment_cigars_py(*args), f"diverged at seed {seed}: {args}"


def test_a_jc_scaffold_has_no_v(seed=0):
    """`t_vend == 0` means the segment is absent — a J+C scaffold. Treating 0 as a real boundary
    would put every position at or below it into V."""
    qaln = taln = "ACGT" * 10
    args = (qaln, taln, 1, 5, 40, 0, 3, 30)
    out = cpp(*args)
    assert out == _segment_cigars_py(*args)
    assert "v_cigar" not in out


def test_an_alignment_wholly_inside_the_np_pad_yields_nothing():
    """Between the V end and the J start there is no germline, so no segment gets a body."""
    qaln = taln = "ACGTACGTAC"
    args = (qaln, taln, 1, 105, 10, 100, 200, 260)
    assert cpp(*args) == _segment_cigars_py(*args) == {}


@pytest.mark.parametrize("seed", range(30))
def test_aln_identity_cpp_matches_python(seed):
    """Second per-column loop moved to C++. The boundary that matters is "nothing covered": the
    Python returns "" and the extension returns -1.0, and conflating that with a real identity of
    0.0 would silently report a perfectly mismatched alignment as an empty field."""
    from arda.annotate.transfer import _aln_identity, _aln_identity_py

    rnd = random.Random(1000 + seed)
    qaln, taln = _random_alignment(rnd, rnd.randint(1, 160))
    tstart = rnd.randint(1, 100)
    t_lo = rnd.randint(1, 200)
    t_hi = t_lo + rnd.randint(0, 100)
    assert _aln_identity(qaln, taln, tstart, t_lo, t_hi) == \
        _aln_identity_py(qaln, taln, tstart, t_lo, t_hi)


def test_aln_identity_distinguishes_uncovered_from_zero_identity():
    from arda.annotate.transfer import _aln_identity, _aln_identity_py

    # Range entirely past the alignment: nothing covered -> "" from both paths.
    assert _aln_identity("ACGT", "ACGT", 1, 900, 999) == "" == \
        _aln_identity_py("ACGT", "ACGT", 1, 900, 999)
    # Covered but every base mismatched -> a real 0.0, not "".
    assert _aln_identity("AAAA", "TTTT", 1, 1, 4) == 0.0
    assert _aln_identity_py("AAAA", "TTTT", 1, 1, 4) == 0.0
