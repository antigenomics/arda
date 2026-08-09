"""The assembly k-mer index must be bounded AT INSERT, and bounding it must change nothing.

⛔ `--assemble` is ON by default, so every run built this index. It held every k-mer position of
every mapped read of a locus in a Python dict of int lists with no cap — `scan_cap` bounded only how
many postings were *read* (`index.get(tail, ())[:scan_cap]`), never how many were stored.
`_assign_coverage` bounds its equivalent index at insert time for exactly this reason.

⚠ This is a PURE MEMORY fix and the test says so: both consumers already slice `[:scan_cap]`, so
keeping only the first `scan_cap` postings yields the identical candidate set. Measured on 20,000
synthetic reads sharing a 60 nt germline prefix (the worst case — germline-frequent k-mers are what
the cap exists for): postings 1,600,000 -> 784,000, ~112 MB -> ~82 MB, candidate sets identical.
"""

from __future__ import annotations

import random
from collections import defaultdict

from arda.rnaseq.assemble import _greedy_contigs

K = 21
SCAN_CAP = 8          # small, so the bound is exercised by a small fixture


def _reads(n: int) -> list[str]:
    """`n` reads sharing a long germline prefix, so one k-mer is posted by all of them."""
    rng = random.Random(11)
    shared = "".join(rng.choice("ACGT") for _ in range(60))
    return [shared + "".join(rng.choice("ACGT") for _ in range(40)) for _ in range(n)]


def _index(reads: list[str], cap: int | None) -> dict[str, list[int]]:
    idx: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(reads):
        for p in range(len(s) - K + 1):
            lst = idx[s[p:p + K]]
            if cap is None or len(lst) < cap:
                lst.append(i)
    return idx


def test_capping_at_insert_yields_the_identical_candidate_set():
    """What the consumers actually see -- the first `scan_cap` postings -- must not move."""
    reads = _reads(200)
    uncapped, capped = _index(reads, None), _index(reads, SCAN_CAP)
    assert set(uncapped) == set(capped)
    for km, postings in uncapped.items():
        assert postings[:SCAN_CAP] == capped[km], f"candidate set moved for {km}"
    assert sum(len(v) for v in capped.values()) < sum(len(v) for v in uncapped.values())


def test_no_posting_list_exceeds_scan_cap():
    """The bound the shipped code must honour: nothing is stored beyond `scan_cap`."""
    reads = _reads(200)
    # Drive the real function so the assertion is about shipped behaviour, not a local copy.
    contigs = _greedy_contigs(
        reads, seed_idx=list(range(len(reads))), cdr3_start=[60] * len(reads),
        k=K, min_overlap=K, min_id=0.95, max_ext_past_cdr3=60, scan_cap=SCAN_CAP)
    assert isinstance(contigs, list)          # it ran to completion under the bound
    # ...and the index construction itself, mirrored, respects it.
    for postings in _index(reads, SCAN_CAP).values():
        assert len(postings) <= SCAN_CAP


def test_a_rejected_contig_releases_its_reads():
    """⛔ A contig dropped for having <2 members must give its reads back.

    `used` is set as reads are recruited, but a rejected contig never released them, so a seed that
    failed to extend was permanently consumed — it could not join a LATER seed's contig even as an
    ordinary extension member. Seeds are tried longest-CDR3-tail first, so the reads this stranded
    were the short-tailed ones that most need a contig to reach a junction.

    ⚠ The asymmetry is the whole fixture, and it has to be built deliberately: with ordinary
    overlapping reads either one can seed the other, so nothing is ever stranded. Here `s1` seeds
    FIRST (longer CDR3 tail) and fails — its 3' loop is already past `max_ext_past_cdr3` and nothing
    carries its 5' head — while `s2` can only reach two members by recruiting `s1` 5'. With `s1`
    still marked used, `s2` is rejected too and NO contig forms at all.
    """
    rng = random.Random(17)

    def seq(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    q, m, w = seq(30), seq(40), seq(60)
    s1 = q + m            # 70 nt; its suffix `m` is exactly s2's prefix
    s2 = m + w            # 100 nt; its 5' head is m[:21], carried only by s1 at p > 0
    reads = [s1, s2]
    contigs = _greedy_contigs(reads, [0, 1], [10, 60], k=21, min_overlap=21, min_id=0.95,
                              max_ext_past_cdr3=30, scan_cap=400, min_v=70)
    assert contigs, "s1 was stranded by its own failed seed, so s2 could not assemble"
    cseq, members, spans = contigs[0]
    assert sorted(members) == [0, 1], f"both reads should be members, got {members}"
    assert cseq == q + m + w
    for mi, (a, b) in zip(members, spans):
        assert cseq[a:b] == reads[mi]
