"""Reference-free per-cell contig assembly.

A droplet V(D)J library gives one barcode per cell and, under it, many UMI-tagged molecules.
An upstream tool (migec, fgbio, UMI-tools) collapses each molecule's reads into one consensus;
what arrives here is therefore a *cell* worth of short, error-corrected fragments. Reads of one
molecule are co-terminal in this chemistry, so a molecule's consensus is a pile at one position
and covers one window of the transcript however deep it is -- but different molecules start at
different positions, so the cell's molecules tile the whole transcript. Measured on 10x's
``sc5p_v2_hs_PBMC_1k`` VDJ-T: every one of Cell Ranger's 943 filtered contigs is >= 99.9%
covered by the 25-mers of its own cell's molecules. The full-length contig is already in the
data; assembling it is what this module does.

The method is overlap-layout with no germline reference anywhere: index every 25-mer of a cell's
molecules, propose an offset wherever two molecules share one, *verify* the implied overlap, and
union-find the survivors into components. Each component's weighted column consensus is a contig.

Two things are load-bearing and both were measured by removing them (CDR3 nt recovered exactly
inside a contig, 244 Cell Ranger chains over the first 120 cells):

===========================================  ======
trimmed + depth-weighted                     0.988
trimmed, every molecule weighted 1           0.967
untrimmed + depth-weighted                   0.898
untrimmed, every molecule weighted 1         0.873
===========================================  ======

* **Never: verify the overlap before joining, and trim the read-through adapter first.** Every
  molecule that ran off the end of its insert carries the same Illumina adapter, so an adapter
  suffix is a shared k-mer *and* a 100%-identical overlap: verification alone cannot refuse it.
  Before verification existed, one component swallowed a cell's TRA and TRB alike. With it, the
  damage is quieter and still 9 points -- the largest component of a cell grows from 70 to 78
  molecules (median over 40 cells) and the adapter columns corrupt the consensus where it does
  merge. Both failures have the same cure, and neither substitutes for the other.
* **Never: weight the consensus by the molecule's read depth.** A one-read molecule is one read;
  giving it the same column vote as an 800-read consensus is what the 0.967 row costs.

Assembly is per cell and cells are independent, so nothing here scales with the library.
Note: it is single-threaded Python at ~54 ms per cell (479 cells in 26 s). A process pool over
cells is the obvious upgrade and was tried and reverted -- see the comment in
:func:`assemble_cells`; on macOS `spawn` re-imports the caller's ``__main__`` in every worker and
`forkserver` will not start, so a pool here crashes library callers to save minutes on a
10k-cell library. The upgrade path is the assembler in C++, not a pool.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path

from .annotate.io import open_text
from .cell import make_parser

__all__ = [
    "Contig", "CellAssemblyReport", "assemble_cell", "assemble_cells", "read_molecules",
    "trim", "pair_chains", "cell_summary", "cell_rank", "threshold_sweep", "HEAVY", "LIGHT",
    "read_called_cells", "read_reference_calls", "read_records", "run", "find_knee",
]

#: Seed length. Long enough that a chance match inside one cell's ~100 kb of sequence is
#: negligible, short enough that a 90 nt fragment still carries 66 of them.
K = 25

#: Illumina TruSeq/Nextera read-through stub. What follows it in a consensus is adapter, not
#: transcript, and it is identical in every molecule that reached it.
ADAPTER = "AGATCGGAAGAGC"

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_HOMOPOLYMER = re.compile(r"(A{12,}|C{12,}|G{12,}|T{12,}).*$")


def _rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def trim(seq: str, *, adapter: str = ADAPTER, min_stub: int = 8) -> str:
    """Cut a read-through adapter and any homopolymer tail off a molecule consensus.

    Both orientations are tried, and a partial stub counts: a consensus that ended eight bases
    into the adapter still shares those eight bases with every other such consensus.
    """
    for probe in (adapter, _rc(adapter)):
        for length in range(len(probe), min_stub - 1, -1):
            cut = seq.find(probe[:length])
            if cut >= 0:
                seq = seq[:cut]
                break
    return _HOMOPOLYMER.sub("", seq)


class _OffsetUnion:
    """Union-find whose members carry an integer offset in their root's coordinates."""

    __slots__ = ("parent", "offset")

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.offset = [0] * n

    def find(self, i: int) -> int:
        path = []
        while self.parent[i] != i:
            path.append(i)
            i = self.parent[i]
        for j in reversed(path):
            self.offset[j] += self.offset[self.parent[j]]
            self.parent[j] = i
        return i

    def join(self, a: int, b: int, offset: int) -> bool:
        """Place ``b`` at ``offset`` relative to ``a``. False if they already share a root."""
        root_a, root_b = self.find(a), self.find(b)
        if root_a == root_b:
            return False
        self.parent[root_b] = root_a
        self.offset[root_b] = self.offset[a] + offset - self.offset[b]
        return True


@dataclass(frozen=True, slots=True)
class Contig:
    """One overlap component of a cell, and what supports it."""

    cell: str
    index: int
    sequence: str
    molecules: int
    reads: int

    @property
    def contig_id(self) -> str:
        """``<cell>_contig_<n>``, which :func:`arda.cell.parse` reads back as the cellranger
        dialect -- so an annotation run over these contigs recovers ``cell_id`` for free."""
        return f"{self.cell}_contig_{self.index}"


@dataclass
class CellAssemblyReport:
    cells: int = 0
    molecules_in: int = 0
    molecules_placed: int = 0
    contigs: int = 0
    contigs_dropped_small: int = 0
    contig_lengths: list[int] = field(default_factory=list)
    #: cell -> molecules READ for it, before assembly. Never derive this by summing contig
    #: membership: a molecule that joined no contig is missing from that sum, and one shared
    #: between two phased haplotypes is counted twice -- so a doublet reads as the deepest cell
    #: in the library and the knee plot's y axis stops being unique molecules.
    molecules_by_cell: dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> dict:
        lengths = sorted(self.contig_lengths, reverse=True)
        total = sum(lengths)
        n50 = 0
        acc = 0
        for length in lengths:
            acc += length
            if acc >= total / 2:
                n50 = length
                break
        return {
            "cells": self.cells,
            "molecules_in": self.molecules_in,
            "molecules_placed": self.molecules_placed,
            "contigs": self.contigs,
            "contigs_dropped_small": self.contigs_dropped_small,
            "contig_n50": n50,
            "contig_max": lengths[0] if lengths else 0,
        }


def _overlap_agrees(x: str, y: str, offset: int, min_overlap: int, min_identity: float) -> bool:
    """Does placing ``y`` at ``offset`` in ``x``'s coordinates survive the overlap it implies?"""
    lo = max(0, offset)
    hi = min(len(x), offset + len(y))
    span = hi - lo
    if span < min_overlap:
        return False
    same = 0
    for p in range(lo, hi):
        if x[p] == y[p - offset]:
            same += 1
    return same >= min_identity * span


def _informative_columns(
    local: list[str],
    local_w: list[int],
    members: list[tuple[int, int]],
    lo: int,
    width: int,
    min_minor_fraction: float,
    min_minor_molecules: int,
) -> list[int]:
    """Columns where a second base holds a real share of the layout, not an error's share."""
    weight: list[dict[str, int]] = [defaultdict(int) for _ in range(width)]
    seen: list[dict[str, int]] = [defaultdict(int) for _ in range(width)]
    for i, off in members:
        seq, w, base = local[i], local_w[i], off - lo
        for p, char in enumerate(seq):
            weight[base + p][char] += w
            seen[base + p][char] += 1
    out = []
    for column in range(width):
        counts = sorted(weight[column].items(), key=lambda kv: -kv[1])
        if len(counts) < 2:
            continue
        minor, minor_weight = counts[1]
        total = sum(w for _, w in counts)
        if minor_weight >= min_minor_fraction * total and seen[column][minor] >= min_minor_molecules:
            out.append(column)
    return out


def _phase_two(
    local: list[str],
    local_w: list[int],
    members: list[tuple[int, int]],
    calls: dict[int, dict[int, str]],
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]] | None:
    """Split a layout into two haplotypes by iterative reassignment, or ``None`` if it is one.

    Never: **first-fit greedy phasing is not good enough here.** Assigning each molecule to the
    first group it does not contradict lets a molecule covering two informative columns fix a
    profile the deep molecules then disagree with, and one component of two TRB transcripts
    splintered into four contigs, none of which carried a callable junction. Two seeds refined
    against the whole assignment converge instead: same input, two contigs, both junctions.

    A molecule that covers no informative column joins BOTH haplotypes -- the shared constant
    region belongs to both chains, and withholding it leaves each one a stub around its CDR3.
    """
    def profile_of(group: list[tuple[int, int]]) -> dict[int, str]:
        weight: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
        for i, _ in group:
            for column, base in calls[i].items():
                weight[column][base] += local_w[i]
        return {c: max(bases.items(), key=lambda kv: kv[1])[0] for c, bases in weight.items()}

    def conflicts(mine: dict[int, str], profile: dict[int, str]) -> int:
        return sum(1 for c, base in mine.items() if profile.get(c, base) != base)

    informative = [m for m in members if calls[m[0]]]
    if len(informative) < 2:
        return None
    seed_a = max(informative, key=lambda m: (len(calls[m[0]]), local_w[m[0]]))
    left = profile_of([seed_a])
    seed_b = max(informative, key=lambda m: (conflicts(calls[m[0]], left), local_w[m[0]]))
    if conflicts(calls[seed_b[0]], left) == 0:
        return None
    right = profile_of([seed_b])

    group_a: list[tuple[int, int]] = []
    group_b: list[tuple[int, int]] = []
    for _ in range(4):
        group_a, group_b, shared = [], [], []
        for i, off in members:
            mine = calls[i]
            if not mine:
                shared.append((i, off))
                continue
            left_conflicts, right_conflicts = conflicts(mine, left), conflicts(mine, right)
            if left_conflicts < right_conflicts:
                group_a.append((i, off))
            elif right_conflicts < left_conflicts:
                group_b.append((i, off))
            else:
                shared.append((i, off))
        if not group_a or not group_b:
            return None
        left, right = profile_of(group_a), profile_of(group_b)
    return group_a + shared, group_b + shared


def _split_members(
    local: list[str],
    local_w: list[int],
    members: list[tuple[int, int]],
    *,
    min_molecules: int,
    min_minor_fraction: float,
    min_minor_molecules: int,
    min_informative_columns: int,
    depth: int = 1,
) -> list[list[tuple[int, int]]]:
    """Phase a component into haplotypes, two ways at a time, up to ``depth`` levels.

    Never: ``depth`` defaults to ONE. A second level split each haplotype of a synthetic doublet
    again into two near-identical halves -- the molecules shared between haplotypes carry their
    own disagreement into the sub-layout, which clears the gate a second time. The junctions
    survived it, so nothing failed loudly; what it cost was four contigs where two are right and
    a doubled molecule count on each. Raise it only with a triplet in front of you.
    """
    lo = min(off for _, off in members)
    width = max(off + len(local[i]) for i, off in members) - lo
    columns = _informative_columns(local, local_w, members, lo, width,
                                   min_minor_fraction, min_minor_molecules)
    if len(columns) < min_informative_columns:
        return [members]

    def call(i: int, off: int) -> dict[int, str]:
        base, seq = off - lo, local[i]
        return {c: seq[c - base] for c in columns if 0 <= c - base < len(seq)}

    pair = _phase_two(local, local_w, members, {i: call(i, off) for i, off in members})
    if pair is None or min(len(pair[0]), len(pair[1])) < min_molecules:
        return [members]
    out: list[list[tuple[int, int]]] = []
    for group in pair:
        if depth > 1 and len(group) >= 2 * min_molecules:
            out.extend(_split_members(
                local, local_w, group, min_molecules=min_molecules,
                min_minor_fraction=min_minor_fraction,
                min_minor_molecules=min_minor_molecules,
                min_informative_columns=min_informative_columns, depth=depth - 1))
        else:
            out.append(group)
    return out


def assemble_cell(
    cell: str,
    sequences: list[str],
    weights: list[int] | None = None,
    *,
    k: int = K,
    min_overlap: int = 30,
    min_identity: float = 0.95,
    min_molecules: int = 2,
    split: bool = True,
    min_minor_fraction: float = 0.2,
    min_minor_molecules: int = 2,
    min_informative_columns: int = 2,
) -> list[Contig]:
    """Assemble one cell's molecule consensuses into contigs, largest component first.

    ``weights`` is each molecule's read depth (``cD:i:`` in a migec consensus); ``None`` weighs
    every molecule 1, which costs about two points of junction recovery -- see the module table.

    Never: **a component must be phased before it is consensed, or a doublet is invisible.** Two
    chains of the SAME locus in one cell share their constant region and most of their V, so the
    overlap layout puts them in one component and the column consensus averages their junctions
    into a sequence that is neither. Measured on two TRB receptors tiled into one synthetic
    cell: one 918 nt contig, no junction called at all, both true junctions lost. With the split
    on, both come back. This is the same failure a de Bruijn graph has on shared germline, and it
    is why ``split`` defaults to on: the layout finds the component, the phasing finds the chain.

    Never: **the split gate is on the minor base's WEIGHT SHARE and its MOLECULE count, both.**
    A share alone promotes one deep molecule's error to an allele; a count alone promotes a
    column that a hundred molecules cover and two disagree on. ``min_informative_columns``
    requires the disagreement to co-segregate across positions rather than stand on one.
    """
    weights = [1] * len(sequences) if weights is None else weights
    usable = [i for i, s in enumerate(sequences) if len(s) >= k]
    if not usable:
        return []
    local = [sequences[i] for i in usable]
    local_w = [weights[i] for i in usable]

    seeds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for i, seq in enumerate(local):
        for p in range(len(seq) - k + 1):
            seeds[seq[p : p + k]].append((i, p))

    union = _OffsetUnion(len(local))
    for hits in seeds.values():
        if len(hits) < 2:
            continue
        # Chain every hit to the first one. A k-mer shared by m molecules gives m-1 candidate
        # joins, not m(m-1)/2: any pair it could relate is related through the anchor already.
        anchor, anchor_pos = hits[0]
        for other, other_pos in hits[1:]:
            if union.find(anchor) == union.find(other):
                continue
            offset = anchor_pos - other_pos
            if _overlap_agrees(local[anchor], local[other], offset, min_overlap, min_identity):
                union.join(anchor, other, offset)

    components: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i in range(len(local)):
        components[union.find(i)].append((i, union.offset[i]))

    def consense(members: list[tuple[int, int]]) -> str:
        lo = min(off for _, off in members)
        hi = max(off + len(local[i]) for i, off in members)
        columns: list[dict[str, int]] = [defaultdict(int) for _ in range(hi - lo)]
        for i, off in members:
            seq, weight, base = local[i], local_w[i], off - lo
            for p, char in enumerate(seq):
                columns[base + p][char] += weight
        return "".join(
            max(col.items(), key=lambda kv: kv[1])[0] if col else "N" for col in columns
        )

    contigs: list[Contig] = []
    for members in sorted(components.values(), key=len, reverse=True):
        if len(members) < min_molecules:
            continue
        parts = [members]
        if split and len(members) >= 2 * min_molecules:
            phased = _split_members(
                local, local_w, members, min_molecules=min_molecules,
                min_minor_fraction=min_minor_fraction,
                min_minor_molecules=min_minor_molecules,
                min_informative_columns=min_informative_columns)
            if len(phased) > 1:
                parts = [g for g in phased if len(g) >= min_molecules]
        for part in sorted(parts, key=len, reverse=True):
            contigs.append(
                Contig(
                    cell=cell,
                    index=len(contigs) + 1,
                    sequence=consense(part),
                    molecules=len(part),
                    reads=sum(local_w[i] for i, _ in part),
                )
            )
    return contigs


_DEPTH = re.compile(r"\bcD:i:(\d+)\b")


def read_records(path: str | Path):
    """Yield ``(name, sequence, depth)`` from a consensus FASTQ or FASTA.

    Never: the read depth is in the FASTQ *comment* (``cD:i:``), which the shared reader drops
    with the rest of the header. It is read here rather than by widening
    :func:`~arda.annotate.io.read_sequences`, whose callers are the bulk hot path. A record with
    no ``cD:i:`` weighs 1.
    """
    with open_text(path) as fh:
        first = fh.readline()
        fh.seek(0)
        if first[:1] == "@":
            for header, seq, _plus, _qual in zip(*[iter(fh)] * 4):
                if header[:1] != "@":
                    continue
                name, _, comment = header[1:].rstrip("\n").partition(" ")
                depth = _DEPTH.search(comment)
                yield name, seq.rstrip("\n").upper(), int(depth.group(1)) if depth else 1
        else:
            name, buf = None, []
            for line in fh:
                if line[:1] == ">":
                    if name is not None:
                        yield name, "".join(buf).upper(), 1
                    name = line[1:].rstrip("\n").split()[0] if line[1:].strip() else ""
                    buf = []
                elif line.strip():
                    buf.append(line.strip())
            if name is not None:
                yield name, "".join(buf).upper(), 1


def read_molecules(
    path: str | Path,
    *,
    cell_from: str = "auto",
    cell_regex: str | None = None,
    sniff_sample: int = 10000,
) -> dict[str, list[tuple[str, int]]]:
    """Group a consensus FASTQ/FASTA by cell barcode, keeping each record's read depth.

    Never: the dialect is sniffed over a RESERVOIR sample of the names, never the first N.
    ``migec assemble`` writes in barcode order, so the head of the file shares its leading bases
    and any statistic taken from it describes one corner of the barcode space.

    Never: TWO SAMPLES IN ONE FILE ARE REFUSED. The grouping key is the cell barcode alone, and a
    droplet barcode comes from a fixed whitelist -- concatenating two samples' consensus FASTQs
    merges every barcode they share, and the result reads as an ordinary doublet rather than as an
    error. migec writes one consensus file per sample; run this once per file.
    """
    records = list(read_records(path))
    if cell_regex is None and cell_from == "auto":
        step = max(1, len(records) // max(1, sniff_sample))
        parser = make_parser("auto", None, sample=[n for n, _, _ in records[::step]],
                             with_sample=True)
    else:
        parser = make_parser(cell_from, cell_regex, with_sample=True)
    out: dict[str, list[tuple[str, int]]] = defaultdict(list)
    samples: set[str] = set()
    for name, seq, depth in records:
        sample, cell = parser(name)
        if cell:
            if sample:
                samples.add(sample)
            out[cell].append((seq, depth))
    if records and not out:
        # `--cell-from auto` refuses when no dialect fits, but a NAMED dialect that parses these
        # names and finds no cell in them does not: a bulk migec consensus is `<sample>.<umi>`,
        # which `parse_migec` reads happily with `cell=None`, so the command succeeds having
        # assembled nothing at all.
        named = "regex" if cell_regex else repr(cell_from)
        raise ValueError(
            f"{path}: none of its {len(records):,} records carry a cell barcode under the "
            f"{named} dialect. A bulk migec consensus is named `<sample>.<umi>` and has no cell "
            "in it -- `arda rnaseq` is the command for that")
    if len(samples) > 1:
        shown = ", ".join(sorted(samples)[:4]) + (", ..." if len(samples) > 4 else "")
        raise ValueError(
            f"{path} holds {len(samples)} samples ({shown}) and the cell barcode is the only "
            "grouping key. Droplet barcodes come from a fixed whitelist, so two samples share "
            "barcodes by design and grouping on the barcode alone merges two cells into one. "
            "Run once per sample")
    return out


def _assemble_one(work: tuple[str, list[str], list[int]], **kw) -> list[Contig]:
    """One cell, at module level so a process pool can pickle it."""
    cell, seqs, weights = work
    return list(assemble_cell(cell, seqs, weights, **kw))


def assemble_cells(
    reads: str | Path,
    output: str | Path,
    *,
    cell_from: str = "auto",
    cell_regex: str | None = None,
    k: int = K,
    min_overlap: int = 30,
    min_identity: float = 0.95,
    min_molecules: int = 2,
    min_length: int = 0,
    trim_adapter: bool = True,
    split: bool = True,
    min_minor_fraction: float = 0.2,
    min_minor_molecules: int = 2,
    min_informative_columns: int = 2,
    cells: set[str] | None = None,
) -> tuple[list[Contig], CellAssemblyReport]:
    """Assemble every cell in a UMI-consensus FASTQ and write the contigs as FASTA.

    Args:
        reads: a per-molecule consensus FASTQ/FASTA whose record names carry the cell barcode --
            what ``migec assemble`` writes, or anything :mod:`arda.cell` can parse.
        output: contigs FASTA, named ``<cell>_contig_<n>``.
        cell_from: a :mod:`arda.cell` dialect, ``"auto"`` to sniff, or ``"regex"`` with
            ``cell_regex``.
        min_molecules: molecules a component needs before it is emitted. Never: 1 is not a
            useful setting -- on the 10x library, an extra chain supported by exactly one
            molecule is contamination 96% of the time (see :func:`pair_chains`).
        cells: restrict to these barcodes -- the called cells. Never: cell CALLING is the
            upstream tool's job (migec's ``refine`` writes ``<sample>.cells.tsv``); assembling
            every empty droplet costs the run and calls chains in ambient RNA.

    Returns:
        The contigs, and a :class:`CellAssemblyReport`.
    """
    grouped = read_molecules(reads, cell_from=cell_from, cell_regex=cell_regex)
    if cells is not None:
        grouped = {c: v for c, v in grouped.items() if c in cells}
    else:
        # Never: an unfiltered droplet library is REFUSED, not warned about. Cell calling is the
        # upstream tool's job, but a run without it does not fail -- it succeeds, having assembled
        # every ambient droplet that happened to carry two molecules, and those barcodes appear in
        # `.cells.tsv` and the chain table looking exactly like cells. On `sc5p_v2_hs_PBMC_1k` the
        # knee is at rank 376 of 136,032 observed barcodes: 99.7% of what would be assembled is
        # ambient RNA. The gate is the knee, which is refused when there is not one (a plate or
        # combinatorial library has no ambient tail, so `knee_rank` is 0 and this never fires).
        knee_rank, _ = find_knee([len(v) for v in grouped.values()])
        if knee_rank and len(grouped) > 2 * knee_rank:
            raise ValueError(
                f"{len(grouped):,} barcodes carry a molecule but the rank curve breaks at "
                f"{knee_rank:,} -- the rest is ambient RNA, and assembling it emits contigs and "
                "chain calls for empty droplets. Pass --cells with the called set (migec's "
                "`refine` writes `<sample>.cells.tsv`, whose `called` column is respected)")
    report = CellAssemblyReport(cells=len(grouped))
    order = sorted(grouped)
    work = []
    for cell in order:
        molecules = grouped[cell]
        report.molecules_in += len(molecules)
        report.molecules_by_cell[cell] = len(molecules)
        work.append((cell,
                     [trim(s) if trim_adapter else s for s, _ in molecules],
                     [w for _, w in molecules]))

    # Never: SERIAL, and a process pool was tried and reverted. The layout is pure Python, so a
    # pool is the obvious upgrade and cells are independent -- but on macOS the default start
    # method is `spawn`, which re-imports the caller's `__main__` in every worker and so crashes
    # any script, notebook or test module calling this without an `if __name__ == "__main__":`
    # guard; and `forkserver`, which does not touch `__main__`, dies here with "did not receive
    # acknowledgement of fd" before a single cell is assembled. 54 ms per cell, 479 cells in 26 s,
    # so the upgrade buys minutes on a 10k-cell library and costs a crash for every library
    # caller. ponytail: revisit when the assembler itself moves to C++, where threads need none
    # of this.
    one = partial(_assemble_one, k=k, min_overlap=min_overlap, min_identity=min_identity,
                  min_molecules=min_molecules, split=split,
                  min_minor_fraction=min_minor_fraction,
                  min_minor_molecules=min_minor_molecules,
                  min_informative_columns=min_informative_columns)
    all_contigs: list[Contig] = []
    results = (one(w) for w in work)
    for contigs in results:
        for contig in contigs:
            if len(contig.sequence) < min_length:
                report.contigs_dropped_small += 1
                continue
            all_contigs.append(contig)
            report.molecules_placed += contig.molecules
            report.contig_lengths.append(len(contig.sequence))
    report.contigs = len(all_contigs)

    output = Path(output)
    with open(output, "w") as fh:
        for contig in all_contigs:
            fh.write(f">{contig.contig_id} molecules={contig.molecules} reads={contig.reads}\n")
            fh.write(f"{contig.sequence}\n")
    return all_contigs, report


#: Which locus of a class carries the rearrangement that is unique per cell. Two of these in one
#: cell is the doublet signal; two light chains are not -- allelic inclusion at TRA and IGK/IGL
#: is a known, real population.
HEAVY = {"TRB", "TRD", "IGH"}
LIGHT = {"TRA", "TRG", "IGK", "IGL"}


def _is_productive(value) -> bool:
    return str(value).strip().upper() in ("T", "TRUE", "1", "YES")


def pair_chains(
    rows: list[dict],
    *,
    min_extra_molecules: int = 3,
    min_extra_fraction: float = 0.2,
    require_productive: bool = True,
) -> list[dict]:
    """Turn annotated contigs into one row per (cell, locus, junction), with a doublet flag.

    ``rows`` need ``cell_id``, ``locus``, ``junction_aa`` and ``molecules``; ``productive``,
    ``reads``, ``sequence`` and the gene calls are carried when present. Within a cell and locus
    the chains are ranked by supporting molecules. Rank 1 is ``primary``. A lower-ranked chain is
    believed only if it is productive and clears BOTH support gates, and then a second HEAVY
    chain is ``doublet_candidate`` while a second LIGHT chain is ``secondary`` -- allelic
    inclusion at TRA and IGK/IGL is a known, real population, not a doublet. Everything else is
    ``extra``, kept in the table with the columns that say why.

    Never: **an extra chain must be PRODUCTIVE, and that is the discriminator, not the count.**
    Measured on ``sc5p_v2_hs_PBMC_1k`` against Cell Ranger's own per-cell calls: of the extra
    chains it agrees with, 60/60 are productive; of those it does not, 20/128 are. Molecule
    support separates the same two groups far less sharply (median 10 against 4).

    Never: **the threshold is on MOLECULES, not reads, and 1 is not a threshold.** Sweeping the
    extra chain's molecule support against Cell Ranger's calls, before the productivity gate:

    ======  ====  =========  ======
    locus   min   precision  recall
    ======  ====  =========  ======
    TRA     1     0.360      1.000
    TRA     2     0.946      0.972
    TRA     3     1.000      0.917
    TRB     1     0.034      1.000
    TRB     3     0.727      1.000
    TRB     5     0.889      1.000
    ======  ====  =========  ======

    An extra chain supported by exactly one molecule is contamination in 96-97% of cases and a
    real second chain in none of the TRB cases; the distribution is bimodal, not a tail. That is
    why :func:`assemble_cells` defaults ``min_molecules`` to 2 and never to 1.
    """
    by_cell: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        cell, locus, junction = row.get("cell_id"), row.get("locus"), row.get("junction_aa")
        if not cell or not locus or not junction:
            continue
        by_cell[cell][locus].append(row)

    out: list[dict] = []
    for cell in sorted(by_cell):
        for locus, chains in sorted(by_cell[cell].items()):
            merged: dict[str, dict] = {}
            for row in chains:
                key = row["junction_aa"]
                seen = merged.get(key)
                if seen is None:
                    merged[key] = dict(row)
                    merged[key]["contigs"] = 1
                    merged[key]["contig_length"] = len(row.get("sequence") or "")
                else:
                    seen["molecules"] = seen.get("molecules", 0) + row.get("molecules", 0)
                    seen["reads"] = seen.get("reads", 0) + row.get("reads", 0)
                    seen["contigs"] += 1
                    seen["contig_length"] = max(seen["contig_length"],
                                                len(row.get("sequence") or ""))
                    # Productive anywhere is productive: a second contig of the same junction is
                    # the same rearrangement seen over a shorter window.
                    if _is_productive(row.get("productive")):
                        seen["productive"] = "T"
            ranked = sorted(merged.values(),
                            key=lambda r: (-r.get("molecules", 0), r["junction_aa"]))
            top = ranked[0].get("molecules", 0) if ranked else 0
            for rank, row in enumerate(ranked, start=1):
                molecules = row.get("molecules", 0)
                fraction = molecules / top if top else 0.0
                productive = _is_productive(row.get("productive"))
                if rank == 1:
                    status = "primary"
                elif require_productive and not productive:
                    status = "extra"
                elif molecules < min_extra_molecules or fraction < min_extra_fraction:
                    status = "extra"
                elif locus in HEAVY:
                    status = "doublet_candidate"
                else:
                    status = "secondary"
                out.append({
                    "cell_id": cell,
                    "locus": locus,
                    "chain_rank": rank,
                    "junction_aa": row["junction_aa"],
                    "junction": row.get("junction", ""),
                    "v_call": row.get("v_call", ""),
                    "j_call": row.get("j_call", ""),
                    "c_call": row.get("c_call", ""),
                    "productive": "T" if productive else "F",
                    "molecules": molecules,
                    "reads": row.get("reads", 0),
                    "contigs": row.get("contigs", 1),
                    "contig_length": row.get("contig_length", 0),
                    "fraction_of_top": round(fraction, 4),
                    "status": status,
                })
    return out


def cell_rank(molecules_per_cell: dict[str, int]) -> list[dict]:
    """The knee curve: cells ranked by how many molecules they carry.

    Never: the y axis is MOLECULES, never reads. One over-amplified molecule otherwise puts an
    empty droplet high on the curve, which is the artefact the plot exists to show.
    """
    ordered = sorted(molecules_per_cell.items(), key=lambda kv: (-kv[1], kv[0]))
    return [{"rank": i, "cell_id": cell, "molecules": n}
            for i, (cell, n) in enumerate(ordered, start=1)]


def find_knee(sizes: list[int], *, min_ratio_to_mean: float = 10.0) -> tuple[int, int]:
    """Where a descending molecules-per-barcode curve breaks. ``(rank, molecules)``, 1-based.

    ``sizes`` must be sorted descending. Returns ``(0, 0)`` when the curve has no knee worth
    reporting, which is a real answer and not a failure.

    This is Kneedle (`Satopaa 2011 <https://raghavan.usc.edu/papers/kneedle-simplex11.pdf>`_)
    taken at the **global maximum** of its difference curve. On a decreasing curve the distance to
    the chord joining the endpoints and Kneedle's normalised ``|y_n - (1 - x_n)|`` are the same
    quantity up to a constant; what differs is the selection rule.

    Never: **Kneedle's published rule -- walk the LOCAL maxima, stop at the first that then falls
    by a sensitivity step -- needs the smoothing spline the paper specifies, and without it the
    answer is degenerate.** Measured on a 136,032-barcode 10x library: unsmoothed, the difference
    curve carries 287 local maxima and the walk stops at rank 15 (1,057 molecules), which would
    call fifteen cells. Smoothed over a 0.01 window in normalised log-rank it gives rank 358 (314
    molecules) -- the global maximum, rank 376 (308 molecules), to within 5% -- and over-smoothing
    drifts it off again (rank 180 at a 0.1 window). The global maximum reaches the same place with
    no smoothing parameter to tune.

    Never: **the guard is not decoration.** A global maximum always exists, so an ambient-only
    library returns a rank too, and a "knee" there is a threshold that calls every barcode a cell.
    The floor is ``min_ratio_to_mean`` times the mean molecules per barcode -- one order of
    magnitude, which is the unit a log-log curve is read in. The mean alone is too weak and was
    tried: an ambient-only library of 1-3 molecules per barcode is a step curve whose corner sits
    at 3 against a mean of 2.0 and clears a mean-only test at 1.5x. On the real library the knee
    is 308 against a mean of 3.44, a factor of 89. Two orders apart, so the boundary is not tuned.
    """
    n = len(sizes)
    if n < 3 or sizes[0] <= 0 or sizes[-1] <= 0:
        return 0, 0

    x0, y0 = 0.0, math.log10(sizes[0])
    x1, y1 = math.log10(n), math.log10(sizes[-1])
    best, rank, molecules = -1.0, 0, 0
    for i in range(1, n - 1):
        x, y = math.log10(i + 1), math.log10(sizes[i])
        distance = abs((y1 - y0) * x - (x1 - x0) * y + x1 * y0 - y1 * x0)
        if distance > best:
            best, rank, molecules = distance, i + 1, sizes[i]
    mean = sum(sizes) / n
    if molecules < min_ratio_to_mean * mean:
        return 0, 0
    return rank, molecules


def cell_summary(
    pairs: list[dict],
    molecules_per_cell: dict[str, int] | None = None,
) -> list[dict]:
    """One row per cell: what chains it holds and what that makes it.

    ``pairs`` is :func:`pair_chains` output. The row carries the top two HEAVY and top two LIGHT
    chains by molecule support, which is what the doublet scatter is drawn from -- a plot of
    ``heavy_molecules_2`` against ``heavy_molecules_1`` separates the ambient cloud along the
    bottom from the multiplets off the axis.

    ``cell_status`` is one of:

    =====================  ==============================================================
    value                  meaning
    =====================  ==============================================================
    ``paired``             one heavy and at least one light chain
    ``heavy_only``         a heavy chain and no light one
    ``light_only``         a light chain and no heavy one
    ``doublet_candidate``  two heavy chains both clearing the support threshold
    ``no_chain``           no chain was called at all
    =====================  ==============================================================
    """
    molecules_per_cell = molecules_per_cell or {}
    by_cell: dict[str, list[dict]] = defaultdict(list)
    for row in pairs:
        by_cell[row["cell_id"]].append(row)

    out: list[dict] = []
    for cell in sorted(set(by_cell) | set(molecules_per_cell)):
        rows = by_cell.get(cell, [])
        heavy = sorted((r for r in rows if r["locus"] in HEAVY),
                       key=lambda r: -r["molecules"])
        light = sorted((r for r in rows if r["locus"] in LIGHT),
                       key=lambda r: -r["molecules"])
        doublet = any(r["status"] == "doublet_candidate" for r in rows)
        if doublet:
            status = "doublet_candidate"
        elif heavy and light:
            status = "paired"
        elif heavy:
            status = "heavy_only"
        elif light:
            status = "light_only"
        else:
            status = "no_chain"
        out.append({
            "cell_id": cell,
            "molecules": molecules_per_cell.get(cell, 0),
            "chains": len(rows),
            "loci": ",".join(sorted({r["locus"] for r in rows})),
            "heavy_locus": heavy[0]["locus"] if heavy else "",
            "heavy_junction_aa": heavy[0]["junction_aa"] if heavy else "",
            "heavy_molecules_1": heavy[0]["molecules"] if heavy else 0,
            "heavy_molecules_2": heavy[1]["molecules"] if len(heavy) > 1 else 0,
            "light_locus": light[0]["locus"] if light else "",
            "light_junction_aa": light[0]["junction_aa"] if light else "",
            "light_molecules_1": light[0]["molecules"] if light else 0,
            "light_molecules_2": light[1]["molecules"] if len(light) > 1 else 0,
            "cell_status": status,
        })
    return out


def threshold_sweep(
    pairs: list[dict],
    reference: set[tuple[str, str, str]],
    *,
    thresholds: tuple[int, ...] = (1, 2, 3, 5, 8, 12),
) -> list[dict]:
    """Precision and recall of the extra-chain filter, per locus, against a reference call set.

    ``reference`` is ``{(cell_id, locus, junction_aa)}`` from a tool you are willing to score
    against -- Cell Ranger's ``filtered_contig_annotations.csv``, a hashtag demultiplex, a
    simulation's truth. Only NON-primary chains are swept: the top chain of a cell is kept at
    every threshold, so including it would put a constant in both columns and flatten the curve.
    Each locus is swept twice, with the productivity gate off and on, so the two levers are
    separable rather than confounded.

    This is the table that picks ``min_extra_molecules`` and ``require_productive``; both are
    emitted rather than baked in, because the right setting is a property of the library's
    ambient load and not of the code.
    """
    extras: dict[tuple[str, bool], list[tuple[int, bool]]] = defaultdict(list)
    for row in pairs:
        if row["chain_rank"] == 1:
            continue
        key = (row["cell_id"], row["locus"], row["junction_aa"])
        hit = key in reference
        extras[(row["locus"], False)].append((row["molecules"], hit))
        if _is_productive(row.get("productive")):
            extras[(row["locus"], True)].append((row["molecules"], hit))

    out: list[dict] = []
    for locus, productive_only in sorted(extras, key=lambda kv: (kv[0], kv[1])):
        items = extras[(locus, productive_only)]
        # Recall is against every reference extra chain of this locus, so turning the gate on
        # can only lower it -- a denominator that shrank with the filter would hide the cost.
        positives = sum(1 for _, ok in extras[(locus, False)] if ok)
        for threshold in thresholds:
            kept = [ok for n, ok in items if n >= threshold]
            true_positive = sum(kept)
            out.append({
                "locus": locus,
                "require_productive": "T" if productive_only else "F",
                "min_molecules": threshold,
                "extra_chains_total": len(items),
                "kept": len(kept),
                "true_positive": true_positive,
                "false_positive": len(kept) - true_positive,
                "precision": true_positive / len(kept) if kept else float("nan"),
                "recall": true_positive / positives if positives else float("nan"),
            })
    return out


def _write_tsv(rows: list[dict], path: Path, columns: tuple[str, ...] | None = None) -> Path:
    columns = columns or (tuple(rows[0]) if rows else ())
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(_fmt(row.get(c, "")) for c in columns) + "\n")
    return path


def _fmt(value) -> str:
    if isinstance(value, float):
        return "NA" if value != value else f"{value:.4f}"
    return str(value)


def read_reference_calls(path: str | Path) -> set[tuple[str, str, str]]:
    """``{(cell_id, locus, junction_aa)}`` from a comparator's own call table.

    Accepts Cell Ranger's ``filtered_contig_annotations.csv`` (``barcode``/``chain``/``cdr3``,
    gem-group suffix stripped, ``is_cell`` respected when present) and any TSV/CSV carrying
    ``cell_id``/``locus``/``junction_aa``. This is only ever a *scoring* input --
    :func:`assemble_cells` never reads a reference of any kind.
    """
    import csv

    path = Path(path)
    with open_text(path) as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delimiter = "," if sample.count(",") > sample.count("\t") else "\t"
        out: set[tuple[str, str, str]] = set()
        for row in csv.DictReader(fh, delimiter=delimiter):
            if str(row.get("is_cell", "true")).strip().lower() not in ("true", "1", "yes"):
                continue
            cell = (row.get("cell_id") or row.get("barcode") or "").split("-", 1)[0]
            locus = (row.get("locus") or row.get("chain") or "").strip()
            junction = (row.get("junction_aa") or row.get("cdr3") or "").strip()
            if cell and locus and junction and junction != "None":
                out.add((cell, locus, junction))
    return out


def read_called_cells(path: str | Path) -> set[str]:
    """The called barcodes, from a bare one-per-line list or from migec's ``<sample>.cells.tsv``.

    A header naming a ``cell`` (or ``cell_id``/``barcode``) column switches to table mode, and a
    ``called`` column is then respected: migec writes every barcode it saw with a boolean, so
    taking column 1 unfiltered would restore the empty droplets this option exists to remove.
    A ``-1`` gem-group suffix is stripped.
    """
    out: set[str] = set()
    with open_text(path) as fh:
        first = fh.readline().rstrip("\n")
        fields = first.split("\t")
        header = [f.strip().lower() for f in fields]
        cell_col = next((i for i, f in enumerate(header)
                         if f in ("cell", "cell_id", "barcode")), None)
        if cell_col is None:
            if first.strip():
                out.add(first.strip().split("\t")[0].split("-", 1)[0])
            for line in fh:
                if line.strip():
                    out.add(line.strip().split("\t")[0].split("-", 1)[0])
            return out
        called_col = next((i for i, f in enumerate(header)
                           if f in ("called", "is_cell")), None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= cell_col:
                continue
            if called_col is not None and len(parts) > called_col:
                if parts[called_col].strip().lower() not in ("true", "1", "yes", "t"):
                    continue
            out.add(parts[cell_col].split("-", 1)[0])
    return out


def run(
    reads: str | Path,
    prefix: str | Path,
    *,
    organism: str = "human",
    cell_from: str = "auto",
    cell_regex: str | None = None,
    k: int = K,
    min_overlap: int = 30,
    min_identity: float = 0.95,
    min_molecules: int = 2,
    min_length: int = 100,
    trim_adapter: bool = True,
    split: bool = True,
    min_minor_fraction: float = 0.2,
    min_minor_molecules: int = 2,
    min_informative_columns: int = 2,
    min_extra_molecules: int = 3,
    min_extra_fraction: float = 0.2,
    require_productive: bool = True,
    cells_file: str | Path | None = None,
    reference: str | Path | None = None,
    threads: int = 0,
    map_d: bool = True,
    plot: str | None = None,
    gnuplot: str | None = None,
) -> dict:
    """Assemble, annotate, pair and diagnose a droplet V(D)J library. Returns the run report.

    Writes, all under ``prefix``:

    =======================  =====================================================================
    suffix                   what is in it
    =======================  =====================================================================
    ``.contigs.fasta``       the assembled per-cell contigs
    ``.contigs.airr.tsv``    their AIRR annotation, ``cell_id`` included
    ``.chains.tsv``          one row per (cell, locus, junction), ranked, with a status
    ``.cells.tsv``           one row per cell, with the doublet scatter's two columns
    ``.cell_rank.tsv``       the knee curve, molecules per cell
    ``.contig_lengths.tsv``  the contig length histogram
    ``.sweep.tsv``           precision/recall of the extra-chain filter (needs ``reference``)
    ``.partition.tsv``       clustering agreement of the clonotypes; see :mod:`arda.partition`
    ``.arda.json``           the run report
    =======================  =====================================================================
    """
    from .annotate.contig import reannotate_contigs

    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    def out(suffix: str) -> Path:
        # Never: `with_suffix` REPLACES the last dotted field, and a prefix may legitimately
        # carry one (`out/PBMC.rep2`). Append.
        return prefix.with_name(prefix.name + suffix)

    called = read_called_cells(cells_file) if cells_file else None
    contigs, report = assemble_cells(
        reads, out(".contigs.fasta"),
        cell_from=cell_from, cell_regex=cell_regex, k=k, min_overlap=min_overlap,
        min_identity=min_identity, min_molecules=min_molecules, min_length=min_length,
        trim_adapter=trim_adapter, split=split, min_minor_fraction=min_minor_fraction,
        min_minor_molecules=min_minor_molecules,
        min_informative_columns=min_informative_columns, cells=called,
    )
    molecules_per_cell = report.molecules_by_cell

    annotated = reannotate_contigs(
        [(c.contig_id, c.sequence) for c in contigs], organism,
        threads=threads, map_d=map_d,
    ) if contigs else []
    support = {c.contig_id: c for c in contigs}
    rows: list[dict] = []
    for record in annotated:
        contig = support.get(record.get("sequence_id", ""))
        if contig is None:
            continue
        rows.append({**record, "cell_id": contig.cell,
                     "molecules": contig.molecules, "reads": contig.reads})
    _write_tsv(rows, out(".contigs.airr.tsv"))

    chains = pair_chains(rows, min_extra_molecules=min_extra_molecules,
                         min_extra_fraction=min_extra_fraction,
                         require_productive=require_productive)
    cells = cell_summary(chains, molecules_per_cell)
    _write_tsv(chains, out(".chains.tsv"), (
        "cell_id", "locus", "chain_rank", "junction_aa", "junction",
        "v_call", "j_call", "c_call", "productive", "molecules", "reads", "contigs",
        "contig_length", "fraction_of_top", "status"))
    _write_tsv(cells, out(".cells.tsv"))
    ranked = cell_rank(molecules_per_cell)
    knee_rank, knee_molecules = find_knee([r["molecules"] for r in ranked])
    # Marked in the table rather than recomputed by the plot, so the figure and the number in the
    # report cannot drift apart.
    for row in ranked:
        row["knee"] = 1 if row["rank"] == knee_rank else 0
    _write_tsv(ranked, out(".cell_rank.tsv"))
    lengths: dict[int, int] = defaultdict(int)
    for length in report.contig_lengths:
        lengths[length // 25 * 25] += 1
    _write_tsv([{"length_bin": b, "contigs": lengths[b]} for b in sorted(lengths)],
               out(".contig_lengths.tsv"))

    summary = report.as_dict()
    summary["chains"] = len(chains)
    for status in ("paired", "heavy_only", "light_only", "doublet_candidate", "no_chain"):
        summary[f"cells_{status}"] = sum(1 for c in cells if c["cell_status"] == status)
    summary["min_extra_molecules"] = min_extra_molecules
    summary["min_extra_fraction"] = min_extra_fraction
    summary["require_productive"] = require_productive
    summary["split"] = split
    summary["knee_rank"] = knee_rank
    summary["knee_molecules"] = knee_molecules

    if reference is not None:
        from .partition import compare_partitions

        truth = read_reference_calls(reference)
        _write_tsv(threshold_sweep(chains, truth), out(".sweep.tsv"))
        called = {(r["cell_id"], r["locus"], r["junction_aa"])
                  for r in chains if r["status"] != "extra"}
        summary["reference_chains"] = len(truth)
        summary["reference_shared"] = len(called & truth)
        summary["reference_recall"] = len(called & truth) / len(truth) if truth else float("nan")

        # Cells clustered by the clonotype they carry, ours against the reference's. Per-chain
        # recall cannot see a split clone or a merged pair; this can.
        truth_label: dict[str, frozenset] = defaultdict(frozenset)
        called_label: dict[str, frozenset] = defaultdict(frozenset)
        for cell, locus, junction in truth:
            truth_label[cell] |= {(locus, junction)}
        for cell, locus, junction in called:
            called_label[cell] |= {(locus, junction)}
        shared_cells = sorted(set(truth_label) & set(called_label))
        if shared_cells:
            scores = compare_partitions(
                [truth_label[c] for c in shared_cells],
                [called_label[c] for c in shared_cells],
            )
            scores["cells_reference"] = len(truth_label)
            scores["cells_called"] = len(called_label)
            _write_tsv([{"metric": k, "value": v} for k, v in scores.items()],
                       out(".partition.tsv"))
            summary["partition"] = scores

    import json
    out(".arda.json").write_text(json.dumps(summary, indent=2) + "\n")
    if plot:
        from .scplot import draw
        summary["figures"] = [str(p) for p in draw(prefix, fmt=plot, gnuplot=gnuplot)]
    return summary
