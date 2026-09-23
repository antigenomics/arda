"""Single-cell QC figures, drawn with gnuplot from the tables ``arda cells`` already wrote.

Nothing here computes anything. Every panel is a gnuplot script over a TSV that is written
whether or not gnuplot exists, which is what lets a figure be redrawn months later by someone who
no longer has the FASTQ. It also means arda gains no plotting dependency: the ``.gp`` scripts are
written either way and gnuplot runs only if it is on PATH.

Colours are ColorBrewer Dark2 (qualitative, colour-blind safe), never hand-picked. Transparent
background and one mid-grey ink, so a single SVG serves a light page, a dark page and print.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

__all__ = ["draw", "PANELS"]

logger = logging.getLogger(__name__)

DARK2 = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]
INK = "#808080"

TERMINALS = {
    "svg": 'set terminal svg size 760,520 font "Helvetica,13" enhanced',
    "png": 'set terminal pngcairo size 1000,680 font "Helvetica,13" transparent',
    "pdf": 'set terminal pdfcairo size 6.4in,4.4in font "Helvetica,11"',
}

_PREAMBLE = f"""
set output "{{out}}"
set border 3 lw 1 lc rgb "{INK}"
set tics nomirror out textcolor rgb "{INK}"
set title textcolor rgb "{INK}"
set xlabel textcolor rgb "{INK}"
set ylabel textcolor rgb "{INK}"
set key inside bottom right box lw 0.5 lc rgb "{INK}" textcolor rgb "{INK}" samplen 2 spacing 1.1
set style fill solid 0.85 noborder
set grid ytics lw 0.5 lc rgb "{INK}" dt 3
set datafile separator "\\t"
set datafile missing "NA"
# Never: every panel addresses its columns by NAME, never by index. These tables gain columns --
# `productive` and `fraction_of_top` landed in the middle of `.chains.tsv` after the first draft
# -- and a positional `$9` then draws the wrong series silently, with no error and a plausible
# figure. `columnheaders` plus `column("name")` fails loudly instead.
set datafile columnheaders
"""

#: ``suffix -> (stem, title, body)``. The body is gnuplot with ``{data}`` for the TSV path.
PANELS: dict[str, tuple[str, str, str]] = {
    ".cell_rank.tsv": ("cell_rank", "Molecules per cell, ranked", f"""
set title "Cell barcode rank plot"
set xlabel "cell barcode rank"
set ylabel "molecules (unique UMIs)"
set logscale xy
set key inside bottom left box lw 0.5 lc rgb "{INK}" textcolor rgb "{INK}"
# The knee, drawn from the same table rather than recomputed: the rank whose molecule count is
# the largest in the file at or above the reported knee. `arda cells` writes the knee into its
# JSON report and marks the row here, so the line and the number cannot drift apart.
plot "{{data}}" using (column("rank")):(column("molecules")) \\
        with lines lw 2 lc rgb "{DARK2[0]}" title "barcodes", \\
     "{{data}}" using (strcol("knee") eq "1" ? column("rank") : 1/0):(column("molecules")) \\
        with points pt 7 ps 1.4 lc rgb "{DARK2[1]}" title "knee"
"""),
    ".cells.tsv": ("doublet_scatter", "Second heavy chain against first", f"""
set title "Doublet scatter: cells carrying two heavy chains"
set xlabel "molecules on the cell's top heavy chain"
set ylabel "molecules on its second heavy chain"
set logscale xy
set key inside top left box lw 0.5 lc rgb "{INK}" textcolor rgb "{INK}"
# Never: only cells that HAVE a second heavy chain are drawn. A cell with one chain has y = 0,
# which a log axis cannot show, and substituting a floor value would invent the datum the plot
# is being read for. The one-chain cells are the knee plot's population, not this one.
h1 = "heavy_molecules_1"; h2 = "heavy_molecules_2"
plot "{{data}}" using (column(h2) > 0 ? column(h1) : 1/0):(column(h2)) \\
        with points pt 7 ps 0.6 lc rgb "{DARK2[7]}" title "two heavy chains", \\
     "{{data}}" using (strcol("cell_status") eq "doublet_candidate" ? column(h1) : 1/0):(column(h2)) \\
        with points pt 7 ps 1.0 lc rgb "{DARK2[1]}" title "doublet candidate"
"""),
    ".contig_lengths.tsv": ("contig_lengths", "Assembled contig lengths", f"""
set title "Assembled contig length"
set xlabel "contig length (nt, 25 nt bins)"
set ylabel "contigs"
set logscale y
unset key
set boxwidth 22
plot "{{data}}" using (column("length_bin")):(column("contigs")) \\
        with boxes lc rgb "{DARK2[2]}"
"""),
    ".sweep.tsv": ("filter_sweep", "Extra-chain filter, precision against recall", f"""
set title "Extra-chain filter: precision and recall against the reference calls"
set xlabel "minimum molecules on the extra chain"
set ylabel "rate"
set yrange [0:1.05]
set key inside bottom left box lw 0.5 lc rgb "{INK}" textcolor rgb "{INK}"
prod(v) = strcol("require_productive") eq v ? column("min_molecules") : 1/0
plot "{{data}}" using (prod("F")):(column("precision")) \\
        with points pt 6 ps 0.9 lc rgb "{DARK2[0]}" title "precision, any chain", \\
     "{{data}}" using (prod("T")):(column("precision")) \\
        with points pt 7 ps 0.9 lc rgb "{DARK2[0]}" title "precision, productive only", \\
     "{{data}}" using (prod("F")):(column("recall")) \\
        with points pt 4 ps 0.9 lc rgb "{DARK2[1]}" title "recall, any chain", \\
     "{{data}}" using (prod("T")):(column("recall")) \\
        with points pt 5 ps 0.9 lc rgb "{DARK2[1]}" title "recall, productive only"
"""),
    ".chains.tsv": ("chain_support", "Molecules per chain, by rank", f"""
set title "Molecules supporting a cell's chains, by rank within the cell"
set xlabel "molecules on the chain (powers of two)"
set ylabel "chains"
set logscale xy
set key inside top right box lw 0.5 lc rgb "{INK}" textcolor rgb "{INK}"
bin(x) = x < 1 ? 1 : 2**(floor(log(x)/log(2)))
plot "{{data}}" using (column("chain_rank") == 1 ? bin(column("molecules")) : 1/0):(1) \\
        smooth frequency with steps lw 2 lc rgb "{DARK2[0]}" title "rank 1", \\
     "{{data}}" using (column("chain_rank") > 1 ? bin(column("molecules")) : 1/0):(1) \\
        smooth frequency with steps lw 2 lc rgb "{DARK2[3]}" title "rank 2+"
"""),
}


def draw(prefix: str | Path, *, fmt: str = "svg", gnuplot: str | None = None) -> list[Path]:
    """Write a ``.gp`` script per available table and render it if gnuplot is present.

    ``gnuplot=""`` writes the scripts and renders nothing, which is also how the tests exercise
    the no-gnuplot path. Returns the figures that were actually rendered.
    """
    if fmt not in TERMINALS:
        raise ValueError(f"unknown format {fmt!r}; expected one of {sorted(TERMINALS)}")
    prefix = Path(prefix)
    binary = shutil.which("gnuplot") if gnuplot is None else (gnuplot or None)
    figures: list[Path] = []
    for suffix, (stem, _title, body) in PANELS.items():
        data = prefix.with_name(prefix.name + suffix)
        if not data.exists():
            continue
        script = prefix.with_name(f"{prefix.name}.{stem}.gp")
        figure = prefix.with_name(f"{prefix.name}.{stem}.{fmt}")
        script.write_text(
            TERMINALS[fmt] + "\n"
            + _PREAMBLE.format(out=figure.name)
            + body.format(data=data.name)
        )
        if binary:
            # Never: a figure never kills the run. The tables above are the deliverable and the
            # panels are derived from them, so a gnuplot that refuses one panel must not throw
            # away the other four and the report -- which is what `check=True` did the first time
            # `--plot` met a library with no doublets in it: `doublet_scatter` filters every cell
            # without a second heavy chain to `1/0`, an empty log-scale plot is "x range is
            # invalid", and a clean run therefore crashed where a dirty one drew. The `.gp` script
            # is on disk either way, which is the same contract as having no gnuplot at all.
            done = subprocess.run([binary, script.name], cwd=str(script.parent),
                                  capture_output=True, text=True)
            if done.returncode:
                # gnuplot opens its output before it evaluates the plot, so a refused panel
                # leaves a ZERO-BYTE file behind. Unlink it: an empty .svg in the output
                # directory reads as a figure that exists, which is worse than one that does not.
                figure.unlink(missing_ok=True)
                tail = (done.stderr or "").strip().splitlines()
                logger.warning("scplot: %s not drawn -- gnuplot said: %s", figure.name,
                               tail[-1] if tail else "(no message)")
                continue
            figures.append(figure)
    return figures
