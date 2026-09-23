# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "polars", "altair", "pyarrow", "pandas"]
# ///
"""Single-cell QC: the knee, the doublet scatter, the filter sweep, the clustering scores.

Run with:  uvx marimo edit notebooks/singlecell_qc.py

Every figure is drawn from a TSV that `arda cells` wrote, never recomputed here, so any of them
can be redrawn from a committed table long after the run. Point `PREFIX` at your own output and
everything below follows; with no output to hand it falls back to a small built-in example so the
notebook still runs and still says what each panel means.
"""

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # Reading a droplet V(D)J run

        `arda cells` assembles each cell's contigs with no germline reference, annotates them,
        ranks the chains within each cell and writes six tables. Four questions:

        1. **Which barcodes are cells?** The knee curve (`.cell_rank.tsv`) — molecules per
           barcode against rank. Never reads: one over-amplified molecule otherwise puts an empty
           droplet high on the curve, which is the artifact the plot exists to show.
        2. **Which cells are doublets?** The second heavy chain's support against the first
           (`.cells.tsv`). Two rearranged heavy chains in one droplet is a doublet; two light
           chains is allelic inclusion and is real.
        3. **Where should the filter sit?** Precision and recall of the extra-chain rule
           (`.sweep.tsv`), swept with the productivity gate off and on.
        4. **Did we recover the same clones, or just the same chains?** The homogeneity-parsimony
           scores (`.partition.tsv`). Recall cannot see a split clone or a merged pair.
        """
    )
    return (mo,)


@app.cell
def _():
    import os
    from pathlib import Path

    import polars as pl

    # Point this at your own run: everything below is `PREFIX + <suffix>`.
    PREFIX = Path(os.environ.get("ARDA_CELLS_PREFIX", "out/PBMC"))

    def load(suffix: str) -> pl.DataFrame | None:
        path = PREFIX.with_name(PREFIX.name + suffix)
        return pl.read_csv(path, separator="\t") if path.exists() else None

    return PREFIX, Path, load, pl


@app.cell
def _(load, mo, pl):
    # A tiny stand-in so the notebook runs with nothing on disk. It is deliberately obvious:
    # ten cells, one of them a doublet. Never mistake it for a result.
    def _example_cells() -> pl.DataFrame:
        return pl.DataFrame({
            "cell_id": [f"CELL{i:02d}" for i in range(10)],
            "molecules": [400, 310, 260, 240, 190, 150, 120, 90, 60, 40],
            "heavy_molecules_1": [90, 70, 64, 58, 44, 35, 28, 20, 14, 9],
            "heavy_molecules_2": [0, 52, 0, 1, 0, 0, 2, 0, 0, 0],
            "cell_status": ["paired", "doublet_candidate", "paired", "paired", "paired",
                            "paired", "paired", "heavy_only", "paired", "light_only"],
        })

    cells = load(".cells.tsv")
    using_example = cells is None
    if using_example:
        cells = _example_cells()
    mo.md(
        "**Using the built-in example** — set `ARDA_CELLS_PREFIX` to your run's prefix."
        if using_example else
        f"Loaded **{cells.height}** cells."
    )
    return cells, using_example


@app.cell
def _(cells, mo):
    counts = cells["cell_status"].value_counts().sort("count", descending=True)
    mo.ui.table(counts, label="cells by status")
    return (counts,)


@app.cell
def _(load, mo, pl):
    import altair as alt

    rank = load(".cell_rank.tsv")
    if rank is None:
        rank = pl.DataFrame({"rank": list(range(1, 11)),
                             "molecules": [400, 310, 260, 240, 190, 150, 120, 90, 60, 40]})
    knee = alt.Chart(rank.to_pandas()).mark_line().encode(
        x=alt.X("rank:Q", scale=alt.Scale(type="log"), title="cell barcode rank"),
        y=alt.Y("molecules:Q", scale=alt.Scale(type="log"), title="molecules (unique UMIs)"),
    ).properties(width=560, height=320, title="Barcode rank plot")
    mo.vstack([mo.md("## 1. The knee"), mo.ui.altair_chart(knee)])
    return alt, knee, rank


@app.cell
def _(alt, cells, mo):
    # Only cells that HAVE a second heavy chain: a cell with one has y = 0, which a log axis
    # cannot show, and substituting a floor value would invent the datum the plot is read for.
    two = cells.filter(cells["heavy_molecules_2"] > 0)
    scatter = alt.Chart(two.to_pandas()).mark_point(filled=True, size=70).encode(
        x=alt.X("heavy_molecules_1:Q", scale=alt.Scale(type="log"),
                title="molecules on the top heavy chain"),
        y=alt.Y("heavy_molecules_2:Q", scale=alt.Scale(type="log"),
                title="molecules on the second heavy chain"),
        color=alt.Color("cell_status:N", title="status"),
        tooltip=["cell_id", "heavy_molecules_1", "heavy_molecules_2", "cell_status"],
    ).properties(width=520, height=380, title="Doublet scatter")
    mo.vstack([
        mo.md(
            "## 2. Doublets\n\n"
            f"**{two.height}** of {cells.height} cells carry a second heavy chain at all. The "
            "ambient cloud sits along the bottom — a second chain on one molecule is "
            "contamination 96-97% of the time. A doublet sits off the axis, at a median 0.63 of "
            "the top chain's support."
        ),
        mo.ui.altair_chart(scatter),
    ])
    return scatter, two


@app.cell
def _(alt, load, mo):
    sweep = load(".sweep.tsv")
    if sweep is None:
        _out = mo.md(
            "## 3. The filter sweep\n\n"
            "_Not written._ Re-run with `--reference <calls.csv>` — Cell Ranger's "
            "`filtered_contig_annotations.csv`, a hashtag demultiplex, or a simulation's truth — "
            "and `.sweep.tsv` reports precision and recall of the extra-chain rule per locus, "
            "with the productivity gate off and on. That is the table the defaults came from."
        )
    else:
        long = sweep.unpivot(
            index=["locus", "require_productive", "min_molecules"],
            on=["precision", "recall"], variable_name="metric", value_name="value")
        chart = alt.Chart(long.to_pandas()).mark_line(point=True).encode(
            x=alt.X("min_molecules:Q", title="minimum molecules on the extra chain"),
            y=alt.Y("value:Q", scale=alt.Scale(domain=[0, 1])),
            color="metric:N",
            strokeDash="require_productive:N",
            column=alt.Column("locus:N"),
        ).properties(width=240, height=260)
        _out = mo.vstack([
            mo.md("## 3. The filter sweep\n\nSolid against dashed is the productivity gate. It "
                  "roughly doubles precision at every threshold and is the discriminator; the "
                  "molecule count is the second gate, not the first."),
            mo.ui.altair_chart(chart),
        ])
    _out
    return (sweep,)


@app.cell
def _(load, mo):
    partition = load(".partition.tsv")
    if partition is None:
        _out = mo.md(
            "## 4. Clustering agreement\n\n_Not written._ Needs `--reference` as well."
        )
    else:
        scores = dict(zip(partition["metric"].to_list(), partition["value"].to_list()))
        singleton = float(scores.get("classes_singleton", 0))
        total = float(scores.get("items", 1))
        _out = mo.vstack([
            mo.md(
                "## 4. Clustering agreement\n\n"
                "Two scores per family: one punishes **merging** two clones, one punishes "
                "**splitting** one. Either alone is trivially maximised, so never quote one.\n\n"
                f"> **{singleton:.0f} of {total:.0f}** reference clonotypes are singletons. "
                "An unexpanded repertoire has almost no clustering to agree about and every "
                "score below is then decided by a handful of cells — read this line first."
            ),
            mo.ui.table(partition, label="partition scores"),
        ])
    _out
    return (partition,)


@app.cell
def _(load, mo):
    chains = load(".chains.tsv")
    if chains is None:
        _out = mo.md("_No `.chains.tsv`; nothing to inspect._")
    else:
        _out = mo.vstack([
            mo.md("## The chain table\n\nOne row per (cell, locus, junction), ranked within the "
                  "cell. `status` is `primary`, `secondary` (allelic inclusion), "
                  "`doublet_candidate` (a second heavy chain) or `extra` (filtered out, kept "
                  "with the columns that say why)."),
            mo.ui.table(chains.head(50), label="chains (first 50 rows)"),
        ])
    _out
    return (chains,)


if __name__ == "__main__":
    app.run()
