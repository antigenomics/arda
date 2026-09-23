# Notebooks

[marimo](https://marimo.io) notebooks. Each is a plain Python file — readable in a diff, runnable
without a notebook server, and with no hidden execution order.

## Running them

Every notebook declares its own dependencies in a [PEP 723](https://peps.python.org/pep-0723/)
header, so `uv` builds the environment and nothing has to be installed first:

```bash
uvx marimo edit notebooks/singlecell_qc.py       # interactive, in a browser
uvx marimo run  notebooks/singlecell_qc.py       # read-only app view
uv run --with marimo --with polars --with altair --with pyarrow --with pandas \
       python notebooks/singlecell_qc.py         # just execute it, no UI
```

## What is here

| notebook | reads | answers |
|---|---|---|
| `singlecell_qc.py` | the tables `arda cells` writes | which barcodes are cells, which are doublets, where the extra-chain filter should sit, and whether the clones agree with a reference |

Point a notebook at your own run with an environment variable rather than editing it:

```bash
ARDA_CELLS_PREFIX=out/PBMC uvx marimo edit notebooks/singlecell_qc.py
```

With no run to hand, `singlecell_qc.py` falls back to a ten-cell built-in example so the notebook
still executes and still says what each panel means. It says so on the page — never mistake it for
a result.

Nothing here computes a pipeline result. Every figure is drawn from a TSV a run already wrote, so
any of them can be redrawn from a committed table long after the FASTQ is gone.
