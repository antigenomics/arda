"""The Nextflow module, checked against the CLI it actually invokes.

`integrations/nextflow/arda/` ships in a directory no test touched, and that is not a theoretical
gap: the module sat pinned to arda 2.7.2 for four releases, and its Dockerfile's acceptance check
grepped `arda rnaseq --help` for four flags that mode has never exposed -- so `docker build`
failed on a perfectly good install, silently, because no CI job builds it.

Nothing here runs Nextflow or Docker. It asserts the two things that can be wrong without anyone
noticing: the version pins drifting away from the package, and the module naming a command or a
flag the CLI does not have. Same argument as `test_snakemake_integration.py`.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from arda import __version__
from arda.cli import app

MOD = Path(__file__).resolve().parents[2] / "integrations/nextflow/arda"
MAIN_NF = (MOD / "main.nf").read_text()
CONFIG = (MOD / "nextflow.config").read_text()
DOCKERFILE = (MOD / "Dockerfile").read_text()
ENVIRONMENT = (MOD / "environment.yml").read_text()

runner = CliRunner()


def _help(*argv: str) -> str:
    """`--help` text at a width typer will not truncate long flag names at."""
    res = runner.invoke(app, [*argv, "--help"], env={"COLUMNS": "200"})
    assert res.exit_code == 0, f"arda {' '.join(argv)} --help exited {res.exit_code}"
    return res.stdout


@pytest.mark.parametrize("name", ["main.nf", "nextflow.config", "environment.yml", "Dockerfile",
                                  "meta.yml", "README.md"])
def test_the_module_ships_every_file_it_documents(name):
    assert (MOD / name).is_file()


def test_every_version_pin_tracks_the_package(): 
    """Never: the pin lives in four places and NO job builds this module.

    It sat on 2.7.2 for four releases exactly because nothing compared it to `arda.__version__`.
    The release checklist says to bump it; this is what makes forgetting loud.
    """
    assert f"arda-mapper=={__version__}" in ENVIRONMENT
    assert f'container "arda-mapper:{__version__}"' in MAIN_NF
    assert f"arda-mapper:{__version__}" in DOCKERFILE
    readme = (MOD / "README.md").read_text()
    assert f"arda {__version__}" in readme or f"arda-mapper:{__version__}" in readme


def _protocol_map() -> dict[str, str]:
    """airrflow's `library_generation_method` -> the arda mode the module runs for it."""
    block = MAIN_NF[MAIN_NF.index("def protocols = ["):]
    block = block[:block.index("\n    ]")]
    return dict(re.findall(r"'([a-z0-9_]+)'\s*:\s*'([a-z]+)'", block))


def test_the_protocol_map_covers_airrflows_bulk_enum_and_nothing_invented():
    """Never: the protocol comes from airrflow's own field, so the keys must be ITS enum.

    A key airrflow cannot emit is dead config; a bulk method missing from the map is a sample the
    module refuses for no reason. `sc_10x_genomics` is deliberately absent -- it is refused, with
    a message, because `arda cells` consumes a per-molecule UMI consensus and not raw 10x reads.
    """
    airrflow_enum = {"specific_pcr_umi", "specific_pcr", "dt_5p_race", "dt_5p_race_umi",
                     "sc_10x_genomics", "trust4"}
    got = _protocol_map()
    assert set(got) == airrflow_enum - {"sc_10x_genomics"}
    assert "sc_10x_genomics" in MAIN_NF, "the single-cell refusal must name the method"


def test_every_mode_the_protocol_map_names_is_a_real_arda_command():
    for method, mode in sorted(_protocol_map().items()):
        assert _help(mode), f"{method} -> arda {mode}"


def test_single_cell_is_refused_rather_than_mapped_to_a_mode():
    """⛔ `arda cells` takes a positional UMI-consensus FASTQ and no --r1/--out-dir at all.

    Routing raw 10x reads at it would fail as a bad repertoire rather than as a wiring error.
    """
    assert "cells" not in _protocol_map().values()
    assert "single_cell" in MAIN_NF and "IllegalArgumentException" in MAIN_NF
    cells_help = _help("cells")
    for flag in ("--r1", "--out-dir", "--out-prefix"):
        assert flag not in cells_help, f"arda cells now has {flag}; revisit the refusal"


def test_every_flag_the_script_block_emits_exists_on_every_mode_it_targets():
    """⛔ The check that caught `--productive-only`, a flag arda has never had.

    The module writes ONE command template and runs it for several modes, so a flag is only safe
    if EVERY targeted mode offers it. Nothing else in this repo compares the two.
    """
    script = MAIN_NF[MAIN_NF.index('    """\n    ${mmseqs_pin}'):]
    script = script[:script.index("END_VERSIONS")]
    flags = sorted(set(re.findall(r"(--[a-z][a-z-]+)", script)))
    assert flags, "no flags found in the script block; update this test"
    for mode in sorted(set(_protocol_map().values())):
        text = _help(mode)
        missing = [f for f in flags if f not in text]
        assert not missing, f"arda {mode} does not offer {missing}"


def test_the_dockerfile_acceptance_check_only_names_flags_the_cli_has():
    """⛔ The live defect this test exists for.

    The check used to grep `arda rnaseq --help` for `--two-pass`, `--fast-segments`,
    `--v-only-on-segment` and `--indel-rescue`. Since 2.16.0 the regime IS the command name and
    `arda rnaseq` exposes none of them -- they live on `arda map` -- so the build failed on every
    correct install. Whatever command the Dockerfile greps, the flags must be ON that command.
    """
    m = re.search(r"COLUMNS=200 arda ([a-z-]+(?: [a-z-]+)*) --help > /tmp/h\.txt", DOCKERFILE)
    assert m, "the Dockerfile no longer greps a --help; update this test with it"
    argv = m.group(1).split()
    text = _help(*argv)
    flags = re.findall(r"(--[a-z][a-z-]+)", DOCKERFILE[m.end():DOCKERFILE.index("rm /tmp/h.txt")])
    assert flags, "no flags are asserted by the Dockerfile"
    missing = [f for f in flags if f not in text]
    assert not missing, f"arda {' '.join(argv)} --help does not offer {missing}"


def test_every_mode_the_dockerfile_asserts_exists():
    m = re.search(r"for m in ([a-z ]+); do", DOCKERFILE)
    assert m, "the Dockerfile no longer asserts the mode names"
    for mode in m.group(1).split():
        assert _help(mode), mode


def test_the_module_names_no_registry_and_no_institution():
    """Never: a vendor-neutral module. A deployment target belongs in the user's own config.

    `<your-registry>` is a placeholder and stays; a real host or organisation name does not.
    """
    for name, text in (("main.nf", MAIN_NF), ("nextflow.config", CONFIG),
                       ("Dockerfile", DOCKERFILE), ("environment.yml", ENVIRONMENT),
                       ("README.md", (MOD / "README.md").read_text())):
        low = text.lower()
        for token in ("ispras", "gamaleya", "гамале"):
            assert token not in low, f"{name} names {token!r}"


def test_no_bare_params_read_survives_in_the_module():
    """Never: Nextflow scans the script STATICALLY, so a bare `params.x` warns even when guarded.

    Under strict mode "Access to undefined parameter" is a hard failure, and this module must stay
    correct when it is included WITHOUT its nextflow.config. Every read goes through
    `params.getOrDefault(...)`; the module's own declared params are the only exceptions.
    """
    declared = set(re.findall(r"^\s{4}([a-z_]+)\s*=", CONFIG, re.M))
    # Comments are prose, not reads -- and this module's comments quote `params.x` while
    # explaining why not to write one. Scan the code.
    code = re.sub(r"//.*$", "", MAIN_NF, flags=re.M)
    bare = {m for m in re.findall(r"\bparams\.([a-zA-Z_]\w*)", code)
            if m not in ("getOrDefault",) and m not in declared}
    assert not bare, f"bare params reads in main.nf: {sorted(bare)}"


def test_conda_and_container_never_sit_at_bare_process_scope():
    """Never: at `process { }` scope they apply to EVERY process in the host pipeline.

    Nextflow then silently builds a different arda *and* a different aligner for tools that never
    asked for one. Every setting belongs inside `withName: 'ARDA' { ... }`.
    """
    head = CONFIG[CONFIG.index("process {"):]
    before_withname = head[:head.index("withName:")]
    for directive in ("conda", "container", "ext.args", "ext.when", "publishDir"):
        assert not re.search(rf"^\s+{re.escape(directive)}\s*=", before_withname, re.M), directive


def test_the_process_is_named_for_what_it_replaces_in_airrflow():
    assert "process ARDA_ASSIGN {" in MAIN_NF
    assert "withName: 'ARDA_ASSIGN'" in CONFIG


def test_the_module_does_not_wire_params_airrflow_applies_itself():
    """Never: a parameter accepted and silently doing nothing is this project's recurring failure.

    arda has no productive filter and stages no per-run germline DB, so `productive_only` and the
    `reference_*` params must be NAMED as unconsumed rather than read into a flag.
    """
    for name in ("productive_only", "reference_igblast"):
        assert name in CONFIG, f"{name} must be named in the config as deliberately unconsumed"
    code = re.sub(r"//.*$", "", MAIN_NF, flags=re.M)
    assert "productive_only" not in code, "productive_only must not be read in the script"
