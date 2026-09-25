# Container for the Nextflow ARDA module: arda-mapper (PyPI) + the mmseqs2 binary.
#
# Vendor-neutral by construction: it pins nothing but arda and its aligner, and names no registry.
# Build it, push it wherever you deploy from, then point the module's `container` directive (or a
# `withName` override in your nextflow.config) at the pushed tag:
#
#   docker build -t arda-mapper:2.29.0 integrations/nextflow/arda
#   docker tag  arda-mapper:2.29.0 <your-registry>/arda-mapper:2.29.0
#   docker push <your-registry>/arda-mapper:2.29.0
#
FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Fail the build early if the install is not usable. `import seqtree` because `arda --version`
# succeeds without it and the failure would otherwise surface mid-pipeline, in `arda rnaseq
# correct`. And assert on the CLI, not on the import: the console script is a separate resolution
# from the import, and a deploy into the wrong environment passes every import-level check while
# the CLI still says "No such option".
#
# ⛔ Assert the MODE NAMES, never a flag string. Since 2.16.0 the regime IS the command name and
# arda owns which flags that implies -- `arda rnaseq` does not expose `--two-pass` /
# `--fast-segments` / `--v-only-on-segment` / `--indel-rescue` at all, and an earlier version of
# this check grepped for exactly those against `arda rnaseq --help` and failed every build. The
# tuning knobs live on `arda map`, which is the A/B surface; that is where they are checked.
RUN arda --version && mmseqs version && python -c "import seqtree" \
 && for m in rnaseq amplicon cells map; do \
        arda "$m" --help > /dev/null || { echo "MISSING MODE: arda $m"; exit 1; }; \
    done \
 && COLUMNS=200 arda map --help > /tmp/h.txt \
 && for f in --two-pass --fast-segments --v-only-on-segment --prefilter --indel-rescue; do \
        grep -q -- "$f" /tmp/h.txt || { echo "MISSING FLAG: arda map $f"; exit 1; }; \
    done \
 && rm /tmp/h.txt
