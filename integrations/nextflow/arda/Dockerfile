# Container for the Nextflow ARDA module: arda-mapper (PyPI) + the mmseqs2 binary.
# Build and push to your own registry (e.g. the ISPRAS private registry), then point the module's
# `container` directive (or a withName override in your nextflow.config) at the pushed tag:
#
#   docker build -t arda-mapper:2.12.0 integrations/nextflow/arda
#   docker tag arda-mapper:2.12.0 <your-registry>/arda-mapper:2.12.0
#   docker push <your-registry>/arda-mapper:2.12.0
#
FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Fail the build early if the install is not usable. `import seqtree` because `arda --version`
# succeeds without it and the failure would otherwise surface mid-pipeline, in `arda rnaseq
# correct`. And assert the tuning flags on the CLI, not on the module: the console script is a
# separate resolution from the import, and a deploy into the wrong environment passes every
# import-level check while the CLI still says "No such option".
# COLUMNS=200 is required -- typer truncates long flag names at 80 columns and the grep then fails
# on a perfectly good build.
RUN arda --version && mmseqs version && python -c "import seqtree" \
 && COLUMNS=200 arda rnaseq run --help > /tmp/h.txt \
 && for f in --two-pass --fast-segments --v-only-on-segment --prefilter --indel-rescue; do \
        grep -q -- "$f" /tmp/h.txt || { echo "MISSING FLAG: $f"; exit 1; }; \
    done \
 && rm /tmp/h.txt
