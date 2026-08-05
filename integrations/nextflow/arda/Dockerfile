# Container for the Nextflow ARDA module: arda-mapper (PyPI) + the mmseqs2 binary.
# Build and push to your own registry (e.g. the ISPRAS private registry), then point the module's
# `container` directive (or a withName override in your nextflow.config) at the pushed tag:
#
#   docker build -t arda-mapper:2.7.0 integrations/nextflow/arda
#   docker tag arda-mapper:2.7.0 <your-registry>/arda-mapper:2.7.0
#   docker push <your-registry>/arda-mapper:2.7.0
#
FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
# fail the build early if arda or mmseqs did not install. `import seqtree` too: it comes from the
# `[rnaseq]` extra, and `arda --version` happily succeeds without it -- the failure would only
# surface mid-pipeline, in `arda rnaseq correct`.
RUN arda --version && mmseqs version && python -c "import seqtree"
