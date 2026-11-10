# Container with all software needed for the popgen tutorial.
#
# Build from the repository root:
#   docker build -t popgen_tutorial -f container/Dockerfile .
FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.source="https://github.com/axeljen/popgen_tutorial"
LABEL org.opencontainers.image.description="Software for the popgen tutorial (plink2, ADMIXTURE, ANGSD, PCAngsd, samtools, bedtools, IQ-TREE 2, newick-utils, ASTRAL-IV, R/tidyverse, fastaparser)"

# fastaparser has no releases, so we pin an exact commit
ARG FASTAPARSER_COMMIT=7e5b5e962e9bb62362556435a4f955ea2eafce5f

# install the conda environment into the base environment of the image
COPY --chown=$MAMBA_USER:$MAMBA_USER container/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# install fastaparser (its dependencies are already installed from environment.yml)
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-deps --no-cache-dir \
        "https://github.com/axeljen/fastaparser/archive/${FASTAPARSER_COMMIT}.tar.gz"

# Apptainer ignores the micromamba entrypoint that normally activates the
# environment, so put it on PATH explicitly
ENV PATH=/opt/conda/bin:$PATH
ENV LANG=C.UTF-8

# small script that checks that all software works:
#   apptainer exec popgen_tutorial.sif check_popgen_software
COPY container/check_popgen_software.sh /usr/local/bin/check_popgen_software
