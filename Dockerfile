# Dockerfile
#
# Builds the Multi-Omics Bioinformatics Portal image entirely from
# environment.yml (this repo's single source of truth for Python, R/
# Bioconductor, and CLI tool dependencies) via micromamba -- so this image
# can never drift out of sync with what a local or HPC conda/mamba install
# would produce from the same environment.yml.
#
# --- Replaces the OLD, hand-maintained Dockerfile (2026-08-17) ---
# The old Dockerfile (rocker/r-ver base, manual apt/curl/conda/R/pip install
# steps scattered across the file) is being retired in favor of this one.
# Reviewing it surfaced several things this file's environment.yml now
# restores that were previously missing here: nextflow + its Java
# dependency, gffread, bioconductor-limma (REQUIRED for DESeq2's batch-
# adjusted PCA -- not optional), and r-duckdb/r-dbi/r-tidyr for the Spatial
# Transcriptomics pipeline (Nextflow -> R spatial core -> DuckDB, see
# spatial_workspace.py / main.nf / bin/process_spatial.R). See
# environment.yml's own header comment for the full list and reasoning.
#
# --- kaleido + Chrome (2026-08-17) ---
# python-kaleido v1+ (see environment.yml -- the old 0.2.* pin was a
# deprecated version) never bundles Chrome, regardless of whether it's
# installed via conda or pip -- a separate one-time `kaleido.get_chrome_sync()`
# call is required to fetch a compatible Chrome/Chromium build. The old
# Dockerfile had already correctly identified this and ran the equivalent
# pip-based prefetch at BUILD time (network access available then) rather
# than leaving it to happen at first real use at RUNTIME (where the
# container may have no outbound network access, silently breaking PDF
# export on its first real use). That same build-time prefetch is
# reproduced below, now against the conda-installed python-kaleido.
#
# --- Repo layout (CONFIRMED, 2026-08-17) ---
# repo root contains environment.yml/Dockerfile/bin/(process_spatial.R)/
# main.nf/nextflow.config, with the Streamlit app itself in an app/
# subfolder (repo/app/app.py, repo/app/.streamlit/config.toml, etc.) --
# matching the OLD Dockerfile's own COPY/chmod paths (COPY . /app; chmod +x
# /app/bin/process_spatial.R), now confirmed against the real checkout.
#
# One deliberate deviation from the old Dockerfile: its CMD launched
# Streamlit as `streamlit run app/app.py` from the REPO ROOT (no `cd` into
# app/ first). This project's own confirmed, HPC-tested fix for a real
# Streamlit upload-limit bug established that `.streamlit/config.toml` is
# ONLY ever read relative to the directory `streamlit run` is invoked
# FROM, and that the working invocation is `cd repo/app && streamlit run
# app.py` (see DEPLOYMENT.md). The old CMD predates that fix and likely
# carried the identical latent bug -- any config.toml at
# repo/app/.streamlit/ would have been silently ignored under a
# "run from repo root" invocation, quite possibly unnoticed simply because
# no single upload through that container ever happened to exceed
# Streamlit's default 200MB limit. This Dockerfile intentionally uses the
# CONFIRMED-working "cd into app/" convention below instead.
#
# BUILD:
#   docker build -t bioportal:latest .
#
# RUN (mounts a persistent projects/reference-cache volume so projects and
# shared reference genomes/indices survive container restarts -- see
# project_manager.py's PROJECTS_ROOT / SHARED_REFERENCES_ROOT):
#   docker run -p 8501:8501 -v "$(pwd)/data:/workspace/repo/app/data" bioportal:latest

FROM mambaorg/micromamba:1.5.8 AS base

# Build the conda/mamba environment from the single-source-of-truth spec.
# Copying ONLY environment.yml first (before the rest of the source) means
# Docker's layer cache is reused for this expensive step on every rebuild
# that doesn't change dependencies -- only touches app code afterward.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# One-time Chrome pre-fetch for kaleido v1+ (see this file's header comment)
# -- done here, at BUILD time while network access is available, rather
# than deferred to the container's first real PDF-export use at runtime.
RUN python -c "import kaleido; kaleido.get_chrome_sync()"

# --- Application code ---
# See this file's "Repo layout (CONFIRMED)" comment above.
WORKDIR /workspace
COPY --chown=$MAMBA_USER:$MAMBA_USER repo/ /workspace/repo/

# The Spatial Transcriptomics pipeline's R script (invoked by main.nf) must
# be executable and reachable on PATH, exactly as the old Dockerfile
# ensured -- restored here since it's an application asset the old file
# handled that had no equivalent step in this file's earlier draft.
RUN chmod +x /workspace/repo/bin/process_spatial.R
ENV PATH="/workspace/repo/bin:${PATH}"

# Streamlit reads .streamlit/config.toml relative to the directory `streamlit
# run` is invoked FROM -- WORKDIR must match that exact directory (repo/app),
# not just contain it, per this project's own confirmed HPC testing (see
# this file's header comment on the repo layout assumption).
WORKDIR /workspace/repo/app

EXPOSE 8501

# --server.maxUploadSize is passed explicitly here too (belt-and-suspenders
# alongside .streamlit/config.toml) since it's cheap insurance against the
# exact config-discovery failure mode described above.
ENTRYPOINT ["streamlit", "run", "app.py", \
            "--server.address=0.0.0.0", \
            "--server.port=8501", \
            "--server.maxUploadSize=2048"]
