# Dockerfile
#
# Builds the Pretty Awesome Transcriptomics Tool image entirely from
# environment.yml (this repo's single source of truth for Python,
# R/Bioconductor, and CLI tool dependencies) via micromamba -- so this
# image can never drift out of sync with what a local or HPC
# conda/mamba install would produce from the same environment.yml.
#
# BUILD (from the repo root -- the directory containing this file):
#   docker build -t bioportal:latest .
#
# RUN -- preferred:
#   docker compose up -d
# (see docker-compose.yml, which also handles the host port and the
# host/container UID mismatch described below)
#
# ===================================================================
# ON PORTS -- read this before changing anything here
# ===================================================================
# The app ALWAYS listens on 8501 INSIDE the container. That is correct
# even if another container on the same host is already using 8501,
# because each container has its own isolated network namespace -- two
# containers can both listen on 8501 internally with no conflict.
#
# To serve on a different HOST port, change ONLY the left-hand side of
# the ports mapping in docker-compose.yml:
#     ports:
#       - "8502:8501"      # host 8502 -> container 8501
#
# ...or set BIOPORTAL_HOST_PORT in your .env file. Do NOT change the
# port in this file -- EXPOSE is documentation only, and the
# --server.port below refers to the container's own namespace.
#
# ===================================================================
# FIXES APPLIED 2026-09-16 -- both of these prevented the image from
# ever building or starting. Documented rather than silently corrected,
# because both are easy to reintroduce.
# ===================================================================
#
# --- FIX 1: the COPY path did not exist in the build context ---
# The previous version did:
#     COPY --chown=... repo/ /workspace/repo/
# The build context is the REPO ROOT itself (`docker build .` run from
# the directory holding environment.yml and Dockerfile). That directory
# contains app/, bin/, environment.yml, main.nf, nextflow.config -- it
# does NOT contain a nested `repo/` subdirectory. The path `repo/` only
# ever existed as a name in this file's own explanatory comments, not on
# disk. Result: the build failed at this layer every time with
# "COPY failed: stat repo/: file does not exist". Corrected to `COPY .`,
# which lands app/ at /workspace/repo/app -- matching both the WORKDIR
# below and the documented volume mount.
#
# --- FIX 2: overriding ENTRYPOINT disabled the conda environment ---
# The previous version ended with:
#     ENTRYPOINT ["streamlit", "run", "app.py", ...]
# mambaorg/micromamba images ship their own
# ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"], and THAT script is the
# only thing that runs `micromamba activate` before exec'ing the command
# after it. Replacing it means the base environment is never activated at
# runtime, so none of its bin/ directory is on PATH and the container
# exits immediately with "executable file not found in $PATH".
#
# ARG MAMBA_DOCKERFILE_ACTIVATE=1 does NOT cover this -- it affects
# build-time RUN layers only and has no effect on the running container.
# (Which is exactly why the kaleido fetch below worked during the build
# while the runtime would still have failed.) Corrected by keeping
# _entrypoint.sh as argv[0] and passing streamlit to it.
#
# --- Repo layout (CONFIRMED) ---
# Repo root holds environment.yml, Dockerfile, bin/process_spatial.R,
# main.nf, nextflow.config, with the Streamlit app in an app/ subfolder
# (app/app.py, app/.streamlit/config.toml, app/single_cell/, ...).
#
# --- Why WORKDIR is repo/app and not repo/ ---
# Two independent reasons, both confirmed against this codebase:
#   1. Streamlit only reads .streamlit/config.toml relative to the
#      directory `streamlit run` is invoked FROM. This project already
#      hit that bug on HPC: a config at the repo root was silently
#      ignored and uploads stayed capped at Streamlit's 200MB default
#      until it was moved to app/.streamlit/.
#   2. app/app_paths.py anchors every data path to the app/ directory via
#      os.path.abspath(__file__), so DATA_ROOT resolves to
#      /workspace/repo/app/data regardless of working directory. The
#      volume mount must match that exact path or every project written
#      inside the container is lost on `docker compose down`.
#
# --- kaleido + Chrome ---
# python-kaleido v1+ never bundles Chrome; a one-time
# kaleido.get_chrome_sync() fetch is required. Done here at BUILD time
# while network access is available, rather than deferred to the
# container's first real PDF export at RUNTIME (where the container may
# have no outbound network and would fail silently on first use).
#
# --- DoubletFinder ---
# Installed below at build time. README.md and deployment_guide.md both
# already told users "the Docker image runs this automatically at build
# time (see Dockerfile)" -- but the previous Dockerfile never did. That
# claim is now true. It is a real RUN layer that will fail the build if
# GitHub is unreachable, deliberately: a silent skip would leave a user
# with an image whose documentation promises a method that isn't there.

FROM mambaorg/micromamba:1.5.8 AS base

# Build the environment from the single-source-of-truth spec.
# Copying ONLY environment.yml first means Docker's layer cache is reused
# for this expensive step (~20-40 min) on every rebuild that doesn't
# change dependencies -- app code changes below won't retrigger it.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activates the base env for the RUN layers BELOW this line (build time
# only -- see FIX 2 in this file's header for the runtime half).
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# One-time Chrome pre-fetch for kaleido v1+ (see header).
RUN python -c "import kaleido; kaleido.get_chrome_sync()"

# DoubletFinder: GitHub-only, no conda/CRAN/Bioconductor release at all.
# r-remotes is in environment.yml specifically to make this possible.
RUN Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never")'

# --- Application code ---
WORKDIR /workspace
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace/repo/

# The Spatial Transcriptomics pipeline's R script (invoked by main.nf)
# must be executable and reachable on PATH.
RUN chmod +x /workspace/repo/bin/process_spatial.R
ENV PATH="/workspace/repo/bin:${PATH}"

# Create the data directory inside the image so a bind mount from the
# host inherits a sane owner, and so the container still starts cleanly
# when run with NO mount at all (projects are then ephemeral, which is
# fine for a smoke test).
RUN mkdir -p /workspace/repo/app/data

# Must MATCH the directory streamlit is launched from -- not merely
# contain it. See "Why WORKDIR is repo/app" in the header.
WORKDIR /workspace/repo/app

# Documentation only -- does not publish anything. The host-side port is
# set in docker-compose.yml. See "ON PORTS" in this file's header.
EXPOSE 8501

# Fails fast and visibly if Streamlit stops serving, rather than leaving
# a container that looks "up" but answers nothing. Uses localhost inside
# the container, so it is unaffected by whatever host port you map to.
HEALTHCHECK --interval=30s --timeout=5s --start-period=40s --retries=3 \
    CMD curl -fsS http://localhost:8501/_stcore/health || exit 1

# _entrypoint.sh MUST stay as argv[0] -- it is what activates the conda
# environment at runtime. See FIX 2 in this file's header.
#
# --server.maxUploadSize is passed explicitly as belt-and-suspenders
# alongside app/.streamlit/config.toml.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", \
            "streamlit", "run", "app.py", \
            "--server.address=0.0.0.0", \
            "--server.port=8501", \
            "--server.maxUploadSize=2048"]
