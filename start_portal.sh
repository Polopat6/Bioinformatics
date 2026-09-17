#!/usr/bin/env bash
#
# start_portal.sh
#
# One-command startup for the Pretty Awesome Transcriptomics Tool on any
# machine: automatically wraps itself in a persistent tmux session (if not
# already running inside one), activates the project's conda environment,
# and launches Streamlit -- collapsing what has repeatedly been a
# multi-step manual dance (tmux new/attach -> conda activate -> cd into
# app/ -> streamlit run, each redone by hand after every SSH disconnect)
# into a single command.
#
# IMPORTANT: this script does NOT make the file-browsing feature (see
# file_browser.py) itself need tmux -- once Streamlit is running inside a
# tmux session, the server-side process is already fully decoupled from
# any one SSH connection; closing a laptop, losing an SSH tunnel, or a
# browser tab crashing has zero effect on the running Streamlit process.
# This script's purpose is just to make STARTING that tmux-protected
# process reliably a single command.
#
# NOTE: this is the NON-Docker path (local machine or HPC). If you are
# running via `docker compose up -d`, you do not need this script at all.
#
# --- Made portable 2026-09-16 ---
# The previous version hardcoded one specific HPC checkout:
#     CONDA_ENV_ROOT="/scratch/bioscratch/Podrab_lab/.../miniforge3"
#     APP_DIR="/scratch/bioscratch/Podrab_lab/.../repo/app"
#     CONDA_ENV_NAME="rnaseq"
# ...so it could not run on any other machine without editing the file,
# and it named an environment ("rnaseq") that environment.yml does not
# create -- environment.yml declares `name: bioportal`. Anyone following
# README.md exactly would create `bioportal`, then find this script
# looking for something else entirely.
#
# Now: APP_DIR is derived from this script's own location, the conda
# installation is auto-detected, and every value can be overridden by an
# environment variable without editing the file.
#
# Usage:
#   ./start_portal.sh                    # normal start (auto-wraps in tmux)
#   ./start_portal.sh --attach           # reattach to a running portal
#   ./start_portal.sh --no-tmux          # run in the foreground (systemd, debugging)
#   BIOPORTAL_ENV=rnaseq ./start_portal.sh          # different env name
#   BIOPORTAL_PORT=8502 ./start_portal.sh           # different port
#   BIOPORTAL_CONDA_ROOT=~/miniforge3 ./start_portal.sh

set -euo pipefail

# --- Configuration (every value overridable via environment variable) ---
TMUX_SESSION_NAME="${BIOPORTAL_SESSION:-bioportal}"

# Must match environment.yml's own `name:` field. Override with
# BIOPORTAL_ENV to point at a differently-named existing environment.
CONDA_ENV_NAME="${BIOPORTAL_ENV:-bioportal}"

# Host port for this (non-Docker) launch. Unrelated to the Docker host
# port in docker-compose.yml -- set BIOPORTAL_PORT if 8501 is already in
# use on this machine (e.g. by a running container).
STREAMLIT_PORT="${BIOPORTAL_PORT:-8501}"
STREAMLIT_MAX_UPLOAD_MB="${BIOPORTAL_MAX_UPLOAD_MB:-2048}"

# Derived from this script's own location rather than hardcoded, so the
# repo can be cloned anywhere. Resolves to <repo>/app given that this
# script sits at the repo root.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="${BIOPORTAL_APP_DIR:-${SCRIPT_DIR}/app}"


# --- Locate a conda/mamba installation ---
# Checked in order: explicit override, already-active conda, then the
# common install locations for miniforge/mambaforge/miniconda/anaconda.
find_conda_root() {
    if [[ -n "${BIOPORTAL_CONDA_ROOT:-}" ]]; then
        echo "${BIOPORTAL_CONDA_ROOT}"; return 0
    fi
    # If a conda env is already active, reuse that installation.
    if [[ -n "${CONDA_EXE:-}" ]]; then
        dirname "$(dirname "${CONDA_EXE}")"; return 0
    fi
    if [[ -n "${MAMBA_EXE:-}" ]]; then
        dirname "$(dirname "${MAMBA_EXE}")"; return 0
    fi
    local candidate
    for candidate in \
        "${HOME}/miniforge3" "${HOME}/mambaforge" \
        "${HOME}/miniconda3" "${HOME}/anaconda3" \
        "/opt/miniforge3" "/opt/conda" "/opt/miniconda3"
    do
        if [[ -f "${candidate}/etc/profile.d/conda.sh" ]]; then
            echo "${candidate}"; return 0
        fi
    done
    return 1
}


activate_env_and_run() {
    local conda_root
    if ! conda_root="$(find_conda_root)"; then
        echo "ERROR: could not find a conda/mamba installation." >&2
        echo "       Set BIOPORTAL_CONDA_ROOT to its path, e.g.:" >&2
        echo "         BIOPORTAL_CONDA_ROOT=\$HOME/miniforge3 $0" >&2
        exit 1
    fi

    # shellcheck disable=SC1091
    source "${conda_root}/etc/profile.d/conda.sh"

    if ! conda activate "${CONDA_ENV_NAME}" 2>/dev/null; then
        echo "ERROR: conda environment '${CONDA_ENV_NAME}' not found in ${conda_root}." >&2
        echo "       Create it from the repo root with:" >&2
        echo "         mamba env create -f environment.yml" >&2
        echo "       ...or point at an existing one:" >&2
        echo "         BIOPORTAL_ENV=<your_env_name> $0" >&2
        exit 1
    fi

    # NOTE on `hash -r`: after a conda environment is activated for the
    # first time in a given shell, bash sometimes still resolves commands
    # like `python3` to a PREVIOUSLY cached location (e.g. the system
    # Python) rather than the just-activated environment's own binaries,
    # until its internal command-location cache is cleared. This was hit
    # directly during real deployment testing -- `which python3` and
    # `sys.executable` both incorrectly pointed at /usr/bin/python3
    # despite `conda activate` having correctly updated $PATH -- and
    # `hash -r` was the fix. Kept defensively so this script never
    # reproduces that confusing failure mode for someone running it fresh.
    hash -r

    if [[ ! -d "${APP_DIR}" ]]; then
        echo "ERROR: app directory not found: ${APP_DIR}" >&2
        echo "       Set BIOPORTAL_APP_DIR if your layout differs." >&2
        exit 1
    fi

    # `cd` into app/ before launching is REQUIRED, not stylistic:
    #   1. Streamlit only reads .streamlit/config.toml relative to the
    #      directory it is invoked FROM. This project already hit that bug
    #      on HPC -- a config at the repo root was silently ignored.
    #   2. It keeps the working directory consistent with what
    #      app_paths.py resolves, which makes the data root obvious when
    #      debugging.
    cd "${APP_DIR}"

    echo "Environment  : ${CONDA_ENV_NAME}  (${conda_root})"
    echo "App directory: ${APP_DIR}"
    echo "Starting Streamlit on port ${STREAMLIT_PORT}..."
    exec streamlit run app.py \
        --server.port "${STREAMLIT_PORT}" \
        --server.headless true \
        --server.maxUploadSize "${STREAMLIT_MAX_UPLOAD_MB}"
}


# --- Argument handling ---
case "${1:-}" in
    --attach)
        if tmux has-session -t "${TMUX_SESSION_NAME}" 2>/dev/null; then
            exec tmux attach -t "${TMUX_SESSION_NAME}"
        fi
        echo "No existing '${TMUX_SESSION_NAME}' tmux session found -- run without --attach to start one."
        exit 1
        ;;
    --no-tmux)
        # Foreground mode, for systemd units or interactive debugging
        # where an extra tmux layer is unwanted.
        activate_env_and_run
        ;;
    "")
        ;;
    *)
        echo "Usage: $0 [--attach | --no-tmux]" >&2
        exit 2
        ;;
esac


# --- Already inside tmux: do the real work directly ---
# ($TMUX is set by tmux itself for any process running inside a tmux
# pane -- the standard way to detect this from a shell script.)
if [[ -n "${TMUX:-}" ]]; then
    echo "Already inside a tmux session -- starting the portal directly."
    activate_env_and_run
fi

# --- Not inside tmux yet ---
if ! command -v tmux >/dev/null 2>&1; then
    echo "NOTE: tmux is not installed -- starting in the foreground instead."
    echo "      The portal will stop if this terminal closes. Install tmux,"
    echo "      or use Docker (see docker-compose.yml) for a persistent service."
    activate_env_and_run
fi

if tmux has-session -t "${TMUX_SESSION_NAME}" 2>/dev/null; then
    echo "An existing '${TMUX_SESSION_NAME}' session is already running -- reattaching."
    echo "(No new Streamlit process will be started.)"
    exec tmux attach -t "${TMUX_SESSION_NAME}"
fi

echo "Starting a new persistent tmux session '${TMUX_SESSION_NAME}'..."
# Re-invoke this exact script (absolute path, so it works regardless of
# the caller's working directory) inside a new tmux session -- that inner
# invocation sees $TMUX set and takes the branch above.
exec tmux new-session -s "${TMUX_SESSION_NAME}" "${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
