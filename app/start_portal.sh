#!/usr/bin/env bash
#
# start_portal.sh
#
# One-command startup for the bioinformatics portal on a shared server /
# HPC node: automatically wraps itself in a persistent tmux session (if
# not already running inside one), activates the project's conda
# environment, and launches Streamlit -- collapsing what has repeatedly
# been a multi-step manual dance (tmux new/attach -> conda activate ->
# cd into app/ -> streamlit run, each redone by hand after every SSH
# disconnect) into a single command.
#
# IMPORTANT: this script does NOT make the file-browsing feature (see
# file_browser.py) itself need tmux -- once Streamlit is running inside
# a tmux session, the server-side process is already fully decoupled
# from any one SSH connection; closing a laptop, losing an SSH tunnel,
# or a browser tab crashing has zero effect on the running Streamlit
# process. This script's actual purpose is just to make STARTING that
# tmux-protected process reliably a single command, rather than
# something that has to be manually re-assembled from memory every time.
#
# Usage:
#   ./start_portal.sh              # normal start (auto-wraps in tmux)
#   ./start_portal.sh --attach     # if already running, just reattach
#
# Configuration (edit these three lines for your own deployment):
TMUX_SESSION_NAME="bioportal"
CONDA_ENV_ROOT="/scratch/bioscratch/Podrab_lab/Pate_is_grate/bioproject/miniforge3"
CONDA_ENV_NAME="rnaseq"
APP_DIR="/scratch/bioscratch/Podrab_lab/Pate_is_grate/bioproject/repo/app"
STREAMLIT_PORT="8501"

set -euo pipefail

# --- Handle --attach: just reattach to an existing session, don't relaunch ---
if [[ "${1:-}" == "--attach" ]]; then
    if tmux has-session -t "$TMUX_SESSION_NAME" 2>/dev/null; then
        exec tmux attach -t "$TMUX_SESSION_NAME"
    else
        echo "No existing '$TMUX_SESSION_NAME' tmux session found -- run without --attach to start a fresh one."
        exit 1
    fi
fi

# --- If we're already inside tmux, just do the actual work directly ---
# ($TMUX is set by tmux itself for any process running inside a tmux
# pane -- this is the standard, documented way to detect "am I already
# inside tmux" from a shell script.)
if [[ -n "${TMUX:-}" ]]; then
    echo "Already inside a tmux session -- starting the portal directly."
    source "${CONDA_ENV_ROOT}/bin/activate" "$CONDA_ENV_NAME"
    hash -r  # clear bash's cached command locations -- see note below
    cd "$APP_DIR"
    echo "Starting Streamlit on port ${STREAMLIT_PORT}..."
    exec streamlit run app.py --server.port "$STREAMLIT_PORT" --server.headless true
fi

# --- Not inside tmux yet: either attach to an existing session, or
# create a new one that re-invokes THIS SAME SCRIPT (which will then
# hit the "already inside tmux" branch above and do the real work) ---
if tmux has-session -t "$TMUX_SESSION_NAME" 2>/dev/null; then
    echo "An existing '$TMUX_SESSION_NAME' session is already running -- reattaching to it"
    echo "(use this to check on / interact with an already-running portal;"
    echo " no new Streamlit process will be started)."
    exec tmux attach -t "$TMUX_SESSION_NAME"
else
    echo "Starting a new persistent tmux session '$TMUX_SESSION_NAME' and launching the portal inside it..."
    # Re-invoke this exact script (with its absolute path, so this works
    # regardless of the caller's current working directory) inside a
    # brand-new tmux session -- that inner invocation will see $TMUX set
    # and fall into the "already inside tmux" branch above, which does
    # the actual conda activate + streamlit run.
    SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    exec tmux new-session -s "$TMUX_SESSION_NAME" "$SCRIPT_PATH"
fi

# NOTE on `hash -r` above: after a conda environment is activated for
# the first time in a given shell, bash sometimes still resolves
# commands like `python3` to a PREVIOUSLY cached location (e.g. the
# system Python) rather than the just-activated environment's own
# binaries, until its internal command-location cache is cleared. This
# was hit directly during real deployment testing -- `which python3`
# and `sys.executable` both incorrectly pointed at /usr/bin/python3
# despite `conda activate` having correctly updated $PATH -- and
# `hash -r` (clear bash's cached command locations, forcing it to
# re-resolve every command against the CURRENT $PATH) was the fix.
# Included here defensively so this script never reproduces that same
# confusing failure mode for someone else running it fresh.
