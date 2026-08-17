"""
deployment_manager.py

Launches and tracks background installs of missing dependencies (Python
packages, CLI tools, R/Bioconductor packages) directly from the "⚙️ Setup &
Deployment" page's Environment & Dependency Check, using the SAME
subprocess.Popen + JSON status-file polling pattern already used elsewhere
in this app for long-running background work (see
advanced_mode_orchestrator.py's background pipeline runs) -- rather than
blocking Streamlit's single-threaded script execution for the several
minutes a real conda/mamba dependency-solve + download can take.

--- Why everything installs through conda/mamba, not pip/BiocManager ---
environment.yml (this project's single source of truth for dependencies)
declares EVERY dependency -- Python packages, CLI tools, AND R/Bioconductor
packages -- as conda packages (conda-forge/bioconda channels). Installing
through the same conda/mamba environment this app is already running in
means this feature can never install a different version than what
environment.yml pins, and avoids maintaining a second, separate pip/
BiocManager-based install path that could drift out of sync with it.

--- Important caveats surfaced to the user in setup_workspace.py's UI ---
- This modifies the ACTIVE conda/mamba environment the app itself is
  running in, live. Installing a package the running app depends on (e.g.
  streamlit, pandas) will not take effect for the CURRENT session --
  restarting the app after a successful install is recommended, and this
  is surfaced explicitly in the UI before the install button is offered.
- Requires the app to actually be running inside a conda/mamba environment
  (CONDA_PREFIX or CONDA_DEFAULT_ENV set) with mamba or conda on PATH --
  see get_active_conda_target(). If either is missing, install is refused
  with a specific, actionable reason rather than silently doing nothing.
"""
import json
import os
import shutil
import subprocess
import threading
from datetime import datetime

INSTALL_STATUS_PATH = "data/setup_install_status.json"
INSTALL_LOG_PATH = "data/setup_install_log.txt"

# --- Atomic-batch-failure design fix (2026-08-17) ---
# A real install attempt bundled ~8 missing packages (including paramiko,
# which resolves fine on its own) into ONE mamba install command. One of
# the OTHER packages in that batch named a version that doesn't exist on
# bioconda, which made mamba fail the ENTIRE solve atomically -- so
# paramiko (and everything else that would have installed correctly)
# never got installed either, even though nothing was wrong with it. A
# single bad spec silently blocking every other, perfectly-good package
# in the same batch is a real design gap, independent of whether any
# particular version pin is itself correct (see environment.yml's own
# "Version-pinning policy correction" for the pinning half of this fix).
#
# _run_install_with_fallback() below tries the full batch first (the fast
# path when nothing is wrong), and if THAT fails, automatically retries
# every package INDIVIDUALLY so one persistently-unresolvable spec can
# only ever block itself, not its batch-mates. Each package's outcome
# (installed vs. failed, with its own error) is tracked separately once
# this fallback path is taken, so the final status/log clearly shows
# exactly which packages succeeded and which didn't, rather than a single
# opaque "the batch failed" message that leaves every package's real
# status a mystery.



def _ensure_parent_dir():
    for path in (INSTALL_STATUS_PATH, INSTALL_LOG_PATH):
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)


def get_conda_or_mamba_executable():
    "Prefer mamba (faster dependency solving) over conda, falling back to conda if mamba isn't on PATH. Returns None if neither is found."
    return shutil.which("mamba") or shutil.which("conda")


def get_active_conda_target():
    """
    Detect which conda/mamba environment this app is currently running in,
    returning (flag, value) e.g. ("-p", "/path/to/env") or ("-n", "envname")
    suitable for splicing directly into an install command's argument list,
    or None if this process doesn't appear to be running inside a conda/
    mamba environment at all (in which case there'd be no way to know
    which environment an install should even target).
    """
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        return ("-p", prefix)
    name = os.environ.get("CONDA_DEFAULT_ENV")
    if name:
        return ("-n", name)
    return None


# --- Missing-channel fix (2026-08-17) ---
# A real install attempt failed on every Bioconductor/bioconda-only package
# (tximport, every org.*.db package) while a conda-forge-only package
# (paramiko) installed fine in the same run -- and the install log showed
# ONLY "conda-forge/linux-64"/"conda-forge/noarch" being fetched/cached;
# bioconda never appeared in the log at all. environment.yml's own
# "channels:" list (conda-forge, bioconda, defaults) is only automatically
# applied to an environment's channel config when that environment is
# CREATED via `mamba env create -f environment.yml` -- it has no effect on
# an ad-hoc `mamba install -p <path> <specs>` call run against an
# environment that was set up some other way (as was the case here), so
# that install silently only searched whatever channels were already
# configured for that specific environment, missing bioconda entirely.
# Every install command built below now explicitly passes -c for each
# channel, in the SAME priority order as environment.yml's own channels
# list (first listed = highest priority for both conda's "channels:" key
# and repeated -c flags), so an install triggered from this page can never
# again depend on the target environment already having the right
# channels configured -- it's self-contained regardless of how that
# environment was originally created.
INSTALL_CHANNELS = ("conda-forge", "bioconda", "defaults")


def _channel_args():
    args = []
    for channel in INSTALL_CHANNELS:
        args += ["-c", channel]
    return args


def build_install_command(package_specs):
    "Build the full mamba/conda install command as an argv list, or None if mamba/conda isn't available, no active environment could be detected, or package_specs is empty."
    exe = get_conda_or_mamba_executable()
    target = get_active_conda_target()
    if not exe or not target or not package_specs:
        return None
    flag, value = target
    return [exe, "install", "-y", flag, value] + _channel_args() + list(package_specs)


def get_install_status():
    "Returns the current/most recent install's status dict, or None if no install has ever been launched."
    if not os.path.exists(INSTALL_STATUS_PATH):
        return None
    with open(INSTALL_STATUS_PATH) as f:
        return json.load(f)


def is_install_in_progress():
    status = get_install_status()
    return bool(status) and status.get("status") == "running"


def read_install_log(n_lines=300):
    "Return the last n_lines of the install log, or '' if it doesn't exist yet."
    if not os.path.exists(INSTALL_LOG_PATH):
        return ""
    with open(INSTALL_LOG_PATH) as f:
        lines = f.readlines()
    return "".join(lines[-n_lines:])


def _write_status(status_dict):
    _ensure_parent_dir()
    with open(INSTALL_STATUS_PATH, "w") as f:
        json.dump(status_dict, f, indent=2)


def _run_and_log(command, log_f):
    "Run command to completion, streaming its combined stdout/stderr into the already-open log_f. Returns the process's exit code. Blocks -- only ever called from within the background watcher thread, never on the main/UI thread."
    process = subprocess.Popen(command, stdout=log_f, stderr=subprocess.STDOUT)
    process.wait()
    return process.returncode


def launch_install(package_specs):
    """
    Launch a background mamba/conda install for package_specs (a list of
    conda package spec strings, e.g. ["fastqc=0.12.*", "bioconductor-deseq2"]).
    Returns (success: bool, message: str). On success, the install runs in
    a background thread; poll get_install_status() / read_install_log()
    for progress (see setup_workspace.py's polling UI) -- this function
    itself returns immediately rather than blocking for the install's full
    duration.

    Tries the full batch first (one solve, fast when nothing's wrong). If
    that fails, automatically falls back to installing each package
    INDIVIDUALLY so a single unresolvable spec can only ever block itself,
    not every other package bundled in the same batch -- see this module's
    "Atomic-batch-failure design fix" note above. The final status's
    "package_results" field records each package's individual outcome
    whenever the fallback path was used, so it's always clear exactly
    which packages installed and which didn't (and why), rather than one
    opaque batch-level failure.
    """
    if is_install_in_progress():
        return False, "An install is already in progress -- wait for it to finish before starting another."

    package_specs = list(package_specs)
    exe = get_conda_or_mamba_executable()
    target = get_active_conda_target()
    if not exe:
        return False, "Neither `mamba` nor `conda` was found on PATH -- cannot install automatically. See DEPLOYMENT.md for manual setup instructions."
    if not target:
        return False, "This app doesn't appear to be running inside a conda/mamba environment (no CONDA_PREFIX/CONDA_DEFAULT_ENV detected) -- cannot determine which environment to install into."
    if not package_specs:
        return False, "No packages were specified to install."

    batch_command = build_install_command(package_specs)
    _ensure_parent_dir()
    log_f = open(INSTALL_LOG_PATH, "w")
    log_f.write(f"$ {' '.join(batch_command)}\n\n")
    log_f.flush()

    _write_status({
        "status": "running",
        "pid": None,
        "command": batch_command,
        "package_specs": package_specs,
        "package_results": None,
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "finished_at": None,
        "returncode": None,
        "used_fallback": False,
    })

    def _watch():
        batch_returncode = _run_and_log(batch_command, log_f)
        if batch_returncode == 0:
            log_f.close()
            status = get_install_status() or {}
            status["status"] = "complete"
            status["finished_at"] = datetime.now().isoformat(timespec="seconds")
            status["returncode"] = 0
            _write_status(status)
            return

        # Full-batch solve failed -- fall back to one-at-a-time so a single
        # bad spec can't block every other, otherwise-installable package.
        log_f.write(
            f"\n--- Batch install failed (exit code {batch_returncode}). "
            "Falling back to installing each package individually so one "
            "bad spec doesn't block the rest. ---\n\n"
        )
        log_f.flush()
        flag, value = target
        package_results = {}
        any_failed = False
        for spec in package_specs:
            individual_command = [exe, "install", "-y", flag, value] + _channel_args() + [spec]
            log_f.write(f"$ {' '.join(individual_command)}\n")
            log_f.flush()
            rc = _run_and_log(individual_command, log_f)
            package_results[spec] = "installed" if rc == 0 else "failed"
            if rc != 0:
                any_failed = True
            log_f.write("\n")
            log_f.flush()
        log_f.close()

        status = get_install_status() or {}
        status["status"] = "error" if any_failed else "complete"
        status["finished_at"] = datetime.now().isoformat(timespec="seconds")
        status["returncode"] = 1 if any_failed else 0
        status["used_fallback"] = True
        status["package_results"] = package_results
        _write_status(status)

    thread = threading.Thread(target=_watch, daemon=True)
    thread.start()
    return True, f"Install started for {len(package_specs)} package(s)."
