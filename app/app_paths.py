"""
app_paths.py

The single source of truth for WHERE this app's on-disk state lives.

--- The problem this fixes ---

Every storage-location constant in this codebase was written as a
RELATIVE path:

    AUTH_DIR              = os.path.join("data", "auth")     # auth_manager.py
    PROJECTS_ROOT         = "data/projects"                  # project_manager.py
    SC_PROJECTS_ROOT      = "data/singlecell_projects"       # sc_project_manager.py
    SHARED_REFERENCES_ROOT= "data/shared_references"         # project_manager.py
    SHARED_WHITELISTS_ROOT= "data/shared_whitelists"         # whitelist_manager.py
    CONNECTIONS_PATH      = "data/hpc_connections.json"      # hpc_manager.py
    INSTALL_STATUS_PATH   = "data/setup_install_status.json" # deployment_manager.py

A relative path resolves against the PROCESS'S CURRENT WORKING
DIRECTORY -- not against the location of the source file that declares
it. So all of the above silently mean "a data/ folder inside whatever
directory the app happened to be launched from."

That is fine exactly as long as the app is ALWAYS launched via
`cd repo/app && streamlit run app.py`. The moment it isn't -- a
different tmux pane, a systemd unit, a cron entry, an HPC job script,
or a user typing `streamlit run app/app.py` from the repo root -- every
path above silently points at a DIFFERENT, empty location, and the app
behaves as though it is a brand-new install.

This project has already been bitten by this exact class of bug once,
with .streamlit/config.toml: the file at repo/.streamlit/config.toml
was never loaded, because Streamlit also resolves it against the launch
directory, and it only started working once it was moved to
repo/app/.streamlit/. Same root cause; the consequences here are worse.

--- Why this is anchored to __file__ ---

os.path.abspath(__file__) is the absolute path of THIS source file,
which does not change based on how or from where the process was
started. Anchoring to its directory makes every data path deterministic
and launch-location-independent.

run_monitor_headless.py already uses exactly this pattern
(sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))) -- the
correct approach was already in the codebase, it just never made it
into the storage-path constants.
"""
import os

# Absolute path to the app/ directory -- i.e. the directory containing
# this file. Independent of the current working directory.
APP_ROOT = os.path.dirname(os.path.abspath(__file__))

# The single data root every other module should build on.
DATA_ROOT = os.path.join(APP_ROOT, "data")


def data_path(*parts):
    """
    Build an absolute path inside this app's data/ directory.

        data_path("auth", "users.json")
        -> /abs/path/to/repo/app/data/auth/users.json

    Use this instead of a bare "data/..." string literal anywhere that
    state is read from or written to disk.
    """
    return os.path.join(DATA_ROOT, *parts)


def describe_data_root():
    """
    Human-readable summary of where this app is actually reading and
    writing state, for display on the Setup & Deployment page.

    Worth surfacing in the UI: the whole point of this module is that
    the answer is no longer "wherever you happened to launch from," and
    being able to SEE the resolved path makes a misconfigured
    deployment obvious immediately instead of presenting as
    mysteriously-missing projects.
    """
    return {
        "app_root": APP_ROOT,
        "data_root": DATA_ROOT,
        "cwd": os.getcwd(),
        "cwd_matches_app_root": os.path.abspath(os.getcwd()) == APP_ROOT,
        "data_root_exists": os.path.isdir(DATA_ROOT),
    }


def find_legacy_data_root():
    """
    Detect a pre-patch data/ directory that was created relative to a
    DIFFERENT launch directory than app/, so the UI can tell a user
    their existing projects are in a stale location rather than
    silently presenting an empty app.

    Returns the absolute path of a legacy data/ directory if one is
    found that is NOT the canonical DATA_ROOT, else None.

    Only checks the two realistic cases: the current working directory,
    and the repository root one level above app/. This is a diagnostic
    helper only -- it deliberately does NOT move or copy anything.
    """
    candidates = [
        os.path.abspath(os.path.join(os.getcwd(), "data")),
        os.path.abspath(os.path.join(APP_ROOT, os.pardir, "data")),
    ]
    for candidate in candidates:
        if candidate == DATA_ROOT:
            continue
        if not os.path.isdir(candidate):
            continue
        # Only report it if it actually looks like OUR data directory,
        # rather than any unrelated folder that happens to be named
        # "data" (e.g. a user's own sequencing data folder).
        markers = ("projects", "singlecell_projects", "auth",
                   "shared_references", "shared_whitelists")
        if any(os.path.isdir(os.path.join(candidate, m)) for m in markers):
            return candidate
    return None
