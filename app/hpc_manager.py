"""
hpc_manager.py

Manages saved SSH connection profiles to remote HPC clusters, and tests
connectivity/detects the remote environment (OS, scheduler, conda/mamba
presence) for a given profile. Used by setup_workspace.py's "HPC
Connections" section.

--- Security design (read this before changing anything here) ---
Connection PROFILES (name, host, port, username, auth method, key path)
are persisted to disk as plain JSON -- none of that is a secret. A
PASSWORD, if a user chooses password-based auth for a one-off test, is
NEVER persisted anywhere: it only ever exists in a Streamlit widget's
session_state for the current session and is passed directly into
paramiko's connect() call, never written to CONNECTIONS_PATH or logged.
save_connection() defensively strips any "password" key that might
accidentally be passed in, as a second line of defense beyond simply
"remembering not to include it".

SSH-key-based (or ssh-agent-based) authentication is the strongly
recommended path for any SAVED connection meant to be reused -- see
setup_workspace.py's UI copy, which surfaces this same guidance to the
user directly.

Host key verification: test_connection() uses paramiko's AutoAddPolicy,
which trusts an unknown host key on first connect (equivalent to SSH's
interactive "yes" prompt for an unrecognized host) rather than either
rejecting unknown hosts outright or silently ignoring host key changes.
This is a reasonable default for an interactive "test this connection"
button a human is actively driving, but is NOT equivalent to strict
known_hosts verification -- a future hardening pass could offer a
"use my local ~/.ssh/known_hosts strictly" toggle via
paramiko.RejectPolicy + load_system_host_keys() for users who want that.

--- paramiko availability ---
paramiko is listed in environment.yml, but this module is written to
degrade gracefully (PARAMIKO_AVAILABLE check) rather than raising an
ImportError at import time, so the rest of the app still works even in
an environment where it genuinely isn't installed yet -- e.g. someone's
existing environment that hasn't been re-synced against environment.yml.
"""
import json
import os
from datetime import datetime

try:
    import paramiko
    PARAMIKO_AVAILABLE = True
except ImportError:
    paramiko = None
    PARAMIKO_AVAILABLE = False

CONNECTIONS_PATH = "data/hpc_connections.json"

AUTH_METHODS = ("key", "agent", "password")

# Read-only remote commands used to build a quick "what does this cluster
# look like" summary after a successful test connection -- deliberately
# non-invasive (no writes, no job submission) since this is a connectivity
# check, not a remote execution feature (that would be a separate,
# explicitly-scoped future addition).
_REMOTE_PROBE_COMMANDS = {
    "os_info": "uname -a",
    "conda_available": "command -v conda || command -v mamba || echo NONE",
    "slurm_available": "command -v sbatch || echo NONE",
    "pbs_available": "command -v qsub || echo NONE",
}


def _ensure_parent_dir():
    parent = os.path.dirname(CONNECTIONS_PATH)
    if parent:
        os.makedirs(parent, exist_ok=True)


def list_connections():
    """
    Return a list of saved connection profile dicts (never containing a
    password -- see module docstring), sorted by profile_name.
    """
    if not os.path.exists(CONNECTIONS_PATH):
        return []
    with open(CONNECTIONS_PATH) as f:
        data = json.load(f)
    return sorted(data.values(), key=lambda c: c.get("profile_name", ""))


def get_connection(profile_name):
    if not os.path.exists(CONNECTIONS_PATH):
        return None
    with open(CONNECTIONS_PATH) as f:
        data = json.load(f)
    return data.get(profile_name)


def connection_exists(profile_name):
    return get_connection(profile_name) is not None


def save_connection(profile):
    """
    Persist a connection profile to disk. `profile` must include at least
    profile_name/host/port/username/auth_method. Any "password" key is
    stripped before writing, regardless of how it got there -- passwords
    are never persisted by this module, full stop.
    """
    profile = dict(profile)
    profile.pop("password", None)
    profile["last_updated"] = datetime.now().isoformat(timespec="seconds")
    _ensure_parent_dir()
    data = {}
    if os.path.exists(CONNECTIONS_PATH):
        with open(CONNECTIONS_PATH) as f:
            data = json.load(f)
    data[profile["profile_name"]] = profile
    with open(CONNECTIONS_PATH, "w") as f:
        json.dump(data, f, indent=2)


def delete_connection(profile_name):
    if not os.path.exists(CONNECTIONS_PATH):
        return False
    with open(CONNECTIONS_PATH) as f:
        data = json.load(f)
    if profile_name not in data:
        return False
    del data[profile_name]
    with open(CONNECTIONS_PATH, "w") as f:
        json.dump(data, f, indent=2)
    return True


def record_test_result(profile_name, success, message):
    """Update a saved profile's last-tested status, shown in the UI so a stale 'last known good' state is never displayed as current without being labeled with when it was actually checked."""
    profile = get_connection(profile_name)
    if not profile:
        return
    profile["last_tested_at"] = datetime.now().isoformat(timespec="seconds")
    profile["last_test_success"] = success
    profile["last_test_message"] = message
    save_connection(profile)


def test_connection(host, port, username, auth_method, key_path=None, password=None, timeout=10):
    """
    Attempt an SSH connection and, if successful, run a small set of
    read-only probe commands to summarize the remote environment.

    Returns (success: bool, message: str, remote_info: dict | None).
    remote_info is None on failure, otherwise a dict with keys matching
    _REMOTE_PROBE_COMMANDS plus a human-readable "scheduler" summary
    derived from the slurm/pbs probes.
    """
    if not PARAMIKO_AVAILABLE:
        return False, (
            "paramiko isn't installed in this environment yet -- it's listed "
            "in environment.yml, so re-syncing your conda/mamba environment "
            "against that file (see DEPLOYMENT.md) will add it."
        ), None

    if auth_method not in AUTH_METHODS:
        return False, f"Unknown auth method '{auth_method}'.", None
    if auth_method == "key" and not key_path:
        return False, "SSH key file auth selected, but no key path was provided.", None
    if auth_method == "key" and not os.path.isfile(key_path):
        return False, f"SSH key file not found at: {key_path}", None
    if auth_method == "password" and not password:
        return False, "Password auth selected, but no password was provided.", None

    client = paramiko.SSHClient()
    # See module docstring's "Host key verification" note.
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    connect_kwargs = {
        "hostname": host,
        "port": int(port),
        "username": username,
        "timeout": timeout,
        "banner_timeout": timeout,
        "auth_timeout": timeout,
    }
    if auth_method == "key":
        connect_kwargs["key_filename"] = key_path
        connect_kwargs["allow_agent"] = False
        connect_kwargs["look_for_keys"] = False
    elif auth_method == "agent":
        connect_kwargs["allow_agent"] = True
        connect_kwargs["look_for_keys"] = True
    else:  # password
        connect_kwargs["password"] = password
        connect_kwargs["allow_agent"] = False
        connect_kwargs["look_for_keys"] = False

    try:
        client.connect(**connect_kwargs)
    except paramiko.AuthenticationException:
        return False, "Authentication failed -- check username/key/password and that this key is authorized on the remote host.", None
    except paramiko.SSHException as e:
        return False, f"SSH negotiation failed: {e}", None
    except (OSError, TimeoutError) as e:
        return False, f"Could not reach {host}:{port} -- {e}", None
    except Exception as e:  # noqa: BLE001 -- surface any unexpected failure to the UI rather than crashing the page
        return False, f"Unexpected error connecting: {e}", None

    remote_info = {}
    try:
        for key, cmd in _REMOTE_PROBE_COMMANDS.items():
            _stdin, stdout, _stderr = client.exec_command(cmd, timeout=timeout)
            remote_info[key] = stdout.read().decode(errors="replace").strip()
    finally:
        client.close()

    schedulers = []
    if remote_info.get("slurm_available") and remote_info["slurm_available"] != "NONE":
        schedulers.append("SLURM")
    if remote_info.get("pbs_available") and remote_info["pbs_available"] != "NONE":
        schedulers.append("PBS/Torque")
    remote_info["scheduler"] = ", ".join(schedulers) if schedulers else "None detected"
    remote_info["conda_or_mamba"] = (
        "Found" if remote_info.get("conda_available") and remote_info["conda_available"] != "NONE" else "Not found on PATH"
    )

    return True, f"Connected to {host} as {username}.", remote_info
