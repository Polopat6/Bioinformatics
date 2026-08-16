"""
monitor_daemon_runner.py

Standalone entry point for ONE Monitor Mode directory-watcher daemon.
Launched by monitor_manager.launch_monitor_daemon(monitor_id) as a
fully detached OS subprocess (NOT imported/run inside the Streamlit app
process) -- mirrors sra_download_runner.py's role for the SRA-download
background job, just for a long-lived watch loop instead of a one-shot
task.

Since monitor_manager.py now supports MULTIPLE independent monitors,
each monitor gets its OWN instance of this script running as its own
OS process, pointed at that specific monitor's own config.json (see
monitor_manager.monitor_config_path(monitor_id)) -- there is no shared
global state between two monitors' daemon processes; each is
completely independent (its own registry, its own activity log, its
own PID file) and can be started/stopped without affecting any other
monitor.

IMPORTANT: this script deliberately does NOT import streamlit anywhere,
directly or indirectly. It runs completely detached from any Streamlit
script-run context, so any st.* call here would either silently no-op
in confusing ways or raise "missing ScriptRunContext" warnings/errors.
The watcher's own status is communicated purely through the plain
JSON/JSONL files monitor_manager.py already reads/writes (this
monitor's own registry + activity log) -- see monitor_mode_workspace.py,
which polls these files from a *separate*, later Streamlit session.

Usage (invoked exactly this way by monitor_manager.launch_monitor_daemon):
    python monitor_daemon_runner.py <this_monitor's_config_json_path>
"""
import json
import sys

import monitor_manager as mm


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <monitor_config_json_path>", file=sys.stderr)
        sys.exit(2)
    config_path = sys.argv[1]
    with open(config_path) as f:
        monitor_config = json.load(f)
    if "monitor_id" not in monitor_config:
        print("Config file is missing required 'monitor_id' field.", file=sys.stderr)
        sys.exit(2)
    # Never returns under normal operation -- this process is stopped
    # externally via SIGTERM (see monitor_manager.stop_monitor(monitor_id)).
    mm.run_monitor_loop(monitor_config)


if __name__ == "__main__":
    main()
