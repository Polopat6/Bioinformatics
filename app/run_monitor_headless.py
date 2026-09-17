#!/usr/bin/env python3
"""
run_monitor_headless.py

Entry point for the Monitor Mode watcher daemon -- mirrors
run_pipeline_headless.py's role for individual pipeline runs, but for
the long-lived directory-watching loop instead. Launched (detached)
by monitor_manager.launch_monitor_daemon(), and left running
indefinitely until explicitly stopped (see monitor_manager.stop_monitor,
which sends SIGTERM to this process's PID).

Usage:
    python run_monitor_headless.py <monitor_config_json_path>
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from monitor_manager import MonitorConfig, run_monitor_loop


def main():
    if len(sys.argv) != 2:
        print("Usage: python run_monitor_headless.py <monitor_config_json_path>", file=sys.stderr)
        sys.exit(2)

    import json
    try:
        with open(sys.argv[1]) as f:
            monitor_config = MonitorConfig.from_dict(json.load(f))
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as e:
        # This process is launched DETACHED, so an uncaught traceback
        # here goes nowhere a user will ever see -- the monitor just
        # silently never starts. Fail with an explicit message and a
        # distinct exit code instead.
        print(f"Could not load monitor config from {sys.argv[1]!r}: {e}", file=sys.stderr)
        sys.exit(3)

    run_monitor_loop(monitor_config)  # runs forever until killed


if __name__ == "__main__":
    main()
