"""
atomic_io.py

Crash-safe JSON read/write helpers shared by every module in this app
that persists small pieces of state to disk (project_info.json,
users.json, roles.json, hpc_connections.json, downstream_recipe.json,
setup_install_status.json, ...).

--- The problem this fixes ---

Every one of those files is currently written like this:

    with open(path, "w") as f:
        json.dump(obj, f, indent=2)

Opening a file with mode "w" TRUNCATES IT TO ZERO BYTES IMMEDIATELY,
before a single byte of new content is written. There is a window --
short, but real -- during which the file on disk is empty or contains
only a partial JSON document.

If the process dies inside that window, the file stays that way
permanently. Causes are not exotic:

  - SIGHUP when an SSH session drops (this app's own deployment notes
    already document Streamlit and its children being killed this way)
  - the OOM killer on a shared HPC node, which targets exactly the kind
    of large-memory process this app spawns
  - `docker stop` / container eviction
  - the user pressing Ctrl-C in the terminal running Streamlit

The damage is then compounded by the read side:

    with open(path) as f:
        return json.load(f)

json.load() raises JSONDecodeError on a truncated file, and nothing
catches it. So a single interrupted write turns into an unhandled
traceback on EVERY page that touches that file -- for
project_info.json, that means the project can never be opened again;
for users.json, it means nobody can log in at all.

--- Why write-temp-then-rename fixes it ---

os.replace() is atomic on POSIX and on Windows: the destination path
either still refers to the OLD file, or it refers to the COMPLETE new
one. There is no intermediate state in which a reader can observe a
half-written file. A crash at any point leaves either the previous
good content or an orphaned temp file -- never a corrupt destination.

Two details matter and are easy to get wrong:

  1. The temp file MUST be created in the SAME DIRECTORY as the
     destination. os.replace() is only atomic within a single
     filesystem; a temp file in /tmp crossing onto an HPC scratch
     mount degrades to a non-atomic copy, reintroducing the bug.

  2. os.fsync() before the rename. Without it the rename can reach the
     disk before the data does, so a machine-level crash (not just a
     process-level one) can leave a correctly-named but empty file.

This is the same temp-then-atomic-rename discipline
reference_manager.ensure_shared_resource() already applies to
multi-GB reference builds, and whitelist_manager.download_whitelist()
applies to downloads. It simply never got applied to the small JSON
files that track whether all of that work succeeded.
"""
import errno
import json
import os
import tempfile


class CorruptStateFile(Exception):
    """
    Raised by read_json() when a state file exists but does not contain
    valid JSON, and the caller asked for on_corrupt="raise".

    Carries the path of the salvaged copy (see read_json) so a UI layer
    can tell the user exactly which file to inspect or restore.
    """

    def __init__(self, path, backup_path, original_error):
        self.path = path
        self.backup_path = backup_path
        self.original_error = original_error
        super().__init__(
            f"{path} exists but is not valid JSON ({original_error}). "
            f"The damaged file has been preserved at {backup_path} for inspection. "
            "This usually means the app was killed partway through writing it."
        )


def atomic_write_json(path, obj, mode=None, indent=2):
    """
    Write obj to path as JSON, atomically.

    path : destination file. Its parent directory is created if needed.
    obj  : any JSON-serializable object.
    mode : optional octal permission bits applied to the file BEFORE it
           is moved into place (e.g. 0o600 for users.json). Setting
           permissions on the temp file rather than after the rename
           means the file is never briefly world-readable at its real
           path.
    """
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)

    # Same directory as the destination -- required for os.replace() to
    # actually be atomic (see module docstring).
    fd, tmp_path = tempfile.mkstemp(
        dir=directory, prefix=f".{os.path.basename(path)}.", suffix=".tmp",
    )
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, indent=indent)
            f.flush()
            os.fsync(f.fileno())   # data hits disk before the rename
        if mode is not None:
            os.chmod(tmp_path, mode)
        os.replace(tmp_path, path)  # atomic
    except Exception:
        # Never leave a stray temp file behind on failure.
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


def read_json(path, default=None, on_corrupt="raise"):
    """
    Read JSON from path with an explicit, deliberate policy for what to
    do when the file exists but is damaged.

    default    : returned when the file simply DOES NOT EXIST. This is
                 the ordinary "nothing saved yet" case and is never an
                 error.

    on_corrupt : what to do when the file EXISTS but does not parse.
                 "raise"   -> preserve a copy and raise CorruptStateFile
                 "default" -> preserve a copy and return `default`

    --- Why on_corrupt defaults to "raise", and why that matters most
        for auth ---

    It is tempting to have every reader quietly fall back to `default`
    on a parse error, because it keeps the app running. For most state
    that is merely lossy. For authentication it is a SECURITY HOLE:

        auth_manager.any_users_exist() is len(_load_users()) > 0, and
        render_login_gate() shows the unauthenticated "create the first
        admin account" bootstrap form whenever that is False.

    So a corrupt users.json that silently reads as {} does not fail
    closed -- it fails OPEN, offering admin creation to whoever loads
    the page. "Missing" and "damaged" must not be collapsed into the
    same outcome. Callers must opt in to lenient behavior, per file,
    with their eyes open.

    In BOTH modes the damaged bytes are copied aside to
    <path>.corrupt-<n> before anything else happens, so a recoverable
    file is never destroyed by a subsequent write.
    """
    if not os.path.exists(path):
        return default

    try:
        with open(path) as f:
            return json.load(f)
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        backup_path = _preserve_corrupt_file(path)
        if on_corrupt == "default":
            return default
        raise CorruptStateFile(path, backup_path, e) from e


def _preserve_corrupt_file(path):
    """
    Copy a damaged state file aside to <path>.corrupt-<n>, choosing the
    first index not already taken so repeated failures never overwrite
    the first (usually most useful) captured copy.

    Uses os.replace() to MOVE rather than copy: the damaged file is
    already unreadable, and moving it means the next atomic_write_json()
    starts from a clean slate instead of layering on top of garbage.
    """
    index = 0
    while True:
        candidate = f"{path}.corrupt-{index}"
        if not os.path.exists(candidate):
            break
        index += 1
    try:
        os.replace(path, candidate)
    except OSError:
        return None
    return candidate


class ExclusiveLock:
    """
    A cross-process mutex built on a lock FILE, for guarding operations
    that must never run twice concurrently.

    --- The problem this fixes ---

    deployment_manager.launch_install() guards itself like this:

        if is_install_in_progress():           # reads a JSON file
            return False, "already in progress"
        ...
        _write_status({"status": "running"})   # writes it

    The check and the write are separate operations with a gap between
    them. Two users on the same deployment clicking "Install Selected
    Packages" at the same moment both read "not running", both proceed,
    and both then run `mamba install` against the SAME live conda
    environment while truncating each other's shared log file. That
    module's own docstring notes these installs modify the environment
    the app itself is running in, so two concurrent solves is a real
    corruption risk rather than a cosmetic glitch.

    fcntl.flock() closes the gap because acquiring the lock is a single
    atomic kernel operation -- there is no window between testing and
    taking it.

    --- Why this also fixes the permanently-wedged install ---

    The JSON status file records "running" but is only ever cleared by
    the same process that set it. If that process is killed, the file
    says "running" forever, is_install_in_progress() returns True
    forever, and no further install can be started through the UI ever
    again. (Both launch functions write "pid": None and never populate
    it, so there is no liveness check available either.)

    An OS-level lock does not have this failure mode: the kernel
    releases it automatically when the holding process exits, however
    it exits. A stale lock is therefore impossible by construction.

    Usage:

        lock = ExclusiveLock(data_path("locks", "install.lock"))
        if not lock.acquire():
            return False, "An install is already running."
        try:
            ...
        finally:
            lock.release()
    """

    def __init__(self, lock_path):
        self.lock_path = lock_path
        self._fd = None

    def acquire(self):
        """
        Try to take the lock without blocking. Returns True on success,
        False if another live process already holds it.
        """
        import fcntl

        os.makedirs(os.path.dirname(os.path.abspath(self.lock_path)), exist_ok=True)
        fd = os.open(self.lock_path, os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as e:
            os.close(fd)
            if e.errno in (errno.EACCES, errno.EAGAIN):
                return False   # someone else holds it -- expected, not an error
            raise
        # Record who holds it, purely for diagnostics in the UI.
        try:
            os.ftruncate(fd, 0)
            os.write(fd, str(os.getpid()).encode())
            os.fsync(fd)
        except OSError:
            pass
        self._fd = fd
        return True

    def release(self):
        if self._fd is None:
            return
        import fcntl
        try:
            fcntl.flock(self._fd, fcntl.LOCK_UN)
        finally:
            os.close(self._fd)
            self._fd = None

    def holder_pid(self):
        """
        Best-effort read of the PID recorded in the lock file, for
        showing "an install started by process N is already running"
        rather than an anonymous refusal. Returns None if unreadable.
        """
        try:
            with open(self.lock_path) as f:
                return int(f.read().strip())
        except (OSError, ValueError):
            return None

    def __enter__(self):
        if not self.acquire():
            raise RuntimeError(f"Could not acquire lock: {self.lock_path}")
        return self

    def __exit__(self, *exc):
        self.release()
        return False
