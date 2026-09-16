"""
auth_manager.py

Lightweight, dependency-free username/password authentication with
ADMIN-CUSTOMIZABLE role-based permissions -- built for a single
researcher or a small lab's own private deployment of this portal, NOT
for public-facing multi-tenant use (see "What this does and does NOT
protect against" below for the honest scope).

--- Why stdlib-only, not a package (2026-08-24) ---
At the scale this is actually needed for (a handful of named users in
one lab), a full third-party auth package is more machinery than the
problem requires. hashlib's pbkdf2_hmac (a real, slow, salted key-
derivation function) plus the `secrets` module for cryptographically-
secure randomness cover everything needed, with zero new entries added
to environment.yml.

--- What this does and does NOT protect against ---
This module gates ACTIONS INSIDE THE APP by checking the CURRENT
session's role's granted permissions. It does NOT, and cannot, protect
the underlying data FILES themselves from anyone who already has SSH/
shell access to the same HPC account, `docker exec` access into a
running container, or direct filesystem/root access. A login screen
stops someone who only has the app's URL; it does nothing for someone
who already has a shell. This module is the right tool for "my
lab-mate with a limited role shouldn't be able to delete my finished
runs" -- it is NOT the right tool for "protect this data from anyone
who can already log into the server."

--- Custom roles + permissions (2026-08-24 redesign) ---
Earlier revision of this module used a fixed, hardcoded two-role model
(admin/tech). Replaced with a genuinely admin-customizable system:

  - PERMISSION_CATALOG defines every individually grantable capability
    this app currently exposes (e.g. "delete_projects",
    "manage_eggnog_database") -- each with a plain-language label and
    description for the role-editor UI. Adding a new gated action
    anywhere in the app in the future means adding ONE new entry here
    plus the corresponding has_permission()/require_permission() check
    at that action's own call site -- nothing else needs to change.

  - Exactly ONE built-in, IMMUTABLE role always exists: ROLE_ADMIN
    ("admin"). It implicitly holds EVERY permission in
    PERMISSION_CATALOG, including any added to the catalog in the
    FUTURE -- implemented as a wildcard check in has_permission(),
    not as an enumerated list that could silently miss a newly-added
    permission. This role's own permission set cannot be edited or
    reduced, and creating/deleting USER ACCOUNTS and ROLES themselves
    (see "Why user/role management stays admin-only" below) is ONLY
    ever available to this exact role -- never delegable to any
    custom role, no matter what permissions that role is granted.

  - Exactly ONE built-in, but EDITABLE, default role also always
    exists: ROLE_TECH ("tech") -- the role new accounts default to.
    Its permission set starts as this app's own previous behavior
    (day-to-day pipeline use, no destructive/admin-flagged actions),
    but an admin can freely ADD permissions to it, same as any custom
    role. It cannot be deleted outright (there must always be at least
    one non-admin default role for ordinary accounts to fall back to),
    but its permission set is fully admin-editable.

  - Any number of ADDITIONAL custom roles can be created by an admin,
    each an arbitrary named subset of PERMISSION_CATALOG's keys (e.g.
    a "senior_tech" role with delete_projects but not
    manage_eggnog_database). Custom roles can be freely edited or
    deleted (deletion is blocked while any user account still holds
    that role, to avoid silently leaving a user with an undefined
    role -- the caller must reassign those users first).

--- Why user/role management stays admin-only, never delegable
    (2026-08-24) ---
If a custom role could be granted permission to create new users or
define new roles, a holder of that role could simply create a brand
new user assigned to the built-in "admin" role (or invent a role that
grants itself every other permission) -- a straightforward privilege-
escalation path. To close this off structurally rather than relying on
careful admin configuration to avoid it, "manage users" and "manage
roles" are NOT permissions in PERMISSION_CATALOG at all -- they are
hardcoded to require the exact built-in ROLE_ADMIN role (checked via
is_admin(), not has_permission()) everywhere they're used, with no
custom-role path to reach them under any configuration.

--- Storage format ---
Users: data/auth/users.json -- {username: {"password_hash", "salt",
"iterations", "role", "created_at", "created_by"}}. Plaintext
passwords are NEVER stored or logged anywhere.
Roles: data/auth/roles.json -- {role_name: {"permissions": [...],
"is_builtin": bool, "created_at", "created_by"}}. The two built-in
roles (admin, tech) are seeded automatically the first time this file
would otherwise not exist, or the first time either is found missing
from an existing file (e.g. upgrading from the prior fixed-role
version of this module) -- see _ensure_builtin_roles().

--- First-run bootstrap ---
A brand-new deployment has no users at all -- any_users_exist() lets
the calling UI show a one-time "create the first admin account" form.
Every subsequent account must be created by an existing admin.
"""
import hashlib
import json
import os
import secrets
from datetime import datetime

import streamlit as st
import app_paths
import atomic_io
# ---------------------------------------------------------------------------
# Storage locations + built-in role constants
# ---------------------------------------------------------------------------

AUTH_DIR = app_paths.data_path("auth")
USERS_PATH = os.path.join(AUTH_DIR, "users.json")
ROLES_PATH = os.path.join(AUTH_DIR, "roles.json")

ROLE_ADMIN = "admin"
ROLE_TECH = "tech"

# The DEFAULT permission set the built-in "tech" role is seeded with --
# deliberately matching exactly what an ordinary (non-admin) user could
# already do BEFORE any of this app's admin-gated actions existed, so
# introducing this permission system is not itself a silent behavior
# regression for existing day-to-day use. An admin is free to add to
# (or remove from) this set at any time after the fact.
_DEFAULT_TECH_PERMISSIONS = ["install_dependencies", "manage_hpc_connections"]

# Every individually-grantable capability this app currently exposes.
# Adding a new gated action elsewhere in the app = add ONE entry here
# + the corresponding has_permission()/require_permission() check at
# that action's own call site. Deliberately does NOT include anything
# related to managing users or roles themselves -- see this module's
# own docstring, "Why user/role management stays admin-only", for why
# that is handled entirely separately (via is_admin()) and can never
# be granted to a custom role under any configuration.
PERMISSION_CATALOG = {
    "delete_projects": {
        "label": "Delete Projects",
        "description": "Permanently delete Bulk RNA-Seq or Single-cell RNA-Seq projects and all their files.",
    },
    "manage_eggnog_database": {
        "label": "Manage eggNOG Database",
        "description": "Download or re-download the large (~49GB), shared eggNOG-mapper orthology database.",
    },
    "install_dependencies": {
        "label": "Install Dependencies",
        "description": "Install missing Python/CLI/R packages into the running environment from the Setup & Deployment page.",
    },
    "manage_hpc_connections": {
        "label": "Manage HPC Connections",
        "description": "Add, test, or delete saved SSH connection profiles to remote HPC clusters.",
    },
    "delete_archives": {
        "label": "Delete Archives",
        "description": "Permanently delete archived projects or monitors (planned feature -- reserved here so it can be gated the moment it ships, without a separate later migration).",
    },
}

# PBKDF2 iteration count -- deliberately high (current, e.g. OWASP-
# ballpark, guidance for PBKDF2-HMAC-SHA256). Stored per-user (see
# verify_login()) so this constant can be safely raised later without
# invalidating existing accounts.
_PBKDF2_ITERATIONS = 600_000
_SALT_BYTES = 16


# ---------------------------------------------------------------------------
# Low-level password hashing
# ---------------------------------------------------------------------------

def _hash_password(password, salt_bytes):
    return hashlib.pbkdf2_hmac(
        "sha256", password.encode("utf-8"), salt_bytes, _PBKDF2_ITERATIONS,
    ).hex()


def _load_users():
    # on_corrupt="raise" is deliberate and load-bearing: a damaged
    # users.json must NEVER read as {} here, because any_users_exist()
    # would then be False and render_login_gate() would show the
    # unauthenticated "create first admin" bootstrap form.
    return atomic_io.read_json(USERS_PATH, default={}, on_corrupt="raise")


def _save_users(users):
    atomic_io.atomic_write_json(USERS_PATH, users, mode=0o600)


# ---------------------------------------------------------------------------
# Role storage -- built-in roles auto-seeded, custom roles admin-managed
# ---------------------------------------------------------------------------

def _load_roles():
    return atomic_io.read_json(ROLES_PATH, default={}, on_corrupt="raise")


def _save_roles(roles):
    atomic_io.atomic_write_json(ROLES_PATH, roles, mode=0o600)


def _ensure_builtin_roles():
    """
    Guarantee both built-in roles (admin, tech) exist in roles.json,
    seeding either one that's missing -- called defensively at the top
    of every role-reading/writing function below, so this module works
    correctly both on a brand-new deployment (roles.json doesn't exist
    yet at all) AND when upgrading an existing deployment that predates
    this custom-role system (roles.json may not exist yet even though
    users.json already does, in which case both built-ins are seeded
    fresh here on first access).
    """
    roles = _load_roles()
    changed = False
    if ROLE_ADMIN not in roles:
        roles[ROLE_ADMIN] = {
            "permissions": "__all__",  # wildcard sentinel -- see has_permission()
            "is_builtin": True,
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "created_by": None,
        }
        changed = True
    if ROLE_TECH not in roles:
        roles[ROLE_TECH] = {
            "permissions": list(_DEFAULT_TECH_PERMISSIONS),
            "is_builtin": True,
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "created_by": None,
        }
        changed = True
    if changed:
        _save_roles(roles)
    return roles


# ---------------------------------------------------------------------------
# Public role-management API (admin-only at the CALL SITE -- see
# render_role_management(), which is the only UI entry point and itself
# requires is_admin() -- these functions do not re-check it themselves,
# exactly mirroring create_user()/change_role()/delete_user()'s own
# existing convention of being plain backend logic with authorization
# enforced by the calling UI, not duplicated in every function)
# ---------------------------------------------------------------------------

def list_roles():
    """
    Return a list of {"role_name", "permissions", "is_builtin",
    "n_users"} dicts for every existing role (both built-in and
    custom) -- for the role-management UI and for populating the role-
    picker dropdown when creating/editing a user. "permissions" is
    either a list of PERMISSION_CATALOG keys, or the literal string
    "__all__" for the built-in admin role (never expanded into an
    enumerated list, so a future PERMISSION_CATALOG addition is
    automatically covered without needing to update every existing
    admin role record).
    """
    roles = _ensure_builtin_roles()
    users = _load_users()
    role_user_counts = {}
    for record in users.values():
        role_user_counts[record["role"]] = role_user_counts.get(record["role"], 0) + 1

    return [
        {
            "role_name": role_name,
            "permissions": record["permissions"],
            "is_builtin": record.get("is_builtin", False),
            "n_users": role_user_counts.get(role_name, 0),
        }
        for role_name, record in sorted(roles.items())
    ]


def role_exists(role_name):
    roles = _ensure_builtin_roles()
    return role_name in roles


def get_role_permissions(role_name):
    "Return the given role's permission list, or the literal string '__all__' for the built-in admin role. Returns [] if role_name doesn't exist at all (a defensive default, not an error -- a caller checking a specific permission against a since-deleted role should simply see 'no permissions', not crash)."
    roles = _ensure_builtin_roles()
    record = roles.get(role_name)
    if record is None:
        return []
    return record["permissions"]


def create_role(role_name, permissions):
    """
    Create a new CUSTOM role (never used for the two built-in roles,
    which are seeded automatically by _ensure_builtin_roles() and can
    never be created/renamed this way).

    role_name: must not already exist (built-in or custom), and must
        not be empty.
    permissions: list of PERMISSION_CATALOG keys this role grants --
        any key not in PERMISSION_CATALOG is silently ignored (not
        stored), so a stale/typo'd permission key can never linger in
        a role's stored permission list.

    Returns (success: bool, message: str).
    """
    roles = _ensure_builtin_roles()
    role_key = role_name.strip()
    if not role_key:
        return False, "Role name cannot be empty."
    if role_key in roles:
        return False, f"A role named '{role_key}' already exists."

    valid_permissions = [p for p in permissions if p in PERMISSION_CATALOG]
    roles[role_key] = {
        "permissions": valid_permissions,
        "is_builtin": False,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "created_by": current_username(),
    }
    _save_roles(roles)
    return True, f"Role '{role_key}' created with {len(valid_permissions)} permission(s)."


def update_role_permissions(role_name, permissions):
    """
    Replace an existing role's permission set. Refuses to modify the
    built-in "admin" role at all (its permissions are the fixed "__all__"
    wildcard, not a real editable list) -- the built-in "tech" role and
    any custom role ARE freely editable this way.

    Returns (success: bool, message: str).
    """
    roles = _ensure_builtin_roles()
    if role_name not in roles:
        return False, f"No such role: '{role_name}'."
    if role_name == ROLE_ADMIN:
        return False, "The built-in 'admin' role always has every permission and cannot be edited."

    valid_permissions = [p for p in permissions if p in PERMISSION_CATALOG]
    roles[role_name]["permissions"] = valid_permissions
    _save_roles(roles)
    return True, f"Role '{role_name}' updated with {len(valid_permissions)} permission(s)."


def delete_role(role_name):
    """
    Permanently delete a CUSTOM role. Refuses to delete either built-in
    role (admin or tech -- see this module's own docstring for why
    "tech" specifically must always exist as a fallback default).
    Refuses to delete a role that any user account currently holds --
    the caller (UI) must reassign those users to a different role
    first, so no user is ever silently left with an undefined role.

    Returns (success: bool, message: str).
    """
    roles = _ensure_builtin_roles()
    if role_name not in roles:
        return False, f"No such role: '{role_name}'."
    if roles[role_name].get("is_builtin"):
        return False, f"'{role_name}' is a built-in role and cannot be deleted."

    users = _load_users()
    holders = [u for u, record in users.items() if record["role"] == role_name]
    if holders:
        return False, (
            f"Cannot delete '{role_name}' -- {len(holders)} user account(s) currently hold it "
            f"({', '.join(holders)}). Reassign them to a different role first."
        )

    del roles[role_name]
    _save_roles(roles)
    return True, f"Role '{role_name}' deleted."


# ---------------------------------------------------------------------------
# Public user-management API
# ---------------------------------------------------------------------------

def any_users_exist():
    return len(_load_users()) > 0


def username_exists(username):
    return username.strip().lower() in _load_users()


def create_user(username, password, role, created_by=None):
    """
    Create a new user account.

    role: must be an EXISTING role name (built-in "admin"/"tech", or
        any admin-created custom role) -- validated via role_exists(),
        NOT against a fixed tuple, since the whole point of this
        module's custom-role redesign is that the valid set of roles
        is admin-defined and open-ended, not a hardcoded pair.

    Returns (success: bool, message: str).
    """
    username_key = username.strip().lower()
    if not username_key:
        return False, "Username cannot be empty."
    if not password:
        return False, "Password cannot be empty."
    if not role_exists(role):
        return False, f"Unknown role {role!r} -- it must already exist (see Role Management)."

    users = _load_users()
    if username_key in users:
        return False, f"A user named '{username_key}' already exists."

    salt_bytes = secrets.token_bytes(_SALT_BYTES)
    users[username_key] = {
        "password_hash": _hash_password(password, salt_bytes),
        "salt": salt_bytes.hex(),
        "iterations": _PBKDF2_ITERATIONS,
        "role": role,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "created_by": created_by,
    }
    _save_users(users)
    return True, f"User '{username_key}' created with role '{role}'."


def verify_login(username, password):
    """
    Check a username/password pair against stored credentials. Returns
    (success: bool, role_or_None: str, message: str). Uses
    secrets.compare_digest() for a constant-time comparison, and
    returns the SAME generic message for both "no such user" and
    "wrong password" (no username-enumeration leak).
    """
    username_key = username.strip().lower()
    users = _load_users()
    record = users.get(username_key)
    if record is None:
        return False, None, "Incorrect username or password."

    salt_bytes = bytes.fromhex(record["salt"])
    iterations = record.get("iterations", _PBKDF2_ITERATIONS)
    computed_hash = hashlib.pbkdf2_hmac(
        "sha256", password.encode("utf-8"), salt_bytes, iterations,
    ).hex()

    if not secrets.compare_digest(computed_hash, record["password_hash"]):
        return False, None, "Incorrect username or password."

    return True, record["role"], f"Welcome, {username_key}."


def list_users():
    "Return a list of {'username', 'role', 'created_at', 'created_by'} dicts -- never includes password_hash/salt."
    users = _load_users()
    return [
        {
            "username": username,
            "role": record["role"],
            "created_at": record.get("created_at"),
            "created_by": record.get("created_by"),
        }
        for username, record in sorted(users.items())
    ]


def change_role(username, new_role):
    "Change an existing user's role to any EXISTING role (built-in or custom). Returns (success, message)."
    if not role_exists(new_role):
        return False, f"Unknown role {new_role!r} -- it must already exist (see Role Management)."
    username_key = username.strip().lower()
    users = _load_users()
    if username_key not in users:
        return False, f"No such user: '{username_key}'."
    users[username_key]["role"] = new_role
    _save_users(users)
    return True, f"'{username_key}' is now role '{new_role}'."


def delete_user(username):
    username_key = username.strip().lower()
    users = _load_users()
    if username_key not in users:
        return False, f"No such user: '{username_key}'."
    del users[username_key]
    _save_users(users)
    return True, f"User '{username_key}' deleted."


def change_password(username, new_password):
    if not new_password:
        return False, "Password cannot be empty."
    username_key = username.strip().lower()
    users = _load_users()
    if username_key not in users:
        return False, f"No such user: '{username_key}'."
    salt_bytes = secrets.token_bytes(_SALT_BYTES)
    users[username_key]["password_hash"] = _hash_password(new_password, salt_bytes)
    users[username_key]["salt"] = salt_bytes.hex()
    users[username_key]["iterations"] = _PBKDF2_ITERATIONS
    _save_users(users)
    return True, "Password updated."


# ---------------------------------------------------------------------------
# Session-state integration
# ---------------------------------------------------------------------------

_SESSION_USERNAME_KEY = "_auth_username"
_SESSION_ROLE_KEY = "_auth_role"


def is_logged_in():
    return _SESSION_USERNAME_KEY in st.session_state


def current_username():
    return st.session_state.get(_SESSION_USERNAME_KEY)


def current_role():
    return st.session_state.get(_SESSION_ROLE_KEY)


def is_admin():
    "True only if the CURRENT session holds the exact built-in ROLE_ADMIN role -- this is the ONLY check used for user/role management (see this module's own docstring, 'Why user/role management stays admin-only') and is never satisfiable by any custom role, however permissive."
    return current_role() == ROLE_ADMIN


def has_permission(permission_key):
    """
    The primary, forward-looking authorization check for every
    OPERATIONAL (non-user/role-management) gated action in this app --
    e.g. has_permission("delete_projects"). True if the current
    session's role is the built-in admin role (wildcard -- see
    PERMISSION_CATALOG's own docstring) OR if that role's own stored
    permission list explicitly includes permission_key.

    Returns False (never raises) if nobody is logged in, or if
    permission_key isn't a real PERMISSION_CATALOG key at all -- a
    typo'd permission key at a call site should fail closed (deny),
    not silently grant access or crash.
    """
    role = current_role()
    if role is None:
        return False
    if role == ROLE_ADMIN:
        return True
    if permission_key not in PERMISSION_CATALOG:
        return False
    return permission_key in get_role_permissions(role)


def logout():
    st.session_state.pop(_SESSION_USERNAME_KEY, None)
    st.session_state.pop(_SESSION_ROLE_KEY, None)


def require_permission(permission_key, message=None):
    """
    Guard clause for gating a specific ACTION behind a permission
    check -- the primary mechanism going forward (see has_permission()
    above). Renders an st.error() and returns False if the check
    fails; does NOT call st.stop() itself, so the caller can decide
    whether to keep rendering the rest of its page.
    """
    if has_permission(permission_key):
        return True
    label = PERMISSION_CATALOG.get(permission_key, {}).get("label", permission_key)
    st.error(message or f"⚠️ This action requires the '{label}' permission.")
    return False


def require_admin(message=None):
    "Guard clause for the small, fixed set of actions that ONLY the exact built-in admin role can ever perform (creating/deleting users, creating/editing/deleting roles) -- see this module's own docstring for why these are never permission-catalog-based or delegable to any custom role."
    if is_admin():
        return True
    st.error(message or "⚠️ This action requires administrator access.")
    return False


# ---------------------------------------------------------------------------
# Login / bootstrap UI
# ---------------------------------------------------------------------------

def render_login_gate():
    if is_logged_in():
        return True
    try:
        users_exist = any_users_exist()
    except atomic_io.CorruptStateFile as e:
        st.error(
            "🔒 The user account file is damaged, so login is disabled. "
            "This is a fail-closed safety measure -- it does NOT mean your "
            "accounts are gone."
        )
        st.caption(f"Damaged file preserved at: {e.backup_path}")
        return False
    if not users_exist:
        _render_bootstrap_admin_form()
        return False
    _render_login_form()
    return False


def _render_bootstrap_admin_form():
    st.title("🔐 Welcome -- First-Time Setup")
    st.markdown(
        "No user accounts exist yet for this deployment. Create the "
        "**first administrator account** below -- only this account "
        "(and any other 'admin' account it creates) can manage users, "
        "define custom roles, and grant permissions to them afterward."
    )
    with st.form("auth_bootstrap_admin_form"):
        username = st.text_input("Admin username:")
        password = st.text_input("Admin password:", type="password")
        password_confirm = st.text_input("Confirm password:", type="password")
        submitted = st.form_submit_button("Create Admin Account", type="primary")

    if submitted:
        if password != password_confirm:
            st.error("Passwords do not match.")
        elif len(password) < 8:
            st.error("Password must be at least 8 characters.")
        else:
            _ensure_builtin_roles()
            success, message = create_user(username, password, ROLE_ADMIN, created_by=None)
            if success:
                st.success(f"{message} You can now log in below.")
                st.rerun()
            else:
                st.error(message)


def _render_login_form():
    st.title("🔐 Log In")
    with st.form("auth_login_form"):
        username = st.text_input("Username:")
        password = st.text_input("Password:", type="password")
        submitted = st.form_submit_button("Log In", type="primary")

    if submitted:
        success, role, message = verify_login(username, password)
        if success:
            st.session_state[_SESSION_USERNAME_KEY] = username.strip().lower()
            st.session_state[_SESSION_ROLE_KEY] = role
            st.rerun()
        else:
            st.error(message)


def render_user_badge():
    username = current_username()
    role = current_role()
    if not username:
        return
    st.sidebar.markdown(f"👤 **{username}** ({role})")
    if st.sidebar.button("🚪 Log Out", key="auth_logout_btn", use_container_width=True):
        logout()
        st.rerun()


# ---------------------------------------------------------------------------
# Admin-only user management page
# ---------------------------------------------------------------------------

def render_user_management():
    """
    Full user-management page: create/delete accounts, change roles,
    reset passwords, PLUS role/permission management (create custom
    roles, edit "tech"/custom role permissions). Both halves require
    the exact built-in admin role (require_admin()) -- see this
    module's own docstring for why role/user management is never
    delegable to any custom role, however permissive.
    """
    st.title("👥 User & Role Management")
    if not require_admin("⚠️ Only administrators can manage users and roles."):
        return

    _render_users_section()
    st.markdown("---")
    _render_roles_section()


def _render_users_section():
    st.subheader("Users")
    st.markdown(
        "Create and manage accounts for this lab's users. Each account "
        "is assigned exactly one role, which determines what it can do -- "
        "see **Roles & Permissions** below to see or customize what each "
        "role actually grants."
    )

    users = list_users()
    role_names = [r["role_name"] for r in list_roles()]

    if not users:
        st.info("No users found -- this shouldn't be possible if you're viewing this page.")
    else:
        for user in users:
            col1, col2, col3, col4 = st.columns([2, 2, 2, 1])
            with col1:
                st.markdown(f"**{user['username']}**" + (" (you)" if user["username"] == current_username() else ""))
            with col2:
                new_role = st.selectbox(
                    "Role:", options=role_names,
                    index=role_names.index(user["role"]) if user["role"] in role_names else 0,
                    key=f"auth_role_select_{user['username']}",
                    label_visibility="collapsed",
                )
                if new_role != user["role"]:
                    if st.button("Update", key=f"auth_role_update_{user['username']}"):
                        success, message = change_role(user["username"], new_role)
                        if success:
                            st.success(message)
                            st.rerun()
                        else:
                            st.error(message)
            with col3:
                st.caption(f"Created: {user.get('created_at', '—')}")
            with col4:
                is_last_admin = (
                    user["role"] == ROLE_ADMIN
                    and sum(1 for u in users if u["role"] == ROLE_ADMIN) == 1
                )
                if st.button(
                    "🗑️", key=f"auth_delete_user_{user['username']}",
                    disabled=is_last_admin,
                    help="Cannot delete the last remaining admin account." if is_last_admin else "Delete this user",
                ):
                    st.session_state[f"auth_confirm_delete_{user['username']}"] = True
                    st.rerun()

            if st.session_state.get(f"auth_confirm_delete_{user['username']}"):
                st.warning(f"Delete user '{user['username']}'? This cannot be undone.")
                cc1, cc2 = st.columns(2)
                with cc1:
                    if st.button("Yes, delete", key=f"auth_confirm_delete_yes_{user['username']}"):
                        success, message = delete_user(user["username"])
                        st.session_state.pop(f"auth_confirm_delete_{user['username']}", None)
                        if success:
                            st.success(message)
                        else:
                            st.error(message)
                        st.rerun()
                with cc2:
                    if st.button("Cancel", key=f"auth_confirm_delete_cancel_{user['username']}"):
                        st.session_state.pop(f"auth_confirm_delete_{user['username']}", None)
                        st.rerun()

    st.markdown("**➕ Add New User**")
    with st.form("auth_add_user_form"):
        new_username = st.text_input("Username:")
        new_password = st.text_input("Password:", type="password")
        new_role = st.selectbox("Role:", options=role_names)
        submitted = st.form_submit_button("Create User", type="primary")

    if submitted:
        if len(new_password) < 8:
            st.error("Password must be at least 8 characters.")
        else:
            success, message = create_user(new_username, new_password, new_role, created_by=current_username())
            if success:
                st.success(message)
                st.rerun()
            else:
                st.error(message)

    st.markdown("---")
    st.markdown("**🔑 Change My Own Password**")
    with st.form("auth_change_own_password_form"):
        current_pw = st.text_input("Current password:", type="password")
        new_pw = st.text_input("New password:", type="password")
        new_pw_confirm = st.text_input("Confirm new password:", type="password")
        submitted_pw = st.form_submit_button("Change Password")

    if submitted_pw:
        verify_success, _, _ = verify_login(current_username(), current_pw)
        if not verify_success:
            st.error("Current password is incorrect.")
        elif new_pw != new_pw_confirm:
            st.error("New passwords do not match.")
        elif len(new_pw) < 8:
            st.error("New password must be at least 8 characters.")
        else:
            success, message = change_password(current_username(), new_pw)
            if success:
                st.success(message)
            else:
                st.error(message)


def _render_roles_section():
    st.subheader("Roles & Permissions")
    st.markdown(
        "The built-in **admin** role always has every permission and cannot be edited. "
        "The built-in **tech** role is the default for new users and starts with a "
        "reasonable baseline -- you can freely add or remove permissions from it below. "
        "You can also create any number of your own custom roles (e.g. a "
        "'senior_tech' role that can delete projects but not manage the eggNOG database)."
    )

    roles = list_roles()
    for role in roles:
        with st.expander(
            f"{'🔒 ' if role['is_builtin'] else ''}{role['role_name']}"
            + (" (built-in, all permissions)" if role["permissions"] == "__all__" else "")
            + f" — {role['n_users']} user(s)",
            expanded=False,
        ):
            if role["permissions"] == "__all__":
                st.caption("This role always has every permission this app defines -- nothing to configure.")
            else:
                current_perms = set(role["permissions"])
                new_perms = []
                for perm_key, perm_info in PERMISSION_CATALOG.items():
                    checked = st.checkbox(
                        perm_info["label"], value=(perm_key in current_perms),
                        key=f"auth_role_perm_{role['role_name']}_{perm_key}",
                        help=perm_info["description"],
                    )
                    if checked:
                        new_perms.append(perm_key)

                if set(new_perms) != current_perms:
                    if st.button("💾 Save Permission Changes", key=f"auth_role_save_{role['role_name']}"):
                        success, message = update_role_permissions(role["role_name"], new_perms)
                        if success:
                            st.success(message)
                            st.rerun()
                        else:
                            st.error(message)

            if not role["is_builtin"]:
                if role["n_users"] > 0:
                    st.caption(f"⚠️ Cannot delete -- {role['n_users']} user(s) currently hold this role.")
                else:
                    if st.button("🗑️ Delete This Role", key=f"auth_role_delete_{role['role_name']}"):
                        success, message = delete_role(role["role_name"])
                        if success:
                            st.success(message)
                            st.rerun()
                        else:
                            st.error(message)

    st.markdown("**➕ Create a New Custom Role**")
    with st.form("auth_create_role_form"):
        new_role_name = st.text_input("Role name (e.g. 'senior_tech'):")
        st.caption("Select which permissions this new role should grant:")
        selected_perms = []
        for perm_key, perm_info in PERMISSION_CATALOG.items():
            if st.checkbox(perm_info["label"], key=f"auth_new_role_perm_{perm_key}", help=perm_info["description"]):
                selected_perms.append(perm_key)
        submitted_role = st.form_submit_button("Create Role", type="primary")

    if submitted_role:
        success, message = create_role(new_role_name, selected_perms)
        if success:
            st.success(message)
            st.rerun()
        else:
            st.error(message)
