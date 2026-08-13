"""
file_browser.py

Reusable "browse files/directories already on this server" widgets --
an alternative to Streamlit's built-in st.file_uploader for situations
where the app is running on a remote host (an HPC node, a shared lab
server, a personal Linux/Mac workstation used as a server) that ALREADY
has the desired file(s) sitting on its own disk. In that situation,
st.file_uploader is actively counterproductive: it always transfers a
file FROM the browser's local machine, over HTTP, even when the file
already exists on the exact machine Streamlit itself is running on --
for a multi-GB genome FASTA, or a whole directory of raw FASTQ files,
this means needlessly downloading them to a laptop and then
re-uploading them, when the app could simply point at the existing
path(s) directly.

Three widgets are provided:
  - render_server_file_browser: navigate to and select a SINGLE file
    (e.g. a genome FASTA, GTF annotation, or metadata spreadsheet).
  - render_server_directory_browser: navigate to and confirm a
    DIRECTORY itself (e.g. a folder already containing a full set of
    FASTQ files for many samples) -- appropriate when the caller needs
    to work with a whole directory's contents rather than one specific
    file.
  - (internal) _render_directory_navigator: the shared navigation UI
    used by both of the above.

Design goals (shared by all widgets):
  - Reusable across every place in the app that currently only offers
    st.file_uploader for a large reference/input file or a batch of
    files (motivating cases: alignment_workspace.py's custom genome/GTF
    upload, bulk_rnaseq_workspace.py's FASTQ upload, and
    bulk_rnaseq_workspace.py's metadata file upload).
  - SANDBOXED: navigation is restricted to a single configured root
    directory (and its subdirectories) -- never outside it. This is a
    hard security boundary, not a UI convenience; see
    _is_path_within_root() for the actual enforcement, which resolves
    symlinks and ".." components before checking, rather than doing a
    naive string-prefix comparison (which a symlink or ".." could
    trivially bypass). The DEFAULT root (see _default_browse_root
    below) is the filesystem root "/" itself -- i.e. by default there
    is effectively no sandboxing restriction beyond the underlying OS's
    own file permissions, since real deployments (e.g. an HPC user
    needing to reach shared /scratch or /disk storage well outside
    their own home directory) were found in practice to need to
    navigate to arbitrary locations on the machine, not just their own
    home directory tree. Callers that DO want a narrower, safer
    restriction (e.g. a genuinely multi-tenant deployment where users
    shouldn't browse each other's files) should pass an explicit,
    narrower `root_dir` themselves.
  - FAST NAVIGATION: since the default root is now the entire
    filesystem, reaching a deeply nested real-world path (e.g. an HPC
    /scratch project directory several levels deep) purely by clicking
    one subdirectory at a time would be slow and tedious. Two features
    address this directly (see _render_directory_navigator):
      1. A "jump directly to a path" text input -- type a full known
         path and go there in one step.
      2. A live-filtered subdirectory list -- type a partial name and
         the subdirectory list narrows to matches as you type (the
         closest practical equivalent to shell tab-completion that
         Streamlit's widget model actually supports -- genuine
         browser-level tab-autocomplete isn't something Streamlit
         exposes an API for, so this live-filter approach is the
         closest available substitute for quickly narrowing a long
         subdirectory listing without scrolling).
  - Returns an existing file's or directory's ABSOLUTE PATH, not its
    contents -- nothing is read into memory or copied anywhere by any
    widget in this module. The whole point is to avoid moving/
    duplicating large files that are already sitting in a perfectly
    good location; the caller is responsible for deciding what to do
    with the returned path (e.g. using it directly, or -- for the
    FASTQ-directory case -- symlinking individual files from it into a
    project's own fastq_dir so existing downstream code that expects
    real files/symlinks at fixed project-relative paths keeps working
    unmodified).
"""

import os

import streamlit as st


def _default_browse_root():
    """
    The default starting root for the file/directory browser, when the
    caller doesn't specify one explicitly.

    Defaults to "/" (the filesystem root) -- i.e. no sandboxing
    restriction by default beyond the OS's own file permissions. This
    was changed from an earlier default of the user's home directory
    after real deployment testing on an HPC found that navigation
    needed to reach shared storage locations (e.g. /scratch or /disk
    mounts) well outside the user's own home directory tree, and a
    home-directory sandbox made those genuinely necessary locations
    unreachable. Callers that specifically want a narrower restriction
    (e.g. a genuinely multi-tenant deployment where per-user isolation
    matters) should pass an explicit, narrower `root_dir` of their own
    rather than relying on this default.
    """
    return "/"


def _is_path_within_root(candidate_path, root_path):
    """
    Security check: verify that candidate_path is genuinely located
    within root_path (or is root_path itself), resolving symlinks and
    any ".." / "." components on BOTH sides first via os.path.realpath.

    This must use realpath (not just os.path.abspath or a plain string
    prefix check) for two reasons:
      1. A path containing ".." (e.g. "root/../../../etc/passwd") would
         pass a naive string-prefix check against "root" if compared
         before normalization -- realpath fully resolves these.
      2. A symlink INSIDE root_path that happens to point somewhere
         else entirely (e.g. "root/escape_hatch -> /etc") would also
         defeat a plain prefix check on the un-resolved path string --
         realpath follows symlinks to their real, final target before
         we compare.

    Returns True only if the resolved candidate_path is root_path itself
    or a genuine descendant of it.
    """
    resolved_root = os.path.realpath(root_path)
    resolved_candidate = os.path.realpath(candidate_path)
    return resolved_candidate == resolved_root or resolved_candidate.startswith(
        resolved_root.rstrip(os.sep) + os.sep
    )


def _list_directory_entries(dir_path, file_extensions=None):
    """
    List a directory's immediate contents, separated into subdirectories
    and files, sorted alphabetically (directories first, since that's
    the conventional order for a file browser and keeps navigation
    predictable).

    file_extensions: optional list of extensions (e.g. [".fa", ".fasta",
        ".fa.gz"]) to filter displayed files to -- directories are never
        filtered by this (the user needs to be able to navigate THROUGH
        a directory to reach a matching file inside it, regardless of
        the directory's own name). Case-insensitive matching. Pass None
        to show every file with no extension filtering.

    Returns (subdirs: list[str], files: list[str]) -- both just the
    entry NAMES (not full paths); the caller joins these with dir_path
    as needed. Silently returns ([], []) if dir_path can't be listed
    (e.g. a permissions error) rather than raising, so a single
    unreadable subdirectory doesn't crash the whole browser -- the
    caller should show an appropriate message when both lists come back
    empty and the directory wasn't actually empty (permission denied).
    """
    try:
        entries = os.listdir(dir_path)
    except OSError:
        return [], []

    subdirs = []
    files = []
    for entry in entries:
        # Skip hidden dotfiles/dotdirs (., .., .cache, .ssh, etc.) --
        # these are essentially never what a user browsing for a
        # reference file is looking for, and hiding them keeps the
        # listing focused and avoids surfacing sensitive-looking
        # dotdirs unnecessarily.
        if entry.startswith("."):
            continue
        full_path = os.path.join(dir_path, entry)
        try:
            is_dir = os.path.isdir(full_path)
            is_file = os.path.isfile(full_path)
        except OSError:
            # Permission denied stat'ing this specific entry (can
            # happen with a "/" root, e.g. certain system directories)
            # -- skip just this one entry rather than failing the
            # entire listing.
            continue
        if is_dir:
            subdirs.append(entry)
        elif is_file:
            if file_extensions is None or any(
                entry.lower().endswith(ext.lower()) for ext in file_extensions
            ):
                files.append(entry)

    return sorted(subdirs), sorted(files)


def _render_directory_navigator(key_prefix, root_dir):
    """
    Shared navigation UI used by both render_server_file_browser and
    render_server_directory_browser below -- factored out since the
    navigation mechanics are identical between the two; they only
    differ in what the final "confirm"/"select" action actually does
    (pick a FILE within the current directory, vs. confirm the CURRENT
    DIRECTORY itself).

    Provides three complementary ways to navigate, since with the
    default root now being "/" (see _default_browse_root), reaching a
    real, often deeply-nested target directory purely by clicking one
    subdirectory at a time would be slow:
      1. "Jump directly to a path" -- type a full known path and go
         there in one step (e.g. a known /scratch project directory).
      2. A LIVE-FILTERED subdirectory picker -- type a partial name in
         the filter box and the subdirectory dropdown narrows to only
         matching names as you type, then pick one and click "Go" to
         descend into it. This is the closest practical equivalent to
         shell tab-completion available within Streamlit's widget
         model -- genuine keystroke-level autocomplete isn't something
         Streamlit exposes an API for, so live-filtering-as-you-type is
         used instead to solve the same real problem (quickly narrowing
         a long subdirectory listing without needing to scroll it).
      3. Plain "go up one level", for simple step-by-step retreat.

    Returns (current_dir, subdirs, files) -- the validated current
    directory path (after handling any navigation action this run) and
    that directory's listing, ready for the caller to render its own
    file-vs-directory-specific confirmation UI underneath.
    """
    current_dir_key = f"{key_prefix}_browse_current_dir"
    filter_key = f"{key_prefix}_browse_subdir_filter"
    filter_reset_pending_key = f"{key_prefix}_browse_filter_reset_pending"

    if current_dir_key not in st.session_state:
        st.session_state[current_dir_key] = root_dir

    current_dir = st.session_state[current_dir_key]

    # Handle a filter-reset requested by a PREVIOUS run (see the "Go"
    # button handling further below) -- must happen here, before the
    # filter text_input widget is created below, since Streamlit does
    # not allow writing to a widget's session_state key after that
    # widget has already been instantiated in the SAME script run (a
    # real StreamlitAPIException confirmed via direct testing when this
    # was originally attempted immediately at the point of navigation,
    # in the same run as the filter widget's own creation). Clearing it
    # here, before the widget below is created, is always safe.
    if st.session_state.pop(filter_reset_pending_key, False):
        st.session_state[filter_key] = ""

    # Defensive re-check every render (not just on navigation) -- if
    # root_dir ever changes between reruns (e.g. a different project
    # selected), or session_state somehow held a stale/invalid path,
    # snap back to root_dir rather than trusting old state blindly.
    if not _is_path_within_root(current_dir, root_dir) or not os.path.isdir(current_dir):
        current_dir = root_dir
        st.session_state[current_dir_key] = root_dir

    # Deliberately rendered with st.code() rather than a disabled
    # st.text_input(): a stateful widget (even one with disabled=True)
    # caches its displayed value in st.session_state under its own key
    # after the FIRST render -- passing a different `value=` on a later
    # rerun is silently ignored by Streamlit unless the widget's key
    # itself changes. Confirmed via direct testing: after navigating
    # into a subdirectory, a text_input-based display kept showing the
    # OLD (root) directory even though the actual navigation state had
    # correctly updated -- a confusing, misleading UI bug. st.code() has
    # no such caching since it isn't a stateful widget at all, so it
    # always reflects current_dir exactly as computed on this run.
    st.code(current_dir, language=None)

    # --- Option 1: jump directly to a known path ---
    jump_col, jump_btn_col = st.columns([4, 1])
    with jump_col:
        jump_path = st.text_input(
            "Jump directly to a path:",
            value="", placeholder="e.g. /scratch/bioscratch/Podrab_lab",
            key=f"{key_prefix}_browse_jump_path_input",
        )
    with jump_btn_col:
        st.markdown("<div style='margin-top: 1.7rem;'></div>", unsafe_allow_html=True)
        jump_clicked = st.button("Jump", key=f"{key_prefix}_browse_jump_btn")

    if jump_clicked and jump_path:
        if _is_path_within_root(jump_path, root_dir) and os.path.isdir(jump_path):
            st.session_state[current_dir_key] = os.path.realpath(jump_path)
            st.rerun()
        else:
            st.error(
                f"⚠️ '{jump_path}' isn't accessible from here -- double "
                "check the path exists and is a directory (not a file)."
            )

    subdirs, _ = _list_directory_entries(current_dir, file_extensions=None)

    # --- Option 2: live-filtered subdirectory picker ---
    # Streamlit reruns the whole script on every keystroke in a text
    # input, so filtering the dropdown's options based on this filter
    # box's CURRENT value happens automatically on every render -- no
    # separate "apply filter" button needed, giving an experience much
    # closer to live-narrowing/autocomplete than a single static
    # dropdown of every subdirectory would.
    filter_text = st.text_input(
        "Filter subdirectories (type to narrow the list below):",
        value="", key=filter_key,
    )
    if filter_text:
        filtered_subdirs = [d for d in subdirs if filter_text.lower() in d.lower()]
    else:
        filtered_subdirs = subdirs

    if filter_text and not filtered_subdirs:
        st.caption(f"No subdirectories match '{filter_text}'.")

    nav_options = ["(stay here)"]
    if current_dir != root_dir:
        nav_options.append(".. (go up one level)")
    nav_options.extend(f"📁 {d}" for d in filtered_subdirs)

    nav_choice = st.selectbox(
        f"Navigate ({len(filtered_subdirs)} of {len(subdirs)} subdirectories shown):",
        options=nav_options,
        key=f"{key_prefix}_browse_nav_select",
    )
    if st.button("Go", key=f"{key_prefix}_browse_nav_go_btn"):
        if nav_choice == ".. (go up one level)":
            new_dir = os.path.dirname(current_dir.rstrip(os.sep)) or os.sep
        elif nav_choice.startswith("📁 "):
            new_dir = os.path.join(current_dir, nav_choice[2:])
        else:
            new_dir = current_dir

        # Re-validate before committing -- belt-and-suspenders on top
        # of the fact that new_dir was only ever built from names we
        # ourselves listed via os.listdir() above (never from free-text
        # user input), so this should always pass; the check is kept
        # anyway since it's cheap and this is a security boundary worth
        # double-checking unconditionally.
        if _is_path_within_root(new_dir, root_dir) and os.path.isdir(new_dir):
            st.session_state[current_dir_key] = new_dir
            # Request a subdirectory-filter reset for the NEXT run
            # (handled at the top of this function, before the filter
            # widget is created) rather than clearing it directly here
            # -- the filter widget has ALREADY been instantiated
            # earlier in this same run (it's created above, before this
            # navigation handling code), and Streamlit does not allow
            # writing to a widget's session_state key after it has been
            # instantiated within the same run. This avoids a filter
            # typed for the PREVIOUS directory (e.g. "cherry") silently
            # carrying over and hiding every subdirectory of the NEW
            # directory that doesn't happen to also match that same
            # text, with no obvious explanation to the user for why the
            # list looks unexpectedly short/empty.
            st.session_state[filter_reset_pending_key] = True
            st.rerun()
        else:
            st.error("⚠️ That location isn't accessible from here.")

    return current_dir


def render_server_file_browser(key_prefix, root_dir=None, file_extensions=None,
                                label="Browse for a file already on this server:"):
    """
    Render a sandboxed, navigable server-side file browser widget, and
    return the absolute path of whatever file the user has selected (or
    None if nothing is selected yet).

    key_prefix: a unique string scoping this browser instance's
        st.session_state keys and widget keys -- REQUIRED when a page
        might render more than one independent file browser (e.g. one
        for a genome FASTA and a separate one for a GTF), so their
        navigation state doesn't collide with each other.

    root_dir: the sandboxed root directory this browser is allowed to
        navigate within -- defaults to _default_browse_root() ("/", the
        filesystem root -- i.e. effectively unrestricted beyond the
        OS's own file permissions) if not given. Pass an explicit,
        narrower path here if a genuinely restricted sandbox is needed
        for a specific deployment/call site.

    file_extensions: optional list of file extensions to filter the
        file listing to (e.g. [".fa", ".fasta", ".fa.gz", ".fna"] for a
        genome FASTA picker, or [".csv", ".xlsx", ".xls"] for a
        metadata file picker). Directories are always shown regardless
        (so the user can navigate through them), only the FILE listing
        at each level is filtered. Pass None to show all files.

    label: the text shown above the current-directory display, so each
        call site can describe what kind of file the user is looking
        for (e.g. "Browse for a genome FASTA file already on this
        server:").

    Returns the selected file's absolute path (str), or None if the
    user hasn't confirmed a selection yet this run. Once a path is
    returned, the CALLER is responsible for using it directly (e.g.
    passing it straight to whatever function previously took an
    uploaded/saved file's path) -- this function never copies, reads,
    or modifies the file in any way.
    """
    if root_dir is None:
        root_dir = _default_browse_root()
    root_dir = os.path.realpath(root_dir)

    selected_file_key = f"{key_prefix}_browse_selected_file"

    st.caption(label)
    current_dir = _render_directory_navigator(key_prefix, root_dir)

    # Clear any stale selection from a previous directory the moment we
    # detect the user has navigated elsewhere since selecting -- avoids
    # a confusing state where an old selection from a different
    # directory is still shown as "selected" after navigating away.
    stale_selection = st.session_state.get(selected_file_key)
    if stale_selection and os.path.dirname(stale_selection) != current_dir:
        st.session_state.pop(selected_file_key, None)

    _, files = _list_directory_entries(current_dir, file_extensions=file_extensions)

    if files:
        file_choice = st.selectbox(
            "Files in this directory:", options=["(none selected)"] + files,
            key=f"{key_prefix}_browse_file_select",
        )
        if file_choice != "(none selected)":
            if st.button("✅ Select this file", key=f"{key_prefix}_browse_select_btn"):
                candidate = os.path.join(current_dir, file_choice)
                if _is_path_within_root(candidate, root_dir) and os.path.isfile(candidate):
                    st.session_state[selected_file_key] = candidate
                    st.rerun()
                else:
                    st.error("⚠️ That file isn't accessible from here.")
    else:
        ext_note = f" matching {', '.join(file_extensions)}" if file_extensions else ""
        st.caption(f"No files{ext_note} in this directory.")

    selected_path = st.session_state.get(selected_file_key)
    if selected_path:
        size_mb = os.path.getsize(selected_path) / (1024 * 1024)
        st.success(f"✅ Selected: `{selected_path}` ({size_mb:,.1f} MB)")
        if st.button("Clear selection", key=f"{key_prefix}_browse_clear_btn"):
            st.session_state.pop(selected_file_key, None)
            st.rerun()
        return selected_path

    return None


def render_server_directory_browser(key_prefix, root_dir=None,
                                     label="Browse for a directory already on this server:",
                                     preview_extensions=None):
    """
    Render a sandboxed, navigable server-side DIRECTORY browser widget
    -- unlike render_server_file_browser above (which selects one
    specific file), this widget lets the user navigate to and confirm
    a DIRECTORY itself, appropriate when the caller needs to work with
    an entire folder's contents (e.g. a directory already containing a
    full set of FASTQ files for many samples) rather than one file.

    key_prefix, root_dir: same meaning as in render_server_file_browser
        above (root_dir also defaults to "/" if not given).

    preview_extensions: optional list of file extensions (e.g.
        [".fastq", ".fastq.gz"]) -- if given, a live count of how many
        matching files exist directly inside the CURRENT directory
        being browsed is shown as the user navigates, so they get
        immediate feedback on whether they've found the right folder
        before confirming (e.g. "12 matching file(s) in this
        directory"). Purely informational -- has no effect on what
        gets returned; the caller is responsible for actually scanning
        the confirmed directory's contents afterward.

    Returns the confirmed directory's absolute path (str), or None if
    the user hasn't confirmed a directory yet this run. This function
    never reads, lists (beyond the optional preview count), copies, or
    modifies anything inside the returned directory -- the caller is
    fully responsible for whatever it does with the directory's
    contents afterward.
    """
    if root_dir is None:
        root_dir = _default_browse_root()
    root_dir = os.path.realpath(root_dir)

    confirmed_dir_key = f"{key_prefix}_browse_confirmed_dir"

    st.caption(label)
    current_dir = _render_directory_navigator(key_prefix, root_dir)

    if preview_extensions:
        _, files = _list_directory_entries(current_dir, file_extensions=preview_extensions)
        ext_note = ", ".join(preview_extensions)
        st.caption(f"📄 {len(files)} matching file(s) ({ext_note}) directly in this directory.")

    if st.button("✅ Use this directory", key=f"{key_prefix}_browse_confirm_dir_btn"):
        if _is_path_within_root(current_dir, root_dir) and os.path.isdir(current_dir):
            st.session_state[confirmed_dir_key] = current_dir
            st.rerun()
        else:
            st.error("⚠️ That directory isn't accessible from here.")

    confirmed_dir = st.session_state.get(confirmed_dir_key)
    if confirmed_dir:
        st.success(f"✅ Using directory: `{confirmed_dir}`")
        if st.button("Clear selection", key=f"{key_prefix}_browse_clear_dir_btn"):
            st.session_state.pop(confirmed_dir_key, None)
            st.rerun()
        return confirmed_dir

    return None
