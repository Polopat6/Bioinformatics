"""
setup_workspace.py

The "⚙️ Setup & Deployment" page -- a standalone sidebar entry (alongside
"📊 Portal Home", outside any pipeline drawer) since it's a cross-cutting
utility page rather than a pipeline step. Three sections:

1. Environment & Dependency Check -- live, on-demand checks of whether
   every Python package / CLI tool / R+Bioconductor package this app
   needs (per environment.yml, the project's single source of truth for
   dependencies -- see DEPLOYMENT.md) is actually present on THIS
   machine, right now. This exists precisely because of the real gap
   found during HPC deployment testing: code comments claimed FastQC/
   MultiQC/SRA Toolkit/several R packages were "in the Dockerfile", but
   requirements.txt only ever tracked 4 Python packages and had no way
   to express the rest -- so the gap went undiscovered until specific
   workflow steps failed downstream. This page surfaces that gap
   immediately and in one place instead.

2. Install Missing Dependencies -- (2026-08-17) once a check finds
   anything missing, offers a one-click "Install Missing Packages"
   button that launches a background `mamba install`/`conda install`
   (via deployment_manager.py) for exactly the conda package specs
   environment.yml itself declares -- so an install triggered from this
   page can never drift from what environment.yml says should be
   installed. Runs as a background process (same subprocess.Popen +
   JSON status-file polling pattern used elsewhere in this app for
   long-running work, e.g. advanced_mode_orchestrator's pipeline runs)
   with a live log tail, since a real dependency-solve + download can
   take several minutes -- long enough that blocking Streamlit's single-
   threaded script execution on it would freeze the whole page.

3. HPC Connections -- save, test, and manage SSH connection profiles to
   remote HPC clusters via hpc_manager.py. Passwords are never persisted
   (see hpc_manager.py's module docstring for the full security design);
   SSH-key or ssh-agent auth is the recommended path for any saved
   connection. This is scoped deliberately to connectivity/environment
   detection only -- NOT remote job submission/execution, which would be
   a separate, explicitly-scoped future feature building on top of a
   working, tested connection profile.
"""
import importlib.util
import shutil
import subprocess

import streamlit as st

import hpc_manager as hpc
import deployment_manager as dm

# ---------------------------------------------------------------------------
# Environment & Dependency Check
# ---------------------------------------------------------------------------
# Kept in sync BY HAND with environment.yml's dependency list -- if you add
# a package there, add its corresponding entry here too (and vice versa), so
# this page never silently drifts out of sync with the actual single source
# of truth. Each entry now also carries its OWN conda package spec (name +
# version pin, copied exactly from environment.yml) so "Install Missing
# Dependencies" below can install precisely what environment.yml declares --
# never a guessed or differently-pinned version.
#
# There is no automated way to parse environment.yml generically into "how
# do I check if this is importable/on PATH" + "what conda package installs
# it", since conda package names, Python import names, R package names, and
# CLI executable names frequently all differ from one another for the same
# logical dependency (e.g. conda package "python-kaleido" -> Python import
# "kaleido"; conda package "star" -> CLI executable "STAR").

_PYTHON_PACKAGES = {
    "streamlit": {"label": "streamlit", "conda_spec": "streamlit=1.37.*"},
    "pandas": {"label": "pandas", "conda_spec": "pandas=2.2.*"},
    "plotly": {"label": "plotly", "conda_spec": "plotly=5.22.*"},
    "kaleido": {"label": "python-kaleido (Plotly PDF/image export)", "conda_spec": "python-kaleido=0.2.*"},
    "duckdb": {"label": "duckdb", "conda_spec": "duckdb=1.0.*"},
    "openpyxl": {"label": "openpyxl (.xlsx metadata upload support)", "conda_spec": "openpyxl=3.1.*"},
    "scipy": {"label": "scipy", "conda_spec": "scipy=1.13.*"},
    "paramiko": {"label": "paramiko (SSH connections, this page's own HPC section)", "conda_spec": "paramiko=3.4.*"},
    "requests": {"label": "requests (NCBI/SRA lookups)", "conda_spec": "requests=2.32.*"},
}

_CLI_TOOLS = {
    "fastqc": {"label": "FastQC (pre/post-trim QC)", "conda_spec": "fastqc=0.12.*"},
    "multiqc": {"label": "MultiQC (combined QC reports)", "conda_spec": "multiqc=1.21.*"},
    "fastp": {"label": "fastp (adapter/quality trimming)", "conda_spec": "fastp=0.23.*"},
    "salmon": {"label": "Salmon (pseudo-alignment/quantification)", "conda_spec": "salmon=1.10.*"},
    "STAR": {"label": "STAR (splice-aware alignment)", "conda_spec": "star=2.7.11b"},
    # prefetch and fasterq-dump both come from the same sra-tools conda
    # package -- sharing the identical conda_spec here means the install
    # list de-duplication in _collect_missing_specs() below correctly
    # installs sra-tools only ONCE even if both are missing.
    "prefetch": {"label": "SRA Toolkit -- prefetch (NCBI/SRA download)", "conda_spec": "sra-tools=3.1.*"},
    "fasterq-dump": {"label": "SRA Toolkit -- fasterq-dump (NCBI/SRA download)", "conda_spec": "sra-tools=3.1.*"},
}

# One combined Rscript call checks every R/Bioconductor package at once,
# rather than spawning a separate R subprocess per package -- R's own
# startup overhead (loading the base R environment) is the dominant cost
# per invocation, so batching avoids paying that cost N times over.
#
# --- Version-pinning policy correction (2026-08-17) ---
# These conda_spec values previously included exact "=X.YY.*" patch pins
# (e.g. "bioconductor-deseq2=1.44.*") copied from an earlier draft of
# environment.yml. A real HPC install failed with mamba reporting several
# of those exact versions as "does not exist" -- Bioconductor packages are
# only rebuilt on Bioconductor's twice-yearly release cadence and not every
# package gets a version bump every cycle, so a hand-picked exact version
# can easily name a release that was simply never published. These specs
# are now UNPINNED, exactly matching environment.yml's own corrected
# entries -- letting mamba/conda's solver pick whichever real version is
# compatible with the pinned r-base, rather than every hardcoded guess
# risking the same failure mode across 13 different Bioconductor packages
# that don't share a synchronized release history.
_R_PACKAGES = {
    "DESeq2": {"label": "DESeq2 (differential expression)", "conda_spec": "bioconductor-deseq2"},
    "jsonlite": {"label": "jsonlite (DESeq2 job-spec I/O)", "conda_spec": "r-jsonlite"},
    "tximport": {"label": "tximport (Salmon transcript->gene count collapsing)", "conda_spec": "bioconductor-tximport"},
    "clusterProfiler": {"label": "clusterProfiler (bitr() ID mapping; GO/KEGG enrichment)", "conda_spec": "bioconductor-clusterprofiler"},
    "ReactomePA": {"label": "ReactomePA (Reactome pathway enrichment)", "conda_spec": "bioconductor-reactompa"},
    "GOSemSim": {"label": "GOSemSim (GO term semantic-similarity simplification)", "conda_spec": "bioconductor-gosemsim"},
}

# One entry per PRESET species in reference_manager.py's REFERENCE_CATALOG
# -- add a new organism's org.*.db package here the same day it's added to
# environment.yml AND to REFERENCE_CATALOG, so all three never drift out of
# sync with each other. Unpinned for the same reason as _R_PACKAGES above.
_R_ORGANISM_PACKAGES = {
    "org.Hs.eg.db": {"label": "org.Hs.eg.db -- human (Homo sapiens)", "conda_spec": "bioconductor-org.hs.eg.db"},
    "org.Mm.eg.db": {"label": "org.Mm.eg.db -- mouse (Mus musculus)", "conda_spec": "bioconductor-org.mm.eg.db"},
    "org.Dm.eg.db": {"label": "org.Dm.eg.db -- fly (Drosophila melanogaster)", "conda_spec": "bioconductor-org.dm.eg.db"},
    "org.Sc.sgd.db": {"label": "org.Sc.sgd.db -- yeast (Saccharomyces cerevisiae)", "conda_spec": "bioconductor-org.sc.sgd.db"},
    "org.Ce.eg.db": {"label": "org.Ce.eg.db -- roundworm (Caenorhabditis elegans)", "conda_spec": "bioconductor-org.ce.eg.db"},
    "org.Dr.eg.db": {"label": "org.Dr.eg.db -- zebrafish (Danio rerio)", "conda_spec": "bioconductor-org.dr.eg.db"},
    "org.EcK12.eg.db": {"label": "org.EcK12.eg.db -- E. coli strain K-12", "conda_spec": "bioconductor-org.eck12.eg.db"},
}


def _check_python_package(module_name):
    "True if module_name can be imported in this environment, without actually importing it (avoids side effects / slow imports)."
    return importlib.util.find_spec(module_name) is not None


def _check_cli_tool(executable_name):
    "True if executable_name is found on PATH."
    return shutil.which(executable_name) is not None


def _check_r_packages_batched(package_names, timeout=30):
    """
    Check every name in package_names in a SINGLE Rscript invocation,
    returning {package_name: bool}. If Rscript itself isn't on PATH, R is
    reported as entirely unavailable rather than failing per-package (a
    single clear reason, not N confusing ones caused by the same root
    issue).
    """
    if not shutil.which("Rscript"):
        return None  # signals "R itself not found" to the caller
    r_lines = "\n".join(
        f'cat("{name}\\t", requireNamespace("{name}", quietly=TRUE), "\\n", sep="")'
        for name in package_names
    )
    script = f"suppressWarnings({{\n{r_lines}\n}})"
    try:
        result = subprocess.run(
            ["Rscript", "-e", script],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {name: False for name in package_names}
    statuses = {}
    for line in result.stdout.strip().splitlines():
        if "\t" not in line:
            continue
        name, status = line.split("\t", 1)
        statuses[name.strip()] = status.strip().upper() == "TRUE"
    return {name: statuses.get(name, False) for name in package_names}


def _render_status_row(label, found, conda_spec):
    icon = "✅" if found else "❌"
    col1, col2 = st.columns([3, 2])
    with col1:
        st.markdown(f"{icon} {label}")
    if not found:
        with col2:
            st.caption(f"`{conda_spec}`")


def _render_environment_check():
    st.subheader("🔎 Environment & Dependency Check")
    st.markdown(
        "Checks whether every dependency this portal needs -- per "
        "`environment.yml` sitting alongside this app, this project's "
        "single source of truth for dependencies -- is actually present "
        "**on this machine, right now**. Run this after any new install "
        "(Docker, local, or HPC) to catch a missing dependency "
        "immediately, rather than discovering it only when a specific "
        "workflow step fails downstream."
    )
    if st.button("🔄 Check Environment", key="setup_env_check_btn", type="primary"):
        st.session_state["setup_env_check_ran"] = True
        st.session_state.pop("setup_missing_specs", None)  # force recompute below

    if not st.session_state.get("setup_env_check_ran"):
        st.info("Click \"Check Environment\" above to run the check.")
        return

    # Collected while rendering each row below, then used both to build the
    # "Install Missing Dependencies" section right after this one, AND
    # stashed in session_state so that section still has something to work
    # with across reruns triggered by ITS OWN buttons (e.g. clicking
    # "Refresh Install Status") without needing to silently re-run the
    # whole check every time.
    missing = {}  # conda_spec -> label (first label wins if the same spec covers multiple checks, e.g. sra-tools)

    st.markdown("**🐍 Python packages**")
    for module_name, spec in _PYTHON_PACKAGES.items():
        found = _check_python_package(module_name)
        _render_status_row(spec["label"], found, spec["conda_spec"])
        if not found:
            missing.setdefault(spec["conda_spec"], spec["label"])

    st.markdown("**🛠️ External CLI tools**")
    for tool, spec in _CLI_TOOLS.items():
        found = _check_cli_tool(tool)
        _render_status_row(spec["label"], found, spec["conda_spec"])
        if not found:
            missing.setdefault(spec["conda_spec"], spec["label"])

    st.markdown("**📊 R / Bioconductor packages**")
    all_r_specs = {**_R_PACKAGES, **_R_ORGANISM_PACKAGES}
    with st.spinner("Checking R packages (this batches every package into one Rscript call)..."):
        r_statuses = _check_r_packages_batched(list(all_r_specs.keys()))
    if r_statuses is None:
        st.error(
            "❌ `Rscript` was not found on PATH at all -- R itself doesn't "
            "appear to be installed in this environment. Re-sync against "
            "environment.yml (see DEPLOYMENT.md) to install r-base along "
            "with every R/Bioconductor package below."
        )
        missing.setdefault("r-base=4.3.*", "r-base (R itself)")
    else:
        for pkg_name, spec in _R_PACKAGES.items():
            found = r_statuses.get(pkg_name, False)
            _render_status_row(spec["label"], found, spec["conda_spec"])
            if not found:
                missing.setdefault(spec["conda_spec"], spec["label"])

        st.markdown("**🧬 Organism annotation packages** (one per preset species in `reference_manager.py`'s `REFERENCE_CATALOG`)")
        for pkg_name, spec in _R_ORGANISM_PACKAGES.items():
            found = r_statuses.get(pkg_name, False)
            _render_status_row(spec["label"], found, spec["conda_spec"])
            if not found:
                missing.setdefault(spec["conda_spec"], spec["label"])

    st.session_state["setup_missing_specs"] = missing

    st.markdown("---")
    st.caption(
        "See **DEPLOYMENT.md** for full Docker / local / HPC setup "
        "instructions. Anything missing above can be installed directly "
        "from this page -- see \"Install Missing Dependencies\" below."
    )


# ---------------------------------------------------------------------------
# Install Missing Dependencies
# ---------------------------------------------------------------------------
def _render_install_status_panel():
    status = dm.get_install_status()
    if not status:
        return
    st.markdown("**📡 Install Status**")
    if status["status"] == "running":
        st.info(f"🔄 Installing -- started {status['started_at']}...")
    elif status["status"] == "complete":
        st.success(f"✅ Install completed successfully at {status['finished_at']}.")
    elif status["status"] == "error":
        if status.get("used_fallback"):
            st.error(
                "❌ The full batch install failed, so each package was "
                "retried individually (see below for which ones succeeded)."
            )
        else:
            st.error(f"❌ Install failed (exit code {status.get('returncode')}) -- see log below.")
    st.markdown(f"Packages: `{', '.join(status.get('package_specs', []))}`")

    # Once the fallback (one-at-a-time) path has run, show each package's
    # own individual outcome -- this is what actually tells the user
    # "paramiko installed fine, only bioconductor-org.ce.eg.db failed"
    # instead of leaving every package's real status a mystery behind one
    # opaque batch-level failure message.
    package_results = status.get("package_results")
    if package_results:
        st.markdown("**Per-package results:**")
        for spec, result in package_results.items():
            icon = "✅" if result == "installed" else "❌"
            st.markdown(f"{icon} `{spec}` — {result}")

    with st.expander("📜 Install log", expanded=(status["status"] != "running")):
        log_text = dm.read_install_log()
        st.code(log_text or "(no output yet)", language="text")
    if status["status"] == "running":
        if st.button("🔄 Refresh Install Status", key="setup_install_refresh_btn"):
            st.rerun()
    else:
        st.caption(
            "Re-run \"Check Environment\" above to confirm what's now "
            "installed. If you installed a package this running app itself "
            "depends on (e.g. streamlit, pandas), **restart the app** for "
            "that change to take effect in this session."
        )


def _render_install_missing_section():
    st.subheader("📦 Install Missing Dependencies")

    if dm.is_install_in_progress():
        st.markdown(
            "An install is currently running (see status below) -- wait "
            "for it to finish before starting another."
        )
        _render_install_status_panel()
        return

    missing = st.session_state.get("setup_missing_specs")
    if missing is None:
        st.info("Run \"Check Environment\" above first to see what's missing.")
        return
    if not missing:
        st.success("✅ Nothing missing -- every checked dependency was found.")
        _render_install_status_panel()  # still show the last install's result, if any, for reference
        return

    st.warning(
        f"**{len(missing)} package(s) missing.** Clicking install below will "
        "run `mamba install` / `conda install` directly in this app's ACTIVE "
        "environment, live -- this modifies your real environment, can take "
        "several minutes (dependency solving + download), and if a package "
        "this running app itself depends on is installed (e.g. streamlit, "
        "pandas), that change won't take effect until you **restart the app**."
    )
    st.markdown("**Missing packages that will be installed:**")
    for spec, label in missing.items():
        st.markdown(f"- {label} — `{spec}`")

    exe = dm.get_conda_or_mamba_executable()
    target = dm.get_active_conda_target()
    if not exe:
        st.error(
            "❌ Neither `mamba` nor `conda` was found on PATH in this "
            "environment -- cannot install automatically here. See "
            "DEPLOYMENT.md for manual install instructions."
        )
        return
    if not target:
        st.error(
            "❌ This app doesn't appear to be running inside a conda/mamba "
            "environment (no `CONDA_PREFIX`/`CONDA_DEFAULT_ENV` detected) -- "
            "cannot determine which environment to install into. See "
            "DEPLOYMENT.md for manual install instructions."
        )
        return

    flag, value = target
    st.caption(f"Will install into: `{flag} {value}` (via `{exe}`)")

    if st.button("🚀 Install Missing Packages", key="setup_install_btn", type="primary"):
        success, message = dm.launch_install(list(missing.keys()))
        if success:
            st.session_state["setup_install_just_launched"] = True
            st.rerun()
        else:
            st.error(f"❌ {message}")

    _render_install_status_panel()


# ---------------------------------------------------------------------------
# HPC Connections
# ---------------------------------------------------------------------------
def _render_connection_test_result(success, message, remote_info):
    if success:
        st.success(f"✅ {message}")
        if remote_info:
            st.markdown("**Detected remote environment:**")
            st.markdown(f"- OS: `{remote_info.get('os_info', '—')}`")
            st.markdown(f"- Scheduler: {remote_info.get('scheduler', '—')}")
            st.markdown(f"- conda/mamba on PATH: {remote_info.get('conda_or_mamba', '—')}")
    else:
        st.error(f"❌ {message}")


def _render_existing_connections():
    connections = hpc.list_connections()
    if not connections:
        st.caption("No saved HPC connections yet -- add one below.")
        return
    st.markdown(f"**Saved Connections ({len(connections)})**")
    for conn in connections:
        name = conn["profile_name"]
        with st.expander(f"🖥️ {name} — {conn.get('username', '?')}@{conn.get('host', '?')}:{conn.get('port', 22)}", expanded=False):
            st.caption(f"Auth method: {conn.get('auth_method', '?')}")
            last_tested = conn.get("last_tested_at")
            if last_tested:
                status_icon = "✅" if conn.get("last_test_success") else "❌"
                st.caption(f"Last tested: {last_tested} — {status_icon} {conn.get('last_test_message', '')}")
            else:
                st.caption("Never tested since being saved.")

            test_password = None
            if conn.get("auth_method") == "password":
                test_password = st.text_input(
                    "Password (required each time -- never saved):",
                    type="password", key=f"setup_conn_test_pw_{name}",
                )

            col1, col2 = st.columns(2)
            with col1:
                if st.button("🔌 Test Connection", key=f"setup_conn_test_btn_{name}"):
                    with st.spinner(f"Connecting to {conn.get('host')}..."):
                        success, message, remote_info = hpc.test_connection(
                            host=conn.get("host"), port=conn.get("port", 22),
                            username=conn.get("username"), auth_method=conn.get("auth_method"),
                            key_path=conn.get("key_path"), password=test_password,
                        )
                    hpc.record_test_result(name, success, message)
                    st.session_state[f"setup_conn_last_result_{name}"] = (success, message, remote_info)
                    st.rerun()
            with col2:
                if st.button("🗑️ Delete", key=f"setup_conn_delete_btn_{name}"):
                    st.session_state[f"setup_conn_confirm_delete_{name}"] = True
                if st.session_state.get(f"setup_conn_confirm_delete_{name}"):
                    st.warning(f"Delete connection profile '{name}'? This only removes the saved profile, not anything on the remote host.")
                    cc1, cc2 = st.columns(2)
                    with cc1:
                        if st.button("Yes, delete", key=f"setup_conn_confirm_delete_yes_{name}"):
                            hpc.delete_connection(name)
                            st.session_state.pop(f"setup_conn_confirm_delete_{name}", None)
                            st.rerun()
                    with cc2:
                        if st.button("Cancel", key=f"setup_conn_confirm_delete_cancel_{name}"):
                            st.session_state.pop(f"setup_conn_confirm_delete_{name}", None)
                            st.rerun()

            last_result = st.session_state.get(f"setup_conn_last_result_{name}")
            if last_result:
                _render_connection_test_result(*last_result)


def _render_new_connection_form():
    st.markdown("**➕ Add New HPC Connection**")
    profile_name = st.text_input("Connection name:", key="setup_new_conn_name", placeholder="e.g. agamede-psu")
    col1, col2 = st.columns(2)
    with col1:
        host = st.text_input("Host:", key="setup_new_conn_host", placeholder="e.g. agamede.rc.pdx.edu")
    with col2:
        port = st.number_input("Port:", min_value=1, max_value=65535, value=22, key="setup_new_conn_port")
    username = st.text_input("Username:", key="setup_new_conn_username")

    auth_method_label = st.radio(
        "Authentication method:",
        ["🔑 SSH key file", "🔗 SSH agent (already loaded keys)", "🔒 Password (not recommended, never saved)"],
        key="setup_new_conn_auth_radio",
    )
    auth_method = {"🔑 SSH key file": "key", "🔗 SSH agent (already loaded keys)": "agent", "🔒 Password (not recommended, never saved)": "password"}[auth_method_label]

    key_path = None
    test_password = None
    if auth_method == "key":
        key_path = st.text_input(
            "Path to private key file (on this machine):",
            key="setup_new_conn_key_path", placeholder="e.g. ~/.ssh/id_ed25519",
            help="This key must already be authorized (in the remote host's ~/.ssh/authorized_keys) -- this page does not copy or generate keys for you.",
        )
    elif auth_method == "password":
        st.warning(
            "⚠️ Password auth is never saved to disk, by design -- you'll need "
            "to re-enter it every time you test or use this connection. "
            "SSH-key or agent-based auth is strongly recommended for any "
            "connection you plan to reuse."
        )
        test_password = st.text_input("Password (only used for this test, never saved):", type="password", key="setup_new_conn_password")

    if st.button("🔌 Test Connection", key="setup_new_conn_test_btn"):
        if not (profile_name and host and username):
            st.error("Please fill in connection name, host, and username before testing.")
        else:
            with st.spinner(f"Connecting to {host}..."):
                success, message, remote_info = hpc.test_connection(
                    host=host, port=port, username=username, auth_method=auth_method,
                    key_path=key_path, password=test_password,
                )
            st.session_state["setup_new_conn_last_result"] = (success, message, remote_info)

    last_result = st.session_state.get("setup_new_conn_last_result")
    if last_result:
        _render_connection_test_result(*last_result)

    st.markdown("---")
    ready = bool(profile_name) and bool(host) and bool(username) and not hpc.connection_exists(profile_name)
    if profile_name and hpc.connection_exists(profile_name):
        st.error(f"A connection named '{profile_name}' already exists -- choose a different name.")
    if not ready:
        st.info("Provide a connection name, host, and username above to save (testing first is recommended, but not required).")
        return
    if st.button("💾 Save Connection", key="setup_new_conn_save_btn", type="primary"):
        hpc.save_connection({
            "profile_name": profile_name,
            "host": host,
            "port": int(port),
            "username": username,
            "auth_method": auth_method,
            "key_path": key_path,
        })
        st.session_state["setup_conn_just_saved"] = profile_name
        for key in list(st.session_state.keys()):
            if key.startswith("setup_new_conn"):
                del st.session_state[key]
        st.rerun()


def _render_hpc_connections():
    st.subheader("🖥️ HPC Connections")
    st.markdown(
        "Save and test SSH connections to remote HPC clusters. This "
        "checks connectivity and detects the remote environment (OS, "
        "job scheduler, conda/mamba availability) -- it does **not** "
        "submit or run any remote jobs; that would be a separate feature "
        "built on top of a working, tested connection here."
    )
    if not hpc.PARAMIKO_AVAILABLE:
        st.error(
            "❌ `paramiko` isn't installed in this environment -- it's "
            "listed in `environment.yml`, so re-syncing your conda/mamba "
            "environment against that file (see DEPLOYMENT.md), or using "
            "\"Install Missing Dependencies\" above, will add it. The rest "
            "of this section won't work until it's available."
        )
        return

    just_saved = st.session_state.pop("setup_conn_just_saved", None)
    if just_saved:
        st.success(f"✅ Connection '{just_saved}' saved. You can test or delete it below.")

    _render_existing_connections()
    st.markdown("---")
    _render_new_connection_form()


def render():
    st.title("⚙️ Setup & Deployment")
    st.markdown(
        "Check whether this environment has everything the portal needs "
        "installed, install anything missing directly from here, and "
        "configure SSH connections to HPC clusters for future remote-"
        "execution features."
    )
    st.markdown("---")
    _render_environment_check()
    st.markdown("---")
    _render_install_missing_section()
    st.markdown("---")
    _render_hpc_connections()
