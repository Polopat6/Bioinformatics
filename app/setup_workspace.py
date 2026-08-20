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

--- Single-cell Phase 2 dependency-detection gap fix (2026-08-17) ---
A real reported bug: the new Cell-level QC packages (DropletUtils,
scuttle, scDblFinder, Seurat, SoupX, celda/DecontX, remotes,
DoubletFinder) were added to environment.yml but never added to this
file's hand-maintained _R_PACKAGES/_PYTHON_PACKAGES dicts -- so the
Environment Check simply never checked for them at all, meaning they
could never show up as "missing" and could therefore never be offered
for install. Fixed by adding entries for all of them below. While
auditing this gap, alevin-fry (CLI) and the Phase 3 Python analysis-
layer packages (scanpy, anndata, python-igraph, leidenalg, harmonypy,
celltypist) were ALSO found missing from detection for the identical
root-cause reason, and the stale, deprecated kaleido pin
("python-kaleido=0.2.*") was corrected to match environment.yml's own
already-unpinned entry.

--- Batched R-check crash-isolation + timeout bytes/str fixes (2026-08-17,
    later same day) ---
Two real reported bugs, both in _check_r_packages_batched(): (1) a hard
R error loading one heavy package (e.g. Seurat, or celda's rstan/
stanheaders dependency) could abort the entire batched Rscript call and
wipe out previously-printed results for OTHER, perfectly-fine packages
checked earlier in the same run -- fixed by wrapping each
requireNamespace() call in its own tryCatch(), plus flush(stdout())
after every line as defense in depth; (2) a genuine timeout on a slow-
loading package crashed with `TypeError: a bytes-like object is
required, not 'str'`, because Python's subprocess.TimeoutExpired.stdout
is raw bytes even when text=True was passed to subprocess.run() -- a
documented quirk where that decoding only applies to a successful
CompletedProcess, not the timeout exception. Fixed by making
_parse_r_check_output() defensively decode bytes input. See that
function's and _check_r_packages_batched()'s own docstrings for the
full detail on both fixes.

--- GitHub-only R package: real one-click install, not just a manual
    command (2026-08-17, later same day) ---
DoubletFinder has no conda/CRAN/Bioconductor release at all (GitHub-
only: chris-mcginnis-ucsf/DoubletFinder). Previously, if missing, this
page only ever showed the manual `remotes::install_github(...)` command
for the user to copy/run themselves. Now, if r-remotes is CONFIRMED
installed (checked via the same batched Rscript call as every other R
package), a real "📥 Install via GitHub" button is offered instead,
launching the actual install in the background via
deployment_manager.launch_github_r_install() -- the SAME background-
thread + JSON-status-file + log-tail polling pattern already used for
every conda-based install on this page (see that function's own
docstring). If r-remotes is NOT yet installed (or its own check
couldn't be confirmed), the manual command is still shown as a fallback,
along with a nudge that installing r-remotes first (via "Install
Missing Dependencies" below) will unlock the one-click button on the
next check. The install command itself is built programmatically (see
deployment_manager.build_github_r_install_command()) from a plain
"owner/repo" string rather than a hand-written shell command, so the
button's actual subprocess call and the displayed fallback command can
never silently drift apart from each other.
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
    # UNPINNED -- corrected 2026-08-17 from the stale, deprecated
    # "python-kaleido=0.2.*" pin, matching environment.yml's own fix.
    "kaleido": {"label": "python-kaleido (Plotly PDF/image export)", "conda_spec": "python-kaleido"},
    "duckdb": {"label": "duckdb", "conda_spec": "duckdb=1.0.*"},
    "openpyxl": {"label": "openpyxl (.xlsx metadata upload support)", "conda_spec": "openpyxl=3.1.*"},
    "scipy": {"label": "scipy", "conda_spec": "scipy=1.13.*"},
    "paramiko": {"label": "paramiko (SSH connections, this page's own HPC section)", "conda_spec": "paramiko=3.4.*"},
    "requests": {"label": "requests (NCBI/SRA lookups)", "conda_spec": "requests=2.32.*"},

    # --- Single-cell RNA-Seq: Python analysis layer (Phase 3) additions ---
    "scanpy": {"label": "scanpy (single-cell normalization/clustering/UMAP)", "conda_spec": "scanpy"},
    "anndata": {"label": "anndata (scanpy's underlying data structure)", "conda_spec": "anndata"},
    "igraph": {"label": "python-igraph (required by leidenalg for clustering)", "conda_spec": "python-igraph"},
    "leidenalg": {"label": "leidenalg (Leiden clustering algorithm)", "conda_spec": "leidenalg"},
    "harmonypy": {"label": "harmonypy (Harmony batch-correction algorithm)", "conda_spec": "harmonypy"},
    "celltypist": {"label": "celltypist (optional automated cell-type annotation)", "conda_spec": "celltypist"},
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
    # --- Single-cell RNA-Seq: alevin-fry (Phase 1 optional alternate) ---
    "alevin-fry": {"label": "alevin-fry (optional faster alternate to STARsolo)", "conda_spec": "alevin-fry"},
}

# One combined Rscript call checks every R/Bioconductor package at once,
# rather than spawning a separate R subprocess per package -- R's own
# startup overhead (loading the base R environment) is the dominant cost
# per invocation, so batching avoids paying that cost N times over. See
# _check_r_packages_batched() for the crash-isolation + timeout-handling
# fixes (2026-08-17) that make this batching safe even when one package's
# namespace load hard-errors or takes too long.
#
# --- Version-pinning policy correction (2026-08-17) ---
# These conda_spec values previously included exact "=X.YY.*" patch pins
# copied from an earlier draft of environment.yml. A real HPC install
# failed with mamba reporting several of those exact versions as "does
# not exist" -- Bioconductor packages are only rebuilt on Bioconductor's
# twice-yearly release cadence and not every package gets a version bump
# every cycle, so a hand-picked exact version can easily name a release
# that was simply never published. These specs are now UNPINNED, exactly
# matching environment.yml's own corrected entries.
_R_PACKAGES = {
    "DESeq2": {"label": "DESeq2 (differential expression)", "conda_spec": "bioconductor-deseq2"},
    "jsonlite": {"label": "jsonlite (DESeq2/cell-QC job-spec I/O)", "conda_spec": "r-jsonlite"},
    "tximport": {"label": "tximport (Salmon transcript->gene count collapsing)", "conda_spec": "bioconductor-tximport"},
    "clusterProfiler": {"label": "clusterProfiler (bitr() ID mapping; GO/KEGG enrichment)", "conda_spec": "bioconductor-clusterprofiler"},
    "ReactomePA": {"label": "ReactomePA (Reactome pathway enrichment)", "conda_spec": "bioconductor-reactompa"},
    "GOSemSim": {"label": "GOSemSim (GO term semantic-similarity simplification)", "conda_spec": "bioconductor-gosemsim"},
    "limma": {"label": "limma (removeBatchEffect() for DESeq2's batch-adjusted PCA view)", "conda_spec": "bioconductor-limma"},

    # --- Single-cell RNA-Seq: R-based cell-level QC (Phase 2) ---
    "DropletUtils": {"label": "DropletUtils (loads STARsolo's 10x-format MTX output)", "conda_spec": "bioconductor-dropletutils"},
    "scuttle": {"label": "scuttle (per-cell QC metrics + adaptive MAD thresholds)", "conda_spec": "bioconductor-scuttle"},
    "scDblFinder": {"label": "scDblFinder (doublet detection -- default method)", "conda_spec": "bioconductor-scdblfinder"},
    "Seurat": {"label": "Seurat (required by DoubletFinder's internal PCA/clustering step)", "conda_spec": "r-seurat"},
    "SoupX": {"label": "SoupX (ambient RNA correction -- alternative method)", "conda_spec": "r-soupx"},
    "celda": {"label": "celda (provides DecontX -- ambient RNA correction, default method)", "conda_spec": "bioconductor-celda"},
    "remotes": {"label": "remotes (needed to install DoubletFinder from GitHub -- see below)", "conda_spec": "r-remotes"},
}

# --- DoubletFinder: special case, NOT installable via conda/mamba ---
# DoubletFinder is not distributed via conda-forge, bioconda, CRAN, or
# Bioconductor at all -- confirmed GitHub-only
# (chris-mcginnis-ucsf/DoubletFinder). Its presence is still CHECKED via
# the same batched, crash-isolated requireNamespace() call as every other
# R package, but is rendered and handled separately in
# _render_environment_check() -- deliberately never added to the
# "missing" dict that feeds deployment_manager.launch_install(), since
# that function only knows how to run conda/mamba installs and there is
# no conda package for this to install.
#
# --- One-click GitHub install (2026-08-17) ---
# github_repo (a plain "owner/repo" string) is used to build the actual
# install command programmatically via
# deployment_manager.build_github_r_install_command(), for BOTH the
# real one-click button (when r-remotes is confirmed present) and the
# manual-fallback display command (when it isn't) -- see
# _render_environment_check()'s R-packages section for how this dict is
# used, and this file's own module docstring for the full rationale.
_R_GITHUB_PACKAGES = {
    "DoubletFinder": {
        "label": "DoubletFinder (optional alternate doublet-detection method)",
        "github_repo": "chris-mcginnis-ucsf/DoubletFinder",
    },
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


def _check_r_packages_batched(package_names, timeout=90):
    """
    Check every name in package_names in a SINGLE Rscript invocation,
    returning (statuses: dict or None, unreached: list[str]).

    statuses: {package_name: bool} for every package whose check line
        actually executed and printed a result -- includes packages
        confirmed present (True) AND packages confirmed absent/broken
        (False). Returns None (instead of a dict) if Rscript itself
        isn't on PATH at all.
    unreached: list of requested package_names that got NO result line
        back at all -- meaning the script terminated (crashed/timed out)
        before that package's check ever ran. Distinguishing "confirmed
        missing" from "never got checked" lets the UI show an honest
        "❓ could not be checked" instead of implying a package needs
        installing when its real status is simply unknown.

    --- Crash-isolation fix (2026-08-17) ---
    Each requireNamespace() call is wrapped in its OWN tryCatch(...,
    error = function(e) FALSE) -- requireNamespace() actually LOADS a
    package's namespace, unlike a simpler "is this installed" check, so
    a package with a broken/partial install (common for heavy compiled
    dependency chains, e.g. Seurat, or celda's rstan/stanheaders
    dependency) can throw a hard, UNCAUGHT R error rather than quietly
    returning FALSE -- which, unwrapped, would abort the ENTIRE batched
    script and (since stdout is block-buffered through a pipe) could
    wipe out ALL previously-printed results, even for packages checked
    successfully earlier in the same run. flush(stdout()) after every
    line is defense in depth against an even more catastrophic failure
    (a true process crash/segfault) a plain tryCatch cannot protect
    against.

    --- TimeoutExpired bytes/str fix (2026-08-17, later same day) ---
    A genuinely slow-loading package can legitimately exceed `timeout`.
    On a real timeout, Python's subprocess.TimeoutExpired.stdout
    attribute is populated with RAW BYTES, not decoded text, regardless
    of whether the original subprocess.run() call used text=True -- a
    documented quirk (the text=True decoding wrapper only applies to a
    successful CompletedProcess, not the exception raised on timeout).
    _parse_r_check_output() now defensively decodes bytes input (see
    its own docstring), so a timeout now correctly falls through to
    "whatever completed before the timeout is preserved, remaining
    requested packages are reported as unreached" instead of crashing.
    """
    if not shutil.which("Rscript"):
        return None, []
    r_lines = "\n".join(
        f'cat("{name}\\t", tryCatch(requireNamespace("{name}", quietly=TRUE), error = function(e) FALSE), "\\n", sep=""); flush(stdout())'
        for name in package_names
    )
    script = f"suppressWarnings({{\n{r_lines}\n}})"
    try:
        result = subprocess.run(
            ["Rscript", "-e", script],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired as e:
        # See this function's own "TimeoutExpired bytes/str fix"
        # docstring section above for why e.stdout may be bytes here
        # even though text=True was passed to subprocess.run() above --
        # _parse_r_check_output() handles either type safely.
        statuses = _parse_r_check_output(e.stdout)
        unreached = [name for name in package_names if name not in statuses]
        return statuses, unreached

    statuses = _parse_r_check_output(result.stdout)
    unreached = [name for name in package_names if name not in statuses]
    return statuses, unreached


def _parse_r_check_output(stdout_data):
    """
    Parse tab-separated 'name\\tTRUE/FALSE' lines from
    _check_r_packages_batched's R script output into a {name: bool}
    dict.

    stdout_data may be either str (the normal case, from a completed
    subprocess.run(..., text=True) call) OR bytes (the case on a
    subprocess.TimeoutExpired -- see _check_r_packages_batched's own
    "TimeoutExpired bytes/str fix" docstring section). Decoded
    defensively here, rather than only at each call site, so this
    function is safe regardless of which caller/code path passes it
    which type.
    """
    if isinstance(stdout_data, bytes):
        stdout_data = stdout_data.decode("utf-8", errors="replace")
    statuses = {}
    for line in (stdout_data or "").strip().splitlines():
        if "\t" not in line:
            continue
        name, status = line.split("\t", 1)
        statuses[name.strip()] = status.strip().upper() == "TRUE"
    return statuses


def _render_status_row(label, found, conda_spec):
    icon = "✅" if found else "❌"
    col1, col2 = st.columns([3, 2])
    with col1:
        st.markdown(f"{icon} {label}")
    if not found:
        with col2:
            st.caption(f"`{conda_spec}`")


def _render_unreached_row(label):
    "For a package whose check never completed (script crashed/timed out before reaching it) -- distinct from a confirmed-missing ❌, since its real status is genuinely unknown."
    st.markdown(f"❓ {label} — *could not be checked (see warning above)*")


def _render_github_package_row(pkg_name, spec, r_statuses, unreached, remotes_available):
    """
    Render one _R_GITHUB_PACKAGES entry -- unlike a regular R package row,
    this offers a REAL one-click install button (not just a manual
    command) when r-remotes is confirmed available, since
    remotes::install_github() is what actually performs the install.

    See this file's own module docstring, "GitHub-only R package: real
    one-click install" section, for the full rationale.
    """
    if pkg_name in unreached:
        _render_unreached_row(spec["label"])
        return

    found = r_statuses.get(pkg_name, False)
    icon = "✅" if found else "❌"
    st.markdown(f"{icon} {spec['label']}")
    if found:
        return

    display_command = " ".join(dm.build_github_r_install_command(spec["github_repo"]))

    if remotes_available:
        st.caption(
            "Not available via conda/mamba (GitHub-only) -- but `r-remotes` is installed, "
            "so this can be installed directly from here:"
        )
        st.code(display_command, language="bash")
        button_key = f"setup_github_install_btn_{pkg_name}"
        if st.button(f"📥 Install {pkg_name} via GitHub", key=button_key):
            success, message = dm.launch_github_r_install(spec["label"], spec["github_repo"])
            if success:
                st.session_state["setup_install_just_launched"] = True
                st.rerun()
            else:
                st.error(f"❌ {message}")
    else:
        st.caption(
            "Not available via conda/mamba (GitHub-only), and `r-remotes` isn't confirmed "
            "installed yet -- install `r-remotes` first (see \"Install Missing Dependencies\" "
            "below), then re-run this check to unlock a one-click install here. Until then, "
            "you can run this manually inside the environment:"
        )
        st.code(display_command, language="bash")


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
    # DoubletFinder's presence is checked in the SAME batched, crash-
    # isolated Rscript call as everything else (cheap to add to the same
    # call, avoids a second R startup) -- but see _R_GITHUB_PACKAGES' own
    # comment for why it's rendered and handled separately below.
    all_r_specs = {**_R_PACKAGES, **_R_ORGANISM_PACKAGES}
    all_r_check_names = list(all_r_specs.keys()) + list(_R_GITHUB_PACKAGES.keys())
    with st.spinner("Checking R packages (this batches every package into one Rscript call)..."):
        r_statuses, unreached = _check_r_packages_batched(all_r_check_names)
    if r_statuses is None:
        st.error(
            "❌ `Rscript` was not found on PATH at all -- R itself doesn't "
            "appear to be installed in this environment. Re-sync against "
            "environment.yml (see DEPLOYMENT.md) to install r-base along "
            "with every R/Bioconductor package below."
        )
        missing.setdefault("r-base=4.3.*", "r-base (R itself)")
    else:
        if unreached:
            st.warning(
                f"⚠️ {len(unreached)} R package check(s) could not complete -- the R process "
                "likely crashed, timed out, or errored while loading one of the packages below "
                "(this can happen with slow/heavy-to-load packages like Seurat or celda on a "
                "busy or memory-constrained machine, or a broken/partial install). Packages "
                "below marked ❓ have an UNKNOWN status, not a confirmed-missing one -- try "
                "running the check again, or check that package individually with "
                "`Rscript -e 'library(<package>)'` to see the actual error."
            )

        for pkg_name, spec in _R_PACKAGES.items():
            if pkg_name in unreached:
                _render_unreached_row(spec["label"])
                continue
            found = r_statuses.get(pkg_name, False)
            _render_status_row(spec["label"], found, spec["conda_spec"])
            if not found:
                missing.setdefault(spec["conda_spec"], spec["label"])

        st.markdown("**🧬 Organism annotation packages** (one per preset species in `reference_manager.py`'s `REFERENCE_CATALOG`)")
        for pkg_name, spec in _R_ORGANISM_PACKAGES.items():
            if pkg_name in unreached:
                _render_unreached_row(spec["label"])
                continue
            found = r_statuses.get(pkg_name, False)
            _render_status_row(spec["label"], found, spec["conda_spec"])
            if not found:
                missing.setdefault(spec["conda_spec"], spec["label"])

        # --- DoubletFinder: rendered separately, never added to `missing` ---
        # See _R_GITHUB_PACKAGES' own comment above for why this can't go
        # through the SAME conda-install path as everything else -- but
        # DOES now offer its own real one-click install path when
        # r-remotes is confirmed available (see _render_github_package_row).
        remotes_available = ("remotes" not in unreached) and r_statuses.get("remotes", False)
        st.markdown("**🔀 GitHub-only R packages** (not installable via conda/mamba)")
        for pkg_name, spec in _R_GITHUB_PACKAGES.items():
            _render_github_package_row(pkg_name, spec, r_statuses, unreached, remotes_available)

    st.session_state["setup_missing_specs"] = missing

    st.markdown("---")
    st.caption(
        "See **DEPLOYMENT.md** for full Docker / local / HPC setup "
        "instructions. Anything missing above can be installed directly "
        "from this page -- see \"Install Missing Dependencies\" below "
        "(GitHub-only packages have their own install option shown "
        "inline above, once `r-remotes` is available)."
    )


# ---------------------------------------------------------------------------
# Install Missing Dependencies
# ---------------------------------------------------------------------------
def _render_install_status_panel():
    status = dm.get_install_status()
    if not status:
        return
    st.markdown("**📡 Install Status**")
    # Distinguish a GitHub-based R package install from a conda/mamba
    # install -- same status dict shape otherwise (status/started_at/
    # finished_at/returncode all mean the same thing for both), just a
    # different label so it's clear at a glance which kind just ran.
    if status.get("install_type") == "github_r_package":
        st.caption("Install type: 🔀 GitHub R package (via `remotes::install_github()`)")
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
            "for it to finish before starting another. This includes any "
            "GitHub-based R package install started above, since it shares "
            "the same lock as a conda/mamba install."
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
