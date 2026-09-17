"""
setup_workspace.py

The "⚙️ Setup & Deployment" page -- a standalone sidebar entry (alongside
"📊 Portal Home", outside any pipeline drawer) since it's a cross-cutting
utility page rather than a pipeline step. Four sections:

1. Environment & Dependency Check -- live, on-demand checks of whether
   every Python package / CLI tool / R+Bioconductor package this app
   needs (per environment.yml, the project's single source of truth for
   dependencies -- see DEPLOYMENT.md) is actually present on THIS
   machine, right now.

2. Install Missing Dependencies -- once a check finds anything missing,
   offers a selectable, per-package "Install Missing Packages" flow that
   launches a background `mamba install`/`conda install` (via
   deployment_manager.py) for exactly the conda package specs
   environment.yml itself declares.

3. eggNOG-mapper Database Setup -- a large (~49GB), one-time, ADMIN-GATED
   shared-resource download for eggNOG-mapper's orthology database (see
   "Admin-gated eggNOG database setup" section below for why this is now
   role-gated, not just UI-separated).

4. HPC Connections -- save, test, and manage SSH connection profiles to
   remote HPC clusters via hpc_manager.py.

--- Permission-gated eggNOG database setup (2026-08-24, updated for
    custom roles) ---
This ~49GB download is a large, slow, disk-heavy, shared-resource action
with real potential to disrupt a shared HPC allocation's available
scratch space for every other user/project if triggered casually or by
mistake -- gated behind the "manage_eggnog_database" permission (see
auth_manager.PERMISSION_CATALOG), on top of (not instead of) the
existing disk-space pre-flight check and confirmation checkbox. This
mirrors project_manager.py's own identical reasoning for gating
permanent project deletion behind a specific permission rather than a
fixed admin-only check: the existing safeguards protect against a
PERMITTED user accidentally triggering this by mistake; the permission
check protects against an unpermitted user being ABLE to trigger it at
all -- while still letting an admin delegate specifically THIS
capability to a custom role without granting every other admin-only
action too. As with project_manager.py's own delete flow,
_render_eggnog_download_controls() independently re-checks
auth.has_permission("manage_eggnog_database") at its own top as defense
in depth.

--- Permission-gated Install Dependencies + HPC Connections
    (2026-08-24) ---
These two sections were previously completely UNGATED (any logged-in
session could use them) -- now checked against "install_dependencies"
and "manage_hpc_connections" respectively. The built-in "tech" role is
seeded with BOTH of these permissions by default (see auth_manager.py's
own _DEFAULT_TECH_PERMISSIONS), specifically so introducing this check
does not silently take away access an ordinary user already had before
this permission system existed -- an admin can still choose to remove
either permission from the tech role (or any custom role) afterward if
tighter control is actually wanted.

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

--- Per-package isolation redesign (2026-08-23) ---
A real HPC deployment run surfaced a real bug: the previous single-
shared-Rscript-call batching design meant a hang on one heavy package
(e.g. Seurat) exhausted the ENTIRE batch's one shared timeout, and
every package listed AFTER the stuck one never got checked at all --
confirmed from a real run where every package from Seurat onward came
back "could not be checked", a clean cutoff exactly at the first heavy
package in the list. _check_r_packages_batched() was replaced with
_check_r_packages_isolated(): each package is checked via its OWN
subprocess call with its OWN independent timeout, run CONCURRENTLY via
a small thread pool (max_workers=3), so a hang or crash on ANY one
package can only ever affect that single package's own result.

--- Phase 3.9 compositional analysis dependency added (2026-08-23) ---
Added "speckle" to _R_PACKAGES -- powers Step 10 (Compositional
Analysis)'s recommended-default propeller method.

--- Selective/opt-in install for large or narrow-use packages
    (2026-08-23) ---
"Install Missing Dependencies" now renders each missing package as its
own individually-selectable checkbox -- checked by default for
ordinary, broadly-needed packages, but UNCHECKED by default for large/
narrow-use packages (eggnog-mapper + all 7 org.*.eg.db packages, see
_OPTIONAL_INSTALL_INFO) since most installs don't need every one.

--- Checkbox-triggered re-check bug fix (2026-08-23) ---
_render_environment_check() previously re-ran the ENTIRE (slow) check
on every Streamlit rerun, including reruns triggered by unrelated
widgets elsewhere on the page (e.g. an install checkbox). Fixed by
splitting the actual checking work into _run_environment_check(),
called ONLY when "Check Environment" is clicked, with results cached
in session_state; _render_environment_check() now purely renders from
that cache.

--- GitHub-only R package: real one-click install (2026-08-17) ---
DoubletFinder has no conda/CRAN/Bioconductor release at all (GitHub-
only). If r-remotes is confirmed installed, a real "📥 Install via
GitHub" button is offered, launching the install via
deployment_manager.launch_github_r_install().

--- Single-cell original-format BAM recovery dependency added
    (2026-08-24) ---
Added "bamtofastq" to _CLI_TOOLS -- see that entry's own inline comment
for the full rationale (recovering usable FASTQ, with intact barcode/
UMI sequences, from an original-format 10x BAM for SRA runs whose
standard FASTQ extraction lacks that data entirely).
"""
import concurrent.futures
import importlib.util
import os
import shutil
import subprocess

import streamlit as st

import auth_manager as auth
import hpc_manager as hpc
import deployment_manager as dm
import eggnog_manager as egm
import whitelist_manager as wlm
import reference_manager as rm
import app_paths

# ---------------------------------------------------------------------------
# Environment & Dependency Check
# ---------------------------------------------------------------------------
_PYTHON_PACKAGES = {
    "streamlit": {"label": "streamlit", "conda_spec": "streamlit=1.37.*"},
    "pandas": {"label": "pandas", "conda_spec": "pandas=2.2.*"},
    "plotly": {"label": "plotly", "conda_spec": "plotly=5.22.*"},
    "kaleido": {"label": "python-kaleido (Plotly PDF/image export)", "conda_spec": "python-kaleido"},
    "duckdb": {"label": "duckdb", "conda_spec": "duckdb=1.0.*"},
    "openpyxl": {"label": "openpyxl (.xlsx metadata upload support)", "conda_spec": "openpyxl=3.1.*"},
    "scipy": {"label": "scipy", "conda_spec": "scipy=1.13.*"},
    "paramiko": {"label": "paramiko (SSH connections, this page's own HPC section)", "conda_spec": "paramiko=3.4.*"},
    "requests": {"label": "requests (NCBI/SRA lookups)", "conda_spec": "requests=2.32.*"},
    "scanpy": {"label": "scanpy (single-cell normalization/clustering/UMAP)", "conda_spec": "scanpy"},
    "anndata": {"label": "anndata (scanpy's underlying data structure)", "conda_spec": "anndata"},
    "igraph": {"label": "python-igraph (required by leidenalg for clustering)", "conda_spec": "python-igraph"},
    "leidenalg": {"label": "leidenalg (Leiden clustering algorithm)", "conda_spec": "leidenalg"},
    "harmonypy": {"label": "harmonypy (Harmony batch-correction algorithm)", "conda_spec": "harmonypy"},
    "celltypist": {"label": "celltypist (optional automated cell-type annotation)", "conda_spec": "celltypist"},
    "sklearn": {"label": "scikit-learn (3D t-SNE fallback, called directly)", "conda_spec": "scikit-learn"},
}

_CLI_TOOLS = {
    "fastqc": {"label": "FastQC (pre/post-trim QC)", "conda_spec": "fastqc=0.12.*"},
    "multiqc": {"label": "MultiQC (combined QC reports)", "conda_spec": "multiqc=1.21.*"},
    "fastp": {"label": "fastp (adapter/quality trimming)", "conda_spec": "fastp=0.23.*"},
    "salmon": {"label": "Salmon (pseudo-alignment/quantification)", "conda_spec": "salmon=1.10.*"},
    "STAR": {"label": "STAR (splice-aware alignment)", "conda_spec": "star=2.7.11b"},
    "prefetch": {"label": "SRA Toolkit -- prefetch (NCBI/SRA download)", "conda_spec": "sra-tools=3.1.*"},
    "fasterq-dump": {"label": "SRA Toolkit -- fasterq-dump (NCBI/SRA download)", "conda_spec": "sra-tools=3.1.*"},
    "alevin-fry": {"label": "alevin-fry (optional faster alternate to STARsolo)", "conda_spec": "alevin-fry"},
    "emapper.py": {"label": "emapper.py (eggNOG-mapper -- orthology-based annotation)", "conda_spec": "eggnog-mapper"},
    "download_eggnog_data.py": {"label": "download_eggnog_data.py (eggNOG-mapper's database-download script)", "conda_spec": "eggnog-mapper"},
    # --- Single-cell RNA-Seq: original-format BAM recovery (2026-08-24) ---
    # Provides the `bamtofastq` executable, used by
    # single_cell/sc_sra_manager.py's run_bamtofastq() to recover usable
    # FASTQ (with intact cell barcode/UMI sequences) from an original-
    # format 10x BAM, for SRA runs whose standard FASTQ extraction has no
    # usable barcode/UMI data at all (confirmed real motivating case:
    # Kang et al. 2018 / GSE96583, a v1-chemistry-era 10x deposition
    # where only the cDNA read was ever uploaded to SRA as FASTQ). See
    # that module's own docstring, "Original-format BAM recovery", for
    # the full rationale, including why this only helps for the subset
    # of runs that still have a free, direct-HTTP original-format copy
    # rather than requiring the user's own paid AWS/GCP cloud account.
    "bamtofastq": {"label": "bamtofastq (10x original-format BAM -> FASTQ recovery)", "conda_spec": "10x_bamtofastq"},
    "pigz": {"label": "pigz (parallel gzip -- fast FASTQ compression)","conda_spec": "pigz",},
    "seqkit": {"label": "seqkit (fast R1/R2 resync in single-cell trimming)","conda_spec": "seqkit",},
    # --- Two confirmed detection gaps (2026-09-16), same root cause as
    # the 2026-08-17 Phase 2 gap in this file's own docstring: both are
    # real, hard dependencies declared in environment.yml but never
    # added to this hand-maintained dict, so the Environment Check could
    # never report them missing or offer to install them.
    #
    # gffread drives reference_manager._download_ensembl_annotation()'s
    # GFF3->GTF fallback path (used whenever a species' direct Ensembl
    # GTF is unavailable) AND eggnog_manager's protein extraction. Its
    # absence currently surfaces as a confusing mid-reference-download
    # failure -- reference_manager's own error message even says
    # "gffread is part of this project's environment.yml" -- rather than
    # a clean ❌ row here, before anything is attempted.
    "gffread": {"label": "gffread (GFF3->GTF conversion; transcript/protein extraction)", "conda_spec": "gffread"},
    # samtools underpins bulk_bam_manager.py's entire BAM->FASTQ
    # recovery path (`samtools fastq`), which is a real, offered
    # ingestion option in the Bulk RNA-Seq workspace.
    "samtools": {"label": "samtools (BAM -> FASTQ recovery, bulk pipeline)", "conda_spec": "samtools"},
}

_R_PACKAGES = {
    "DESeq2": {"label": "DESeq2 (differential expression)", "conda_spec": "bioconductor-deseq2"},
    "jsonlite": {"label": "jsonlite (DESeq2/cell-QC job-spec I/O)", "conda_spec": "r-jsonlite"},
    "tximport": {"label": "tximport (Salmon transcript->gene count collapsing)", "conda_spec": "bioconductor-tximport"},
    "clusterProfiler": {"label": "clusterProfiler (bitr() ID mapping; GO/KEGG enrichment)", "conda_spec": "bioconductor-clusterprofiler"},
    "ReactomePA": {"label": "ReactomePA (Reactome pathway enrichment)", "conda_spec": "bioconductor-reactompa"},
    "GOSemSim": {"label": "GOSemSim (GO term semantic-similarity simplification)", "conda_spec": "bioconductor-gosemsim"},
    "limma": {"label": "limma (removeBatchEffect() for DESeq2's batch-adjusted PCA view)", "conda_spec": "bioconductor-limma"},
    "DropletUtils": {"label": "DropletUtils (loads STARsolo's 10x-format MTX output)", "conda_spec": "bioconductor-dropletutils"},
    "scuttle": {"label": "scuttle (per-cell QC metrics + adaptive MAD thresholds)", "conda_spec": "bioconductor-scuttle"},
    "scDblFinder": {"label": "scDblFinder (doublet detection -- default method)", "conda_spec": "bioconductor-scdblfinder"},
    "Seurat": {"label": "Seurat (required by DoubletFinder's internal PCA/clustering step)", "conda_spec": "r-seurat"},
    "SoupX": {"label": "SoupX (ambient RNA correction -- alternative method)", "conda_spec": "r-soupx"},
    "celda": {"label": "celda (provides DecontX -- ambient RNA correction, default method)", "conda_spec": "bioconductor-celda"},
    "remotes": {"label": "remotes (needed to install DoubletFinder from GitHub -- see below)", "conda_spec": "r-remotes"},
    "speckle": {"label": "speckle (Step 10 Compositional Analysis -- propeller method, recommended default)", "conda_spec": "bioconductor-speckle"},
}

_R_GITHUB_PACKAGES = {
    "DoubletFinder": {
        "label": "DoubletFinder (optional alternate doublet-detection method)",
        "github_repo": "chris-mcginnis-ucsf/DoubletFinder",
    },
}

_R_ORGANISM_PACKAGES = {
    "org.Hs.eg.db": {"label": "org.Hs.eg.db -- human (Homo sapiens)", "conda_spec": "bioconductor-org.hs.eg.db"},
    "org.Mm.eg.db": {"label": "org.Mm.eg.db -- mouse (Mus musculus)", "conda_spec": "bioconductor-org.mm.eg.db"},
    "org.Dm.eg.db": {"label": "org.Dm.eg.db -- fly (Drosophila melanogaster)", "conda_spec": "bioconductor-org.dm.eg.db"},
    "org.Sc.sgd.db": {"label": "org.Sc.sgd.db -- yeast (Saccharomyces cerevisiae)", "conda_spec": "bioconductor-org.sc.sgd.db"},
    "org.Ce.eg.db": {"label": "org.Ce.eg.db -- roundworm (Caenorhabditis elegans)", "conda_spec": "bioconductor-org.ce.eg.db"},
    "org.Dr.eg.db": {"label": "org.Dr.eg.db -- zebrafish (Danio rerio)", "conda_spec": "bioconductor-org.dr.eg.db"},
    "org.EcK12.eg.db": {"label": "org.EcK12.eg.db -- E. coli strain K-12", "conda_spec": "bioconductor-org.eck12.eg.db"},
}

_OPTIONAL_INSTALL_INFO = {
    "eggnog-mapper": (
        "Pulls in a fairly large bioinformatics toolchain (DIAMOND, HMMER, MMseqs2, "
        "Prodigal) -- only needed for orthology-based functional annotation on "
        "non-model organisms. Skip this if you don't work with non-model organisms. "
        "Note: this only controls the tool/package itself -- the much larger ~49GB "
        "eggNOG *database* remains its own separate, admin-gated action further "
        "down this page regardless of this checkbox."
    ),
    "bioconductor-org.hs.eg.db": "Human (Homo sapiens) gene annotation package -- only needed if you work with human data.",
    "bioconductor-org.mm.eg.db": "Mouse (Mus musculus) gene annotation package -- only needed if you work with mouse data.",
    "bioconductor-org.dm.eg.db": "Fly (Drosophila melanogaster) gene annotation package -- only needed if you work with fly data.",
    "bioconductor-org.sc.sgd.db": "Yeast (Saccharomyces cerevisiae) gene annotation package -- only needed if you work with yeast data.",
    "bioconductor-org.ce.eg.db": "Roundworm (Caenorhabditis elegans) gene annotation package -- only needed if you work with roundworm data.",
    "bioconductor-org.dr.eg.db": "Zebrafish (Danio rerio) gene annotation package -- only needed if you work with zebrafish data.",
    "bioconductor-org.eck12.eg.db": "E. coli strain K-12 gene annotation package -- only needed if you work with E. coli data.",
}


def _check_python_package(module_name):
    return importlib.util.find_spec(module_name) is not None


def _check_cli_tool(executable_name):
    return shutil.which(executable_name) is not None


def _check_one_r_package(package_name, timeout):
    script = (
        f'cat(tryCatch(requireNamespace("{package_name}", quietly = TRUE), '
        f'error = function(e) FALSE))'
    )
    try:
        result = subprocess.run(
            ["Rscript", "-e", script],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return package_name, None, "timeout"
    except Exception as e:
        return package_name, None, str(e)

    output = (result.stdout or "").strip().upper()
    return package_name, (output == "TRUE"), None


def _check_r_packages_isolated(package_names, per_package_timeout=120, max_workers=3):
    if not shutil.which("Rscript"):
        return None, []

    statuses = {}
    unreached = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_check_one_r_package, name, per_package_timeout): name
            for name in package_names
        }
        for future in concurrent.futures.as_completed(futures):
            name = futures[future]
            try:
                _, found, error = future.result()
            except Exception:
                unreached.append(name)
                continue
            if error is not None:
                unreached.append(name)
            else:
                statuses[name] = found

    return statuses, unreached


def _render_status_row(label, found, conda_spec):
    icon = "✅" if found else "❌"
    col1, col2 = st.columns([3, 2])
    with col1:
        st.markdown(f"{icon} {label}")
    if not found:
        with col2:
            st.caption(f"`{conda_spec}`")


def _render_unreached_row(label):
    st.markdown(f"❓ {label} — *could not be checked in time (see warning above) -- try again, or check it individually with `Rscript -e 'library(<package>)'`*")


def _render_github_package_row(pkg_name, spec, r_statuses, unreached, remotes_available):
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


def _run_environment_check():
    """
    Perform every ACTUAL dependency check exactly once, returning a
    plain results dict for _render_environment_check() to render from.
    This is the ONLY place that does real checking work -- see this
    file's own module docstring, "Checkbox-triggered re-check bug fix".
    """
    python_found = {name: _check_python_package(name) for name in _PYTHON_PACKAGES}
    cli_found = {name: _check_cli_tool(name) for name in _CLI_TOOLS}

    all_r_specs = {**_R_PACKAGES, **_R_ORGANISM_PACKAGES}
    all_r_check_names = list(all_r_specs.keys()) + list(_R_GITHUB_PACKAGES.keys())
    r_statuses, r_unreached = _check_r_packages_isolated(all_r_check_names)

    return {
        "python_found": python_found,
        "cli_found": cli_found,
        "r_statuses": r_statuses,
        "r_unreached": r_unreached,
    }


def _render_environment_check():
    st.subheader("🔎 Environment & Dependency Check")
    st.markdown(
        "Checks whether every dependency this portal needs -- per "
        "`environment.yml` sitting alongside this app -- is actually "
        "present **on this machine, right now**."
    )
    if st.button("🔄 Check Environment", key="setup_env_check_btn", type="primary"):
        with st.spinner(
            "Checking Python packages, CLI tools, and R packages (each R package is "
            "checked in its own isolated process -- this may take a little while)..."
        ):
            st.session_state["setup_env_results"] = _run_environment_check()
        st.session_state["setup_env_check_ran"] = True

    if not st.session_state.get("setup_env_check_ran"):
        st.info("Click \"Check Environment\" above to run the check.")
        return

    results = st.session_state.get("setup_env_results")
    if results is None:
        st.info("Click \"Check Environment\" above to run the check.")
        return

    missing = {}

    st.markdown("**🐍 Python packages**")
    for module_name, spec in _PYTHON_PACKAGES.items():
        found = results["python_found"].get(module_name, False)
        _render_status_row(spec["label"], found, spec["conda_spec"])
        if not found:
            missing.setdefault(spec["conda_spec"], spec["label"])

    st.markdown("**🛠️ External CLI tools**")
    for tool, spec in _CLI_TOOLS.items():
        found = results["cli_found"].get(tool, False)
        _render_status_row(spec["label"], found, spec["conda_spec"])
        if not found:
            missing.setdefault(spec["conda_spec"], spec["label"])

    st.markdown("**📊 R / Bioconductor packages**")
    r_statuses = results["r_statuses"]
    unreached = results["r_unreached"]
    if r_statuses is None:
        st.error("❌ `Rscript` was not found on PATH at all.")
        missing.setdefault("r-base=4.3.*", "r-base (R itself)")
    else:
        if unreached:
            st.warning(
                f"⚠️ {len(unreached)} R package check(s) could not complete in time -- each "
                "package is checked in its own fully independent process, so this means "
                "each package listed below with a ❓ genuinely, individually timed out."
            )

        for pkg_name, spec in _R_PACKAGES.items():
            if pkg_name in unreached:
                _render_unreached_row(spec["label"])
                continue
            found = r_statuses.get(pkg_name, False)
            _render_status_row(spec["label"], found, spec["conda_spec"])
            if not found:
                missing.setdefault(spec["conda_spec"], spec["label"])

        st.markdown("**🧬 Organism annotation packages**")
        for pkg_name, spec in _R_ORGANISM_PACKAGES.items():
            if pkg_name in unreached:
                _render_unreached_row(spec["label"])
                continue
            found = r_statuses.get(pkg_name, False)
            _render_status_row(spec["label"], found, spec["conda_spec"])
            if not found:
                missing.setdefault(spec["conda_spec"], spec["label"])

        remotes_available = ("remotes" not in unreached) and r_statuses.get("remotes", False)
        st.markdown("**🔀 GitHub-only R packages**")
        for pkg_name, spec in _R_GITHUB_PACKAGES.items():
            _render_github_package_row(pkg_name, spec, r_statuses, unreached, remotes_available)

    st.session_state["setup_missing_specs"] = missing
    st.markdown("---")
    st.caption("See **DEPLOYMENT.md** for full setup instructions.")


def _render_install_status_panel():
    status = dm.get_install_status()
    if not status:
        return
    st.markdown("**📡 Install Status**")
    if status.get("install_type") == "github_r_package":
        st.caption("Install type: 🔀 GitHub R package (via `remotes::install_github()`)")
    if status["status"] == "running":
        st.info(f"🔄 Installing -- started {status['started_at']}...")
    elif status["status"] == "complete":
        st.success(f"✅ Install completed successfully at {status['finished_at']}.")
    elif status["status"] == "error":
        if status.get("used_fallback"):
            st.error("❌ The full batch install failed, so each package was retried individually.")
        else:
            st.error(f"❌ Install failed (exit code {status.get('returncode')}) -- see log below.")
    st.markdown(f"Packages: `{', '.join(status.get('package_specs', []))}`")

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
        st.caption("Re-run \"Check Environment\" above to confirm what's now installed.")


def _render_install_missing_section():
    st.subheader("📦 Install Missing Dependencies")

    if not auth.has_permission("install_dependencies"):
        st.info(
            "🔒 Your current role does not include permission to install dependencies. "
            "Contact your lab's admin if a package needs to be installed."
        )
        return

    if dm.is_install_in_progress():
        st.markdown("An install is currently running -- wait for it to finish before starting another.")
        _render_install_status_panel()
        return

    missing = st.session_state.get("setup_missing_specs")
    if missing is None:
        st.info("Run \"Check Environment\" above first to see what's missing.")
        return
    if not missing:
        st.success("✅ Nothing missing -- every checked dependency was found.")
        _render_install_status_panel()
        return

    st.warning(f"**{len(missing)} package(s) missing.**")
    st.markdown("**Select which missing packages to install:**")
    st.caption(
        "Everything below is checked by default, EXCEPT eggNOG-mapper and "
        "per-organism annotation packages (org.*.eg.db) -- those are left "
        "**unchecked** by default since they can be sizable and most "
        "installs don't need every one of them."
    )

    selected_specs = []
    for spec, label in missing.items():
        optional_note = _OPTIONAL_INSTALL_INFO.get(spec)
        is_optional = optional_note is not None
        checkbox_key = f"setup_missing_checkbox_{spec}"
        checkbox_label = f"{label} — `{spec}`" + ("  *(optional)*" if is_optional else "")
        checked = st.checkbox(checkbox_label, value=not is_optional, key=checkbox_key)
        if is_optional:
            st.caption(f"↳ {optional_note}")
        if checked:
            selected_specs.append(spec)

    n_skipped = len(missing) - len(selected_specs)
    if n_skipped:
        st.caption(f"ℹ️ {n_skipped} package(s) currently unchecked and will be skipped.")

    if not selected_specs:
        st.info("Nothing selected above -- check at least one package to enable installing.")
        return

    exe = dm.get_conda_or_mamba_executable()
    target = dm.get_active_conda_target()
    if not exe:
        st.error("❌ Neither `mamba` nor `conda` was found on PATH in this environment.")
        return
    if not target:
        st.error("❌ This app doesn't appear to be running inside a conda/mamba environment.")
        return

    flag, value = target
    st.caption(f"Will install into: `{flag} {value}` (via `{exe}`)")

    if st.button(f"🚀 Install Selected Packages ({len(selected_specs)})", key="setup_install_btn", type="primary"):
        success, message = dm.launch_install(selected_specs)
        if success:
            st.session_state["setup_install_just_launched"] = True
            st.rerun()
        else:
            st.error(f"❌ {message}")

    _render_install_status_panel()


# ---------------------------------------------------------------------------
# eggNOG-mapper Database Setup (admin-gated shared resource, 2026-08-24)
# ---------------------------------------------------------------------------

def _render_eggnog_database_setup():
    st.subheader("🧬 eggNOG-mapper Database Setup (admin action)")
    st.markdown(
        "eggNOG-mapper provides orthology-based functional annotation for ANY organism. "
        "This is a **large (~49GB), one-time, shared** resource -- like a reference genome "
        "download, but bigger -- so it's restricted to administrator accounts."
    )

    if not auth.has_permission("manage_eggnog_database"):
        st.info(
            "🔒 Your current role does not include permission to set up or modify "
            "the eggNOG database. If you need this, ask your lab's admin -- it "
            "only needs to be done once, and every project/user shares the same "
            "installed database afterward."
        )
        return

    if not egm.eggnog_mapper_available() or not egm.download_eggnog_data_script_available():
        st.error(
            "❌ The `eggnog-mapper` package isn't installed yet -- install it first "
            "(see \"Install Missing Dependencies\" above), then return here to set up "
            "its database."
        )
        return

    db_dir = app_paths.data_path("shared_resources", "eggnog_database")

    if egm.eggnog_database_is_installed(db_dir):
        st.success(f"✅ The eggNOG database is already installed at `{db_dir}`.")
        with st.expander("🔁 Re-download / repair (advanced)"):
            st.caption("Only do this if you have a specific reason to believe the existing database is corrupted or incomplete.")
            _render_eggnog_download_controls(db_dir, force=True)
        return

    st.warning(
        "⚠️ The eggNOG database has not been installed on this system yet -- this download "
        "is approximately 49 GB (45 GB core annotation database + 4 GB DIAMOND search database)."
    )
    _render_eggnog_download_controls(db_dir, force=False)


def _render_eggnog_download_controls(db_dir, force):
    """
    Independently re-checks auth.has_permission("manage_eggnog_database")
    at its own top, in ADDITION to _render_eggnog_database_setup() above
    already refusing to call this function at all for an unpermitted
    session -- defense in depth, exactly mirroring project_manager.py's
    own identical pattern for its own permission-gated action (see that
    module's own docstring for the full rationale on why BOTH layers
    matter).
    """
    if not auth.has_permission("manage_eggnog_database"):
        st.error("⚠️ Your current role does not include this permission.")
        return

    disk_check = egm.check_disk_space(db_dir)
    st.markdown(disk_check["message"])

    if not disk_check["sufficient"]:
        st.error(
            "❌ Insufficient disk space detected -- proceeding is very likely to fail "
            "partway through."
        )
        confirmed = st.checkbox(
            "I understand the risk and want to attempt this anyway", key="eggnog_db_force_confirm",
        )
        if not confirmed:
            return

    button_label = "🔄 Re-download eggNOG Database" if force else "📥 Download eggNOG Database (~49 GB, one-time)"
    if st.button(button_label, key="eggnog_db_download_btn", type="primary"):
        with st.spinner("Downloading eggNOG database... this can take a long time (large download)."):
            success, message, built = egm.download_eggnog_database(
                db_dir, ensure_shared_resource_fn=rm.ensure_shared_resource,
            )
        if success:
            st.success(f"✅ {message}")
            st.rerun()
        else:
            st.error(f"❌ {message}")
def _render_whitelist_setup():
    st.subheader("🧬 Barcode Whitelist Setup")
    st.markdown(
        "Barcode whitelist files (10x Genomics' official per-chemistry cell-barcode "
        "inclusion lists) are a **shared resource** -- like a reference genome, but much "
        "smaller -- placed once here, then reused by every project/user on this server "
        "for the matching chemistry. Missing a whitelist doesn't block alignment, but it "
        "does mean chemistry selection in Step 1 can only be confirmed by read length "
        "alone, rather than by directly checking real barcode sequences against a known "
        "list."
    )

    if not auth.has_permission("install_dependencies"):
        st.info(
            "🔒 Your current role does not include permission to install dependencies "
            "(which also covers downloading whitelist files). Contact your lab's admin "
            "if a whitelist needs to be installed."
        )
        return

    statuses = wlm.get_all_whitelist_statuses()
    missing = [fname for fname, status in statuses.items() if not status["present"]]

    for fname, status in statuses.items():
        spec = wlm.WHITELIST_CATALOG[fname]
        chem_note = ", ".join(spec["chemistry_keys"])
        if status["present"]:
            st.markdown(f"✅ `{fname}` -- {spec['description']} ({status['size_mb']}MB) -- used by: `{chem_note}`")
        else:
            st.markdown(f"❌ `{fname}` -- {spec['description']} (~{spec['approx_decompressed_mb']}MB when downloaded) -- used by: `{chem_note}`")

    if not missing:
        st.success("✅ Every known barcode whitelist is already installed on this system.")
        return

    st.warning(f"**{len(missing)} whitelist file(s) missing.**")



    col1, col2 = st.columns(2)
    with col1:
        if st.button(f"📥 Download All Missing Whitelists ({len(missing)})", key="setup_whitelist_download_all_btn", type="primary"):
            _run_whitelist_download(missing)
    with col2:
        chosen_single = st.selectbox("Or download just one:", options=missing, key="setup_whitelist_single_select")
        if st.button("📥 Download Selected", key="setup_whitelist_download_one_btn"):
            _run_whitelist_download([chosen_single])

    _render_whitelist_install_status_panel()


# --- Minimal, self-contained progress/status tracking for whitelist
# downloads -- deliberately NOT reusing deployment_manager.py's own
# background-subprocess install-status machinery (that module's pattern
# is built around launching and polling an external `mamba`/`conda`
# subprocess with its own log file; whitelist downloads are plain
# in-process Python urllib calls with no subprocess involved at all, so
# a lighter-weight session_state-based progress tracker is a better fit
# here than forcing this into that subprocess-oriented pattern).
def _run_whitelist_download(filenames):
    progress_lines = []
    progress_area = st.empty()

    def _progress_cb(msg, _lines=progress_lines, _area=progress_area):
        _lines.append(msg)
        _area.code("\n".join(_lines), language="text")

    results = {}
    with st.spinner(f"Downloading {len(filenames)} whitelist file(s)..."):
        for fname in filenames:
            success, message = wlm.download_whitelist(fname, progress_callback=_progress_cb)
            results[fname] = (success, message)

    st.session_state["setup_whitelist_last_results"] = results
    st.rerun()


def _render_whitelist_install_status_panel():
    results = st.session_state.get("setup_whitelist_last_results")
    if not results:
        return
    st.markdown("**Last download result(s):**")
    for fname, (success, message) in results.items():
        icon = "✅" if success else "❌"
        st.markdown(f"{icon} `{fname}`: {message}")


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
        )
    elif auth_method == "password":
        st.warning("⚠️ Password auth is never saved to disk, by design.")
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
        st.info("Provide a connection name, host, and username above to save.")
        return
    if st.button("💾 Save Connection", key="setup_new_conn_save_btn", type="primary"):
        hpc.save_connection({
            "profile_name": profile_name, "host": host, "port": int(port),
            "username": username, "auth_method": auth_method, "key_path": key_path,
        })
        st.session_state["setup_conn_just_saved"] = profile_name
        for key in list(st.session_state.keys()):
            if key.startswith("setup_new_conn"):
                del st.session_state[key]
        st.rerun()


def _render_hpc_connections():
    st.subheader("🖥️ HPC Connections")
    st.markdown("Save and test SSH connections to remote HPC clusters.")

    if not auth.has_permission("manage_hpc_connections"):
        st.info(
            "🔒 Your current role does not include permission to manage HPC "
            "connections. Contact your lab's admin if you need a connection "
            "added, tested, or removed."
        )
        return

    if not hpc.PARAMIKO_AVAILABLE:
        st.error("❌ `paramiko` isn't installed in this environment.")
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
        "installed, install anything missing directly from here, set up "
        "the (optional, large, admin-only) eggNOG-mapper database, and "
        "configure SSH connections to HPC clusters."
    )
    st.markdown("---")
    _render_environment_check()
    st.markdown("---")
    _render_install_missing_section()
    st.markdown("---")
    _render_whitelist_setup()
    st.markdown("---")
    _render_eggnog_database_setup()
    st.markdown("---")
    _render_hpc_connections()
