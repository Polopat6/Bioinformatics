# 🧬 Pretty Awesome Transcriptomics Tool

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

A production-grade, multi-pipeline bioinformatics platform for **Bulk RNA-Seq** and **Single-cell RNA-Seq** analysis — from raw FASTQ/SRA ingestion through alignment, quality control, differential expression, and functional enrichment — built as a single, cohesive Streamlit application designed to run identically on a laptop, an HPC cluster, or in Docker.

A bonus **Spatial Transcriptomics Atlas Viewer** module is also included as a stretch-goal proof of concept, built on a separate R/Nextflow/DuckDB architecture (see below).

---

## 🚀 What's Inside

This portal implements two complete, end-to-end RNA-seq pipelines behind a single web interface:

### 📈 Bulk RNA-Seq Pipeline

`FASTQ / SRA / BAM ingestion → FastQC + MultiQC → fastp trimming → STAR or Salmon alignment/quantification → DESeq2 differential expression → GO / KEGG / Reactome ontology analysis`

- Flexible sample ingestion: browser upload, server-side directory browse, direct SRA/NCBI fetch, or BAM→FASTQ recovery (`samtools fastq`)
- STAR (genome alignment) or Salmon (pseudo-alignment) quantification paths
- DESeq2 analysis supporting multivariate designs, batch covariates, interaction terms, Wald or LRT tests, and multiple simultaneous contrasts — with model-fit diagnostics (dispersion estimates, size factors, sample-distance QC, Cook's-distance/independent-filtering explanations, p-value histogram and MA-plot bias checks)
- Full ontology analysis suite (ORA, GSEA, compareCluster) across GO/KEGG/Reactome with GOSemSim-based term simplification

### 🧫 Single-cell RNA-Seq Pipeline

`FASTQ / SRA / BAM ingestion + chemistry detection → trimming/QC → STARsolo alignment + cell calling → Cell-level QC → Downstream Analysis (Phase 3)`

- Supports 10x Genomics (3′/5′, all chemistry versions), Drop-seq, inDrops, and BD Rhapsody
- Original-format BAM recovery (10x `bamtofastq`) for legacy/public depositions missing barcode+UMI data in their standard SRA FASTQ export — including automatic resolution of v1-chemistry's 4-file-per-lane barcode/UMI split
- Cell-level QC: doublet detection (scDblFinder), ambient RNA correction (DecontX/SoupX), adaptive per-cell filtering
- Full multi-sample downstream workflow: normalization (shifted-log or Pearson residuals) → HVG selection → PCA → Harmony batch correction → Leiden/Louvain clustering → UMAP/t-SNE (2D & 3D) → cell-type annotation (manual marker scoring + CellTypist) → **pseudobulk aggregation feeding directly into the Bulk pipeline's own DESeq2 workflow** → propeller/compositional (cell-type proportion) analysis

Both pipelines share common infrastructure for reference genome management, SRA ingestion, and gene ID/symbol resolution.

### 🤖 Advanced Mode & Monitor Mode (Bulk RNA-Seq)

For unattended, hands-off operation, the Bulk RNA-Seq pipeline supports two headless execution modes, both built on the same resumable, detached-background-process runner (FASTQ → QC → trimming → reference → quantification → counts matrix, then stopping — DESeq2/Ontology remain deliberate, separate interactive steps):

- **Advanced Mode** — define a project upfront (samples, metadata, genome, options) and launch it as a single unattended background run.
- **Monitor Mode** — run any number of independent watch-folder daemons simultaneously, each with its own configuration, that automatically validate (folder stability, metadata presence, sample/metadata count reconciliation) and launch a new Advanced-Mode-style run the moment a complete, well-formed sample folder appears.

Both modes support optional email and free-carrier-gateway text notifications (via a configurable SMTP relay) on run success/failure. Every completed run also produces a shareable **QC certificate** (JSON + standalone HTML) summarizing every stage's quality checks — readable by a collaborator without opening MultiQC. See [**advanced_monitor.md**](advanced_monitor.md) for full configuration details.

### 🗺️ Bonus Module: Spatial Transcriptomics Atlas Viewer

A separate proof-of-concept pipeline for 10x Visium spatial data, intentionally built on a different architecture to demonstrate cross-language pipeline design:

```
[ Nextflow ] ──► [ R: marker selection + coordinate formatting ] ──► [ DuckDB ] ──► [ Streamlit + Plotly ]
```

R performs spatial statistics and writes results to disk; DuckDB serves as an embedded, config-free bridge letting Python query R's output directly via SQL; Streamlit/Plotly renders the interactive tissue atlas viewer.

---

## 🏗️ Architecture

```
                     ┌─────────────────────────────┐
                     │   Streamlit UI (app.py)     │
                     │  Bulk RNA-Seq | Single-cell │
                     └──────────────┬──────────────┘
                                    │
        ┌───────────────────────────┴───────────────────────────┐
        │                                                       │
┌───────▼──────────┐                                  ┌─────────▼──────────┐
│ Python managers  │  subprocess calls to external    │   R (via Rscript)  │
│ (scanpy/anndata, │  tools                           │ DESeq2, propeller, │
│  pandas, etc.)   │ ────────────────────────────────►│ clusterProfiler,   │
│                  │  FastQC/MultiQC, fastp,          │ ReactomePA, limma  │
│                  │  STAR/Salmon, STARsolo,          │                    │
│                  │  SRA Toolkit, samtools,          │                    │
│                  │  bamtofastq                      │                    │
└──────────────────┘                                  └────────────────────┘
```

Every external tool (FastQC, fastp, STAR, Salmon, STARsolo, SRA Toolkit, samtools, R/Bioconductor) is invoked as a managed subprocess with structured logging, job-spec JSON hand-off (for R), and explicit error surfacing rather than silent failure.

**Two cross-cutting modules underpin everything on disk:**

- **`app_paths.py`** — every storage location is anchored to the `app/` directory via `os.path.abspath(__file__)`, so data paths never depend on which directory the app was launched from.
- **`atomic_io.py`** — all persistent JSON state (project tracking, user accounts, monitor configs, analysis recipes) is written temp-file-then-`os.replace()`, so an interrupted write can never leave a truncated, unreadable file behind.

---

## 📂 Repository Structure

```
Bioinformatics/
├── environment.yml              # Single source of truth for every Python, R/Bioconductor,
│                                #   and CLI dependency (Docker, local, and HPC all build from this)
├── Dockerfile                   # Container build (all pipelines)
├── docker-compose.yml           # Preferred way to run on a persistent server
├── .dockerignore                # Keeps data/ and credentials out of the build context
├── start_portal.sh              # One-command tmux-wrapped local/HPC startup
├── ruff.toml                    # Linter config (F821/F811 — undefined names, duplicate defs)
├── deployment_guide.md          # Full Docker / local / HPC deployment instructions
├── advanced_monitor.md          # Advanced Mode / Monitor Mode design + configuration
├── main.nf, nextflow.config     # Spatial Transcriptomics pipeline orchestration
├── bin/
│   └── process_spatial.R        # Spatial module's R statistical core
└── app/
    ├── app.py                            # Main Streamlit entry point / routing
    ├── app_paths.py                      # Launch-independent data path anchoring
    ├── atomic_io.py                      # Crash-safe JSON read/write + cross-process locking
    ├── auth_manager.py                   # Login, roles, permissions (PBKDF2-HMAC-SHA256)
    ├── project_manager.py                # Bulk project creation, step tracking, shared resources
    ├── project_actions.py                # Project export / archive / delete
    ├── setup_workspace.py                # ⚙️ Setup & Deployment page (env check + installer)
    ├── deployment_manager.py             # Background conda/mamba install runner
    ├── hpc_manager.py                    # SSH connection profiles (paramiko)
    ├── whitelist_manager.py              # Shared 10x barcode whitelist downloads
    ├── bulk_rnaseq_workspace.py          # Bulk pipeline: ingestion UI
    ├── bulk_bam_manager.py               # Bulk BAM→FASTQ recovery (samtools)
    ├── trimming_workspace.py             # fastp trimming + post-trim QC UI
    ├── alignment_workspace.py            # STAR/Salmon reference + alignment UI
    ├── differential_expression_workspace.py   # DESeq2 UI
    ├── deseq2_manager.py                 # DESeq2 R subprocess backend
    ├── ontology_workspace.py / ontology_manager.py   # GO/KEGG/Reactome UI + backend
    ├── eggnog_manager.py                 # Orthology-based annotation (non-model organisms)
    ├── advanced_mode_orchestrator.py / advanced_mode_workspace.py   # Advanced Mode
    ├── monitor_manager.py / monitor_mode_workspace.py / monitor_daemon_runner.py  # Monitor Mode
    ├── qc_certificate_manager.py         # Shareable run QC certificate (JSON + HTML)
    ├── notification_manager.py           # Email / carrier-gateway text notifications
    ├── reference_manager.py, sra_manager.py, fastqc_manager.py, fastp_manager.py,
    │   quantification_manager.py, counts_matrix_manager.py, gene_id_mapper.py,
    │   ingestion_manager.py, file_browser.py
    ├── spatial_workspace.py              # Spatial module's Streamlit UI
    ├── .streamlit/config.toml            # Streamlit server config (upload size, etc.)
    └── single_cell/                      # All Single-cell RNA-Seq pipeline modules
        ├── singlecell_workspace.py       # Phase 1 (ingestion→alignment) + Phase 2 (Cell QC) UI
        ├── sc_downstream_workspace.py    # Phase 3 downstream analysis UI
        ├── sc_deseq2_workspace.py        # Pseudobulk → DESeq2 bridge UI
        ├── sc_ontology_workspace.py / sc_comparecluster_workspace.py
        ├── sc_downstream_manager.py, sc_cellqc_manager.py, sc_project_manager.py
        ├── starsolo_manager.py, sc_sra_manager.py, sc_genetic_demux_manager.py
        ├── sc_marker_panel_io.py, sc_marker_confidence.py, marker_library_manager.py
        └── chemistry_manager.py, gff3_gene_name_resolver.py,
            singlecell_ingestion_manager.py, singlecell_trim_manager.py
```

`app.py` dynamically adds `single_cell/` to its own import path at startup, so every module in that subfolder imports as a plain top-level module — no separate packaging or install step needed for it.

---

## ⚙️ How to Run the Project

This project runs consistently across **Docker, local machines, and HPC clusters**, using the single shared [`environment.yml`](./environment.yml) at the repo root as the source of truth for every Python package, R/Bioconductor package, and external CLI tool all three pipelines need. See [**deployment_guide.md**](deployment_guide.md) for the full walkthrough (including HPC-specific tmux/port-forwarding instructions and troubleshooting), and [**advanced_monitor.md**](advanced_monitor.md) for unattended Advanced/Monitor Mode setup.

### Docker (recommended — fully self-contained)

```bash
# One-time: match the container's user to yours so it can write to ./data
echo "UID=$(id -u)" >  .env
echo "GID=$(id -g)" >> .env

docker compose up -d --build
docker compose logs -f
```

Then open **http://localhost:8501** (or `http://<server-ip>:8501` on a home/lab server).

The `./data` bind mount persists `projects/`, `singlecell_projects/`, `shared_references/`, `shared_whitelists/`, `marker_libraries/`, `auth/`, and `monitor/` outside the container. Verify it's wired correctly before doing real work:

```bash
docker compose exec bioportal python -c "import app_paths; print(app_paths.describe_data_root())"
```

`data_root` should read `/workspace/repo/app/data`. Create a test project, then confirm it appears on the host with `ls data/projects/`.

> **⚠️ Don't expose port 8501 directly to the internet.** The built-in login is a deliberate-but-basic gate — there's no rate limiting, account lockout, or session expiry. On a LAN behind a router it's fine. For remote access, bind to localhost only (`"127.0.0.1:8501:8501"` in `docker-compose.yml`) and use an SSH tunnel:
> ```bash
> ssh -L 8501:localhost:8501 you@your-server
> ```

The first build takes roughly 20–40 minutes (the conda environment layer dominates) and is cached on every rebuild that doesn't change `environment.yml`.

### Local machine

```bash
mamba env create -f environment.yml
mamba activate bioportal
./start_portal.sh
```

`start_portal.sh` auto-detects your conda installation, wraps itself in a persistent tmux session, and launches Streamlit from the correct directory. Override anything without editing the file:

```bash
BIOPORTAL_ENV=my_env BIOPORTAL_PORT=9000 ./start_portal.sh
./start_portal.sh --attach     # reattach to a running portal
./start_portal.sh --no-tmux    # foreground (systemd, debugging)
```

Or launch manually:

```bash
cd app && streamlit run app.py
```

⚠️ **Two one-time manual steps are required after environment creation** (the Docker image performs both automatically at build time — this applies only to local/HPC installs):

```bash
# PDF/image plot export (Differential Expression workspace)
python -c "import kaleido; kaleido.get_chrome_sync()"

# DoubletFinder (optional Single-cell doublet-detection method — GitHub-only, no conda/CRAN release)
Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
```

Skipping either doesn't break the app generally — only that specific feature fails on first use, with a clear inline error. scDblFinder (the default doublet method) needs no extra setup.

### HPC cluster

```bash
mamba env create -f environment.yml -p $SCRATCH/envs/bioportal
mamba activate $SCRATCH/envs/bioportal
# (run the same two one-time steps shown above)

BIOPORTAL_CONDA_ROOT=$SCRATCH/envs ./start_portal.sh
# detach from tmux: Ctrl-b then d
```

Then port-forward from your local machine:

```bash
ssh -L 8501:localhost:8501 you@your-hpc-hostname
```

> **⚠️ The Streamlit config-location gotcha.** Streamlit only reads `.streamlit/config.toml` relative to the directory `streamlit run` is invoked *from*. This app is always launched from `app/`, so the config lives at **`app/.streamlit/config.toml`** — one at the repo root is silently never read. If uploads are unexpectedly capped at Streamlit's 200 MB default, this is almost always why. `start_portal.sh` and the Dockerfile both also pass `--server.maxUploadSize` explicitly as belt-and-suspenders.

### 🩺 Built-in environment checker: the "⚙️ Setup & Deployment" page

Once running, the portal's own **Setup & Deployment** page audits every Python, R/Bioconductor, and CLI dependency live on whatever machine it's running on, and can install anything missing via a background conda/mamba install — with `-c conda-forge -c bioconda -c defaults` applied automatically, and a per-package fallback so one bad spec can't block an entire batch. GitHub-only packages without a conda equivalent (currently just DoubletFinder) are flagged separately with the exact manual command to run. The same page also configures and tests SSH connections to an HPC cluster.

---

## 🧰 Tech Stack

- **UI**: Streamlit, Plotly
- **Single-cell analysis**: scanpy, anndata, harmonypy, CellTypist
- **Single-cell QC (R)**: DropletUtils, scuttle, scDblFinder, Seurat, SoupX, celda, DoubletFinder (GitHub-only)
- **Bulk & single-cell alignment**: STAR / STARsolo, Salmon, fastp, FastQC, MultiQC, SRA Toolkit, samtools, 10x `bamtofastq`
- **Statistics (R, via Rscript subprocess)**: DESeq2, limma, speckle (propeller), clusterProfiler, ReactomePA, GOSemSim
- **Spatial module**: Nextflow, R, DuckDB
- **Notifications**: smtplib (email + free carrier email-to-SMS gateways)
- **Deployment**: Docker / Docker Compose, conda/mamba (single shared `environment.yml`), HPC-compatible via tmux session management

---

## 🤖 A Note on AI-Assisted Development

This project was built with heavy use of AI coding assistance (Claude) throughout the design, implementation, debugging, and documentation process, used as an active engineering collaborator rather than a one-off code generator — planning architecture, writing and reviewing implementation code, diagnosing real bugs against actual runtime output and tracebacks, and validating fixes against real data before considering an issue closed.

Every feature in this repository was still driven by explicit domain requirements, tested against real datasets (e.g. Kang et al. 2018 GSE96583, Thompson et al. 2021 GSE166992), and iterated on based on genuine execution failures encountered on real HPC/Docker/local deployments — not accepted as-is from a single generation pass. I take full ownership of the architecture, the correctness of the science, and every design decision in this codebase.

I'm sharing this openly because I believe effective, transparent use of AI tooling is quickly becoming a core professional skill in software and computational biology alike, not something to obscure — and because the resulting engineering discipline (systematic debugging, real-data validation, honest documentation of limitations) is itself representative of how I approach this work.

---

## 🗺️ Roadmap

- Non-model organism annotation support via eggNOG-mapper (pending a stable non-beta release with working database downloads)
- Sample demultiplexing (genetic and barcode-based) — currently an intentional, documented limitation
- Reusable cross-project marker libraries (`marker_library_manager.py` — backend implemented, UI not yet wired)
- Expanded single-cell technology support (10x Flex, CITE-seq, Parse Evercode, Smart-seq)
- Queryable multi-omics data warehouse (TileDB-SOMA) spanning bulk, single-cell, and spatial pipeline runs
- Deeper feature-branch development of the Spatial Transcriptomics module

---

## 👤 Author

Built by [Patrick Clouser](https://github.com/Polopat6) as a demonstration of end-to-end bioinformatics pipeline engineering across bulk and single-cell transcriptomics.

## 📄 License

This project is licensed under the [MIT License](./LICENSE).
