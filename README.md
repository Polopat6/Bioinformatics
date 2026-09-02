# 🧬 Bioinformatics Portal

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

A production-grade, multi-pipeline bioinformatics platform for **Bulk RNA-Seq** and **Single-cell RNA-Seq** analysis — from raw FASTQ/SRA ingestion through alignment, quality control, differential expression, and functional enrichment — built as a single, cohesive Streamlit application designed to run identically on a laptop, an HPC cluster, or in Docker.

A bonus **Spatial Transcriptomics Atlas Viewer** module is also included as a stretch-goal proof of concept, built on a separate R/Nextflow/DuckDB architecture (see below).

---

## 🚀 What's Inside

This portal implements two complete, end-to-end RNA-seq pipelines behind a single web interface:

### 📈 Bulk RNA-Seq Pipeline
`FASTQ / SRA / BAM ingestion → FastQC + MultiQC → fastp trimming → STAR or Salmon alignment/quantification → DESeq2 differential expression → GO / KEGG / Reactome ontology analysis`

- Flexible sample ingestion: browser upload, direct SRA/NCBI fetch, or BAM→FASTQ recovery
- STAR (genome alignment) or Salmon (pseudo-alignment) quantification paths
- DESeq2 analysis supporting multivariate designs, batch covariates, interaction terms, Wald or LRT tests, and multiple simultaneous contrasts — with model-fit diagnostics (dispersion estimates, size factors, sample-distance QC, Cook's-distance/independent-filtering explanations, p-value histogram and MA-plot bias checks)
- Full ontology analysis suite (ORA, GSEA, compareCluster) across GO/KEGG/Reactome with GOSemSim-based term simplification

### 🧫 Single-cell RNA-Seq Pipeline
`FASTQ / SRA / BAM ingestion + chemistry detection → trimming/QC → STARsolo alignment + cell calling → Cell-level QC → Downstream Analysis (Phase 3)`

- Supports 10x Genomics (3′/5′, all chemistry versions), Drop-seq, inDrops, and BD Rhapsody
- Original-format BAM recovery (10x `bamtofastq`) for legacy/public depositions missing barcode+UMI data in their standard SRA FASTQ export
- Cell-level QC: doublet detection (scDblFinder), ambient RNA correction (DecontX/SoupX), adaptive per-cell filtering
- Full multi-sample downstream workflow: normalization (shifted-log or Pearson residuals) → HVG selection → PCA → Harmony batch correction → Leiden/Louvain clustering → UMAP/t-SNE (2D & 3D) → cell-type annotation (manual marker scoring + CellTypist) → **pseudobulk aggregation feeding directly into the Bulk pipeline's own DESeq2 workflow** → propeller/compositional (cell-type proportion) analysis

Both pipelines share common infrastructure for reference genome management, SRA ingestion, and gene ID/symbol resolution.

### 🤖 Advanced Mode & Monitor Mode (Bulk RNA-Seq)
For unattended, hands-off operation, the Bulk RNA-Seq pipeline supports two headless execution modes, both built on the same resumable, detached-background-process runner (FASTQ → QC → trimming → reference → quantification → counts matrix, then stopping — DESeq2/Ontology remain deliberate, separate interactive steps):

- **Advanced Mode** — define a project upfront (samples, metadata, genome, options) and launch it as a single unattended background run.
- **Monitor Mode** — run any number of independent watch-folder daemons simultaneously, each with its own configuration, that automatically validate (folder stability, metadata presence, sample/metadata count reconciliation) and launch a new Advanced-Mode-style run the moment a complete, well-formed sample folder appears.

Both modes support optional email and free-carrier-gateway text notifications (via a configurable SMTP relay) on run success/failure. See **[advanced_monitor.md](advanced_monitor.md)** for full configuration details.

### 🗺️ Bonus Module: Spatial Transcriptomics Atlas Viewer
A separate proof-of-concept pipeline for 10x Visium spatial data, intentionally built on a different architecture to demonstrate cross-language pipeline design:

```
[ Nextflow ] ──► [ R: marker selection + coordinate formatting ] ──► [ DuckDB ] ──► [ Streamlit + Plotly ]
```

R performs spatial statistics and writes results to disk; DuckDB serves as an embedded, config-free bridge letting Python query R's output directly via SQL; Streamlit/Plotly renders the interactive tissue atlas viewer.

---

## 🏗️ Architecture

```
                     ┌────────────────────────────┐
                     │   Streamlit UI (app.py)     │
                     │  Bulk RNA-Seq | Single-cell │
                     └──────────────┬──────────────┘
                                    │
        ┌───────────────────────────┴───────────────────────────┐
        │                                                        │
┌───────▼────────┐                                     ┌─────────▼─────────┐
│ Python managers │  subprocess calls to external tools │   R (via Rscript)  │
│ (scanpy/anndata, │ ───────────────────────────────────►│ DESeq2, propeller, │
│  pandas, etc.)   │  FastQC/MultiQC, fastp, STAR/Salmon, │ clusterProfiler,   │
│                  │  STARsolo, SRA Toolkit, bamtofastq   │ ReactomePA, limma  │
└──────────────────┘                                      └────────────────────┘
```

Every external tool (FastQC, fastp, STAR, Salmon, STARsolo, SRA Toolkit, R/Bioconductor) is invoked as a managed subprocess with structured logging, job-spec JSON hand-off (for R), and explicit error surfacing rather than silent failure.

---

## 📂 Repository Structure

```
Bioinformatics/
├── environment.yml              # Single source of truth for every Python, R/Bioconductor,
│                                 #   and CLI dependency (Docker, local, and HPC all build from this)
├── Dockerfile                   # Container build (all pipelines)
├── deployment_guide.md          # Full Docker / local / HPC deployment instructions
├── advanced_monitor.md          # Advanced Mode / Monitor Mode design + configuration
├── main.nf, nextflow.config     # Spatial Transcriptomics pipeline orchestration
├── bin/
│   └── process_spatial.R        # Spatial module's R statistical core
└── app/
    ├── app.py                            # Main Streamlit entry point / routing
    ├── bulk_rnaseq_workspace.py          # Bulk pipeline: ingestion UI
    ├── alignment_workspace.py            # STAR/Salmon reference + alignment UI
    ├── differential_expression_workspace.py   # DESeq2 UI
    ├── deseq2_manager.py                 # DESeq2 R subprocess backend
    ├── ontology_workspace.py / ontology_manager.py   # GO/KEGG/Reactome UI + backend
    ├── advanced_mode_orchestrator.py / advanced_mode_workspace.py   # Advanced Mode
    ├── monitor_manager.py / monitor_mode_workspace.py / monitor_daemon_runner.py  # Monitor Mode
    ├── notification_manager.py           # Email / carrier-gateway text notifications
    ├── reference_manager.py, sra_manager.py, fastqc_manager.py, fastp_manager.py, ...
    ├── spatial_workspace.py              # Spatial module's Streamlit UI
    ├── .streamlit/config.toml            # Streamlit server config (upload size, etc.)
    └── single_cell/                      # All Single-cell RNA-Seq pipeline modules
        ├── singlecell_workspace.py       # Phase 1 (ingestion→alignment) + Phase 2 (Cell QC) UI
        ├── sc_downstream_workspace.py    # Phase 3 downstream analysis UI
        ├── sc_downstream_manager.py, sc_cellqc_manager.py, sc_project_manager.py
        ├── starsolo_manager.py, sc_sra_manager.py, sc_genetic_demux_manager.py
        └── chemistry_manager.py, gff3_gene_name_resolver.py, ...
```

`app.py` dynamically adds `single_cell/` to its own import path at startup, so every module in that subfolder imports as a plain top-level module — no separate packaging or install step needed for it.

---

## ⚙️ How to Run the Project

This project is designed to run consistently across **Docker, local machines, and HPC clusters**, using the single shared [`environment.yml`](./environment.yml) at the repo root as the source of truth for every Python package, R/Bioconductor package, and external CLI tool all three pipelines need. See **[deployment_guide.md](deployment_guide.md)** for the full walkthrough (including HPC-specific tmux/port-forwarding instructions and troubleshooting), and **[advanced_monitor.md](advanced_monitor.md)** for unattended Advanced/Monitor Mode setup.

### Docker (recommended — fully self-contained)
```bash
docker build -t bioportal:latest .
docker run -p 8501:8501 -v "$(pwd)/data:/workspace/repo/app/data" bioportal:latest
```
The `-v` mount persists `data/projects/`, `data/singlecell_projects/`, and `data/shared_references/` outside the container. Open **http://localhost:8501**.

### Local machine
```bash
mamba env create -f environment.yml
mamba activate bioportal
cd app
streamlit run app.py
```

> ⚠️ **Two one-time manual steps are required after environment creation** (Docker performs both automatically at build time — this only applies to local/HPC installs):
> ```bash
> # PDF/image plot export (Differential Expression workspace)
> python -c "import kaleido; kaleido.get_chrome_sync()"
> # DoubletFinder (optional Single-cell doublet-detection method — GitHub-only, no conda/CRAN release)
> Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
> ```
> Skipping either step doesn't break the app generally — only that specific feature will fail on first use, with a clear inline error. See **[deployment_guide.md](deployment_guide.md)** for details.

### HPC cluster
```bash
mamba env create -f environment.yml -p $SCRATCH/envs/bioportal
mamba activate $SCRATCH/envs/bioportal
# (run the same two one-time steps shown above)
tmux new -s bioportal
cd app
streamlit run app.py --server.maxUploadSize=2048
# detach: Ctrl-b then d, then port-forward from your local machine:
# ssh -L 8501:localhost:8501 you@your-hpc-hostname
```

### 🩺 Built-in environment checker: the "⚙️ Setup & Deployment" page
Once running, the portal's own **Setup & Deployment** page audits every Python, R/Bioconductor, and CLI dependency live on whatever machine it's running on, and can install anything missing directly via a background conda/mamba install (with correct `-c conda-forge -c bioconda -c defaults` channels applied automatically). GitHub-only packages without a conda equivalent (currently just DoubletFinder) are flagged separately with the exact manual command to run. The same page also supports configuring and testing SSH connections to an HPC cluster for remote execution.

---

## 🧰 Tech Stack

- **UI**: Streamlit, Plotly
- **Single-cell analysis**: scanpy, anndata, harmonypy, CellTypist
- **Single-cell QC (R)**: DropletUtils, scuttle, scDblFinder, Seurat, SoupX, celda, DoubletFinder (GitHub-only)
- **Bulk & single-cell alignment**: STAR / STARsolo, Salmon, fastp, FastQC, MultiQC, SRA Toolkit, 10x `bamtofastq`
- **Statistics (R, via Rscript subprocess)**: DESeq2, limma, speckle (propeller), clusterProfiler, ReactomePA, GOSemSim
- **Spatial module**: Nextflow, R, DuckDB
- **Notifications**: smtplib (email + free carrier email-to-SMS gateways)
- **Deployment**: Docker, conda/mamba (single shared `environment.yml`), HPC-compatible via tmux session management

---

## 🤖 A Note on AI-Assisted Development

This project was built with heavy use of AI coding assistance (Claude) throughout the design, implementation, debugging, and documentation process, used as an active engineering collaborator rather than a one-off code generator — planning architecture, writing and reviewing implementation code, diagnosing real bugs against actual runtime output and tracebacks, and validating fixes against real data before considering an issue closed.

Every feature in this repository was still driven by explicit domain requirements, tested against real datasets (e.g. Kang et al. 2018 GSE96583, Thompson et al. 2021 GSE166992), and iterated on based on genuine execution failures encountered on real HPC/Docker/local deployments — not accepted as-is from a single generation pass. I take full ownership of the architecture, the correctness of the science, and every design decision in this codebase.

I'm sharing this openly because I believe effective, transparent use of AI tooling is quickly becoming a core professional skill in software and computational biology alike, not something to obscure — and because the resulting engineering discipline (systematic debugging, real-data validation, honest documentation of limitations) is itself representative of how I approach this work.

---

## 🗺️ Roadmap

- Non-model organism annotation support via eggNOG-mapper (pending a stable non-beta release with working database downloads)
- Sample demultiplexing (genetic and barcode-based) — currently an intentional, documented limitation
- Expanded single-cell technology support (10x Flex, CITE-seq, Parse Evercode, Smart-seq)
- Queryable multi-omics data warehouse (TileDB-SOMA) spanning bulk, single-cell, and spatial pipeline runs
- Deeper feature-branch development of the Spatial Transcriptomics module

---

## 👤 Author

Built by [Patrick Clouser](https://github.com/Polopat6) as a portfolio project demonstrating end-to-end bioinformatics pipeline engineering across bulk and single-cell transcriptomics.

## 📄 License

This project is licensed under the [MIT License](./LICENSE).
