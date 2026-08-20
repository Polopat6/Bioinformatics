## Deployment Guide

This portal supports three deployment modes -- **Docker**, **local machine**,
and **HPC cluster** -- all built from the exact same
[environment.yml](./environment.yml), the single source of truth for every
Python package, R/Bioconductor package, and external CLI tool the Bulk
RNA-Seq, Single-cell RNA-Seq, or Spatial Transcriptomics pipeline needs.
environment.yml lives at the **repo root** (alongside this file and the
Dockerfile), not inside app/ -- it describes the whole runtime stack, not
just the Streamlit app's own Python imports:

```
repo/
├── environment.yml
├── Dockerfile
├── DEPLOYMENT.md
├── main.nf, nextflow.config      <- Spatial Transcriptomics pipeline
├── bin/
│   └── process_spatial.R
└── app/
    ├── app.py
    ├── ... (all other Bulk RNA-Seq / Spatial / shared .py modules)
    ├── single_cell/               <- Single-cell RNA-Seq pipeline modules
    │   ├── singlecell_workspace.py
    │   ├── sc_project_manager.py
    │   ├── sc_cellqc_manager.py   (Phase 2: doublet detection, ambient
    │   │                           RNA correction, per-cell QC)
    │   ├── starsolo_manager.py
    │   └── ... (chemistry_manager.py, singlecell_ingestion_manager.py, etc.)
    └── .streamlit/config.toml
```

`app.py` adds `single_cell/` onto its own import path at startup (see that
file's own module docstring) so its modules import as plain top-level
modules, exactly like everything else in `app/` -- no separate packaging
or install step is needed for this subfolder specifically.

### 1. Docker
```
docker build -t bioportal:latest .
docker run -p 8501:8501 -v "$(pwd)/data:/workspace/repo/app/data" bioportal:latest
```
Open http://localhost:8501. The `-v` mount persists `data/projects/`,
`data/singlecell_projects/`, and `data/shared_references/` outside the
container.

### 2. Local machine
```
mamba env create -f environment.yml
mamba activate bioportal
python -c "import kaleido; kaleido.get_chrome_sync()"   # one-time, see below
Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'   # one-time, see below
cd repo/app
streamlit run app.py
```

#### ⚠️ One-time kaleido Chrome fetch
python-kaleido (used for the Differential Expression workspace's PDF/image
plot export) no longer bundles Chrome as of v1 -- it needs a compatible
Chrome/Chromium build fetched separately, **once**, after installing the
environment:
```
python -c "import kaleido; kaleido.get_chrome_sync()"
```
Skipping this step doesn't break the app generally -- only PDF/image export
specifically will fail the first time it's used, with an error pointing at
the missing Chrome install. The Docker image runs this automatically at
build time (see Dockerfile), so this step is Docker-only-unnecessary; it's
still required once for local and HPC installs.

#### ⚠️ One-time DoubletFinder install (Single-cell Cell-level QC)
DoubletFinder (an optional, selectable doublet-detection method in the
Single-cell RNA-Seq pipeline's Phase 2 Cell-level QC step -- see
`sc_cellqc_manager.py`) is **not available via conda, mamba, CRAN, or
Bioconductor at all** -- it's distributed only from its author's GitHub
repository. `environment.yml` includes `r-remotes` specifically so this
one-time manual step can be run right after environment creation:
```
Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
```
Skipping this step doesn't break the app generally -- scDblFinder (the
default doublet-detection method) works with no additional setup. Only
selecting DoubletFinder specifically as the alternate method will fail
until this one-time install is run. The Docker image runs this
automatically at build time (see Dockerfile), so this step is
Docker-only-unnecessary; it's still required once for local and HPC
installs. The portal's own **"⚙️ Setup & Deployment"** page detects whether
this specific package is installed and, if not, shows this exact command
inline (see "Built-in environment checker & installer" below) -- it cannot
run this install for you automatically the way it can for every other
dependency, since there's no conda package for it to install.

### 3. HPC cluster
```
mamba env create -f environment.yml -p $SCRATCH/envs/bioportal
mamba activate $SCRATCH/envs/bioportal
python -c "import kaleido; kaleido.get_chrome_sync()"   # one-time, see above
Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'   # one-time, see above
```
Run inside a persistent tmux/screen session so it survives an SSH
disconnect:
```
tmux new -s bioportal
cd repo/app
streamlit run app.py --server.maxUploadSize=2048
# detach: Ctrl-b then d
```
Then port-forward from your local machine:
```
ssh -L 8501:localhost:8501 you@your-hpc-hostname
```

#### ⚠️ The Streamlit config-location gotcha
Streamlit only loads `.streamlit/config.toml` **relative to the directory
`streamlit run` is invoked from.** This app is always launched via
`cd repo/app && streamlit run app.py`, so the correct config location is
**`repo/app/.streamlit/config.toml`** -- one sitting at the repo root is
silently never read. If uploads are unexpectedly capped at Streamlit's
default 200 MB, this is almost always why.

#### ⚠️ The "only conda-forge got searched" gotcha (ad-hoc installs)
Running a plain `mamba install <package>` by hand against an environment
that wasn't created via `mamba env create -f environment.yml` may silently
only search whichever channels that specific environment already has
configured -- which may not include bioconda, where every Bioconductor R
package this project uses actually lives. **Symptom:** mamba reports a real,
correctly-named package as *does not exist* (perhaps a typo or a missing
channel), with the install log never once mentioning bioconda. **Fix:**
always pass every channel explicitly:
```
mamba install -y -p $SCRATCH/envs/rnaseq -c conda-forge -c bioconda -c defaults <package>
```
The portal's own **"⚙️ Setup & Deployment"** page does this automatically for
every install it launches.

#### ⚠️ The "GitHub-only R package" gotcha
Not every R dependency this portal uses is available via conda/mamba --
DoubletFinder (see "One-time DoubletFinder install" above) is a GitHub-only
package with no conda, CRAN, or Bioconductor release at all. **Symptom:**
the "⚙️ Setup & Deployment" page's Environment Check shows it as missing,
but clicking "Install Missing Dependencies" never attempts to install it
(this is intentional, not a bug -- that button only knows how to run
conda/mamba installs, and there is no conda package for this one). **Fix:**
the Environment Check page shows the exact manual `remotes::install_github(...)`
command to run instead, right next to that package's status -- this is a
one-time step, same as kaleido's Chrome fetch above.

#### What to do if a dependency audit turns up something missing
- Use the **⚙️ Setup & Deployment** page to see exactly what's missing on
  *this* environment, right now.
- Click **"Install Missing Dependencies"** (installs with correct channels
  automatically), or add the dependency to environment.yml and re-sync.
- If the missing item is a GitHub-only R package (currently: DoubletFinder
  only), see the "GitHub-only R package" gotcha above -- it needs the
  manual command shown on that page instead.

### Built-in environment checker & installer: the "⚙️ Setup & Deployment" page
Checks every Python package, CLI tool, and R/Bioconductor package (including
all 7 organism annotation packages, and the Single-cell RNA-Seq pipeline's
own R packages -- DropletUtils, scuttle, scDblFinder, Seurat, SoupX, celda,
and DoubletFinder) live on whatever machine it's running on, and can install
anything missing directly via a background mamba/conda install -- with
explicit `-c conda-forge -c bioconda -c defaults` channels (avoiding the
channel gotcha above) and an automatic per-package fallback if a batch
install fails, so one bad/incompatible spec can never block every other
package in the same batch. GitHub-only packages (DoubletFinder) are checked
too, but shown separately with their own manual install command rather than
being offered through the same automatic installer (see the "GitHub-only R
package" gotcha above).

That same page also lets you configure and test **SSH connections to an HPC
cluster** (SSH-key or agent-based auth -- passwords are never written to
disk) as a first step toward HPC-backed remote execution.
