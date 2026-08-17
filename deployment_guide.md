# Deployment Guide

This portal supports three deployment modes -- **Docker**, **local machine**,
and **HPC cluster** -- all built from the exact same
[`environment.yml`](./environment.yml), the single source of truth for every
Python package, R/Bioconductor package, and external CLI tool either the Bulk
RNA-Seq or Spatial Transcriptomics pipeline needs. `environment.yml` lives at
the **repo root** (alongside this file and the `Dockerfile`), not inside
`app/` -- it describes the whole runtime stack, not just the Streamlit app's
own Python imports:

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
    ├── ... (all other .py modules)
    └── .streamlit/config.toml
```

---

## 1. Docker

```bash
docker build -t bioportal:latest .
docker run -p 8501:8501 -v "$(pwd)/data:/workspace/repo/app/data" bioportal:latest
```

Open `http://localhost:8501`. The `-v` mount persists `data/projects/` and
`data/shared_references/` outside the container.

---

## 2. Local machine

```bash
mamba env create -f environment.yml
mamba activate bioportal
python -c "import kaleido; kaleido.get_chrome_sync()"   # one-time, see below
cd repo/app
streamlit run app.py
```

### ⚠️ One-time kaleido Chrome fetch

`python-kaleido` (used for the Differential Expression workspace's PDF/image
plot export) no longer bundles Chrome as of v1 -- it needs a compatible
Chrome/Chromium build fetched separately, **once**, after installing the
environment:

```bash
python -c "import kaleido; kaleido.get_chrome_sync()"
```

Skipping this step doesn't break the app generally -- only PDF/image export
specifically will fail the first time it's used, with an error pointing at
the missing Chrome install. The Docker image runs this automatically at
build time (see `Dockerfile`), so this step is Docker-only-unnecessary; it's
still required once for local and HPC installs.

---

## 3. HPC cluster

```bash
mamba env create -f environment.yml -p $SCRATCH/envs/bioportal
mamba activate $SCRATCH/envs/bioportal
python -c "import kaleido; kaleido.get_chrome_sync()"   # one-time, see above
```

Run inside a persistent `tmux`/`screen` session so it survives an SSH
disconnect:

```bash
tmux new -s bioportal
cd repo/app
streamlit run app.py --server.maxUploadSize=2048
# detach: Ctrl-b then d
```

Then port-forward from your local machine:

```bash
ssh -L 8501:localhost:8501 you@your-hpc-hostname
```

### ⚠️ The Streamlit config-location gotcha

Streamlit only loads `.streamlit/config.toml` **relative to the directory
`streamlit run` is invoked from.** This app is always launched via
`cd repo/app && streamlit run app.py`, so the correct config location is
**`repo/app/.streamlit/config.toml`** -- one sitting at the repo root is
silently never read. If uploads are unexpectedly capped at Streamlit's
default 200 MB, this is almost always why.

### ⚠️ The "only conda-forge got searched" gotcha (ad-hoc installs)

Running a plain `mamba install <package>` by hand against an environment
that wasn't created via `mamba env create -f environment.yml` may silently
only search whichever channels that specific environment already has
configured -- which may not include `bioconda`, where every Bioconductor R
package this project uses actually lives. **Symptom:** mamba reports a real,
correctly-named package as `does not exist (perhaps a typo or a missing
channel)`, with the install log never once mentioning `bioconda`. **Fix:**
always pass every channel explicitly:

```bash
mamba install -y -p $SCRATCH/envs/rnaseq -c conda-forge -c bioconda -c defaults <package>
```

The portal's own **"⚙️ Setup & Deployment"** page does this automatically for
every install it launches.

### What to do if a dependency audit turns up something missing

1. Use the **⚙️ Setup & Deployment** page to see exactly what's missing on
   *this* environment, right now.
2. Click **"Install Missing Dependencies"** (installs with correct channels
   automatically), or add the dependency to `environment.yml` and re-sync.

---

## Built-in environment checker & installer: the "⚙️ Setup & Deployment" page

Checks every Python package, CLI tool, and R/Bioconductor package (including
all 7 organism annotation packages) live on whatever machine it's running
on, and can install anything missing directly via a background
`mamba`/`conda install` -- with explicit `-c conda-forge -c bioconda -c
defaults` channels (avoiding the gotcha above) and an automatic per-package
fallback if a batch install fails, so one bad/incompatible spec can never
block every other package in the same batch.

That same page also lets you configure and test **SSH connections to an HPC
cluster** (SSH-key or agent-based auth -- passwords are never written to
disk) as a first step toward HPC-backed remote execution.
