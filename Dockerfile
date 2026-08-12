# Use a stable Ubuntu base image with R pre-installed
FROM rocker/r-ver:4.3.2

# 1. Install system dependencies for spatial genomics, databases, and
#    bulk RNA-seq tooling (Nextflow requires Java; unzip is needed by the
#    Nextflow installer).
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-pip \
    python3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    procps \
    curl \
    default-jre-headless \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Nextflow itself (required to run main.nf inside this
#    container rather than relying on the host machine).
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow

# 3. Install fastp (adapter/quality trimming tool used by the Trimming &
#    Post-Trim QC workspace). Downloaded as a precompiled static binary
#    directly from the official fastp release site, since fastp is a
#    single compiled executable rather than a pip/conda/R package.
RUN curl -L -o /usr/local/bin/fastp http://opengene.org/fastp/fastp.0.23.4 \
    && chmod +x /usr/local/bin/fastp

# 4. Install Miniforge (conda) to get Salmon, FastQC, and MultiQC with
#    pinned, reproducible versions rather than fighting apt package
#    availability.
RUN curl -L -o /tmp/miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh
ENV PATH="/opt/conda/bin:${PATH}"

RUN conda install -y -c bioconda -c conda-forge \
    salmon=1.10.* \
    fastqc=0.12.* \
    multiqc=1.19.* \
    star=2.7.* \
    gffread=0.12.* \
    sra-tools=3.1.* \
    && conda clean -afy

# 5. Install R packages inside the container
#    - DBI/duckdb/tidyr: existing spatial pipeline dependencies (tidyr
#      kept installed just in case any script still references it, even
#      though it was found to be unused in the cleaned-up version)
#    - jsonlite: reads the job-spec JSON file passed into the DESeq2
#      Rscript (deseq2_manager.py's _DESEQ2_R_SCRIPT) -- required
#      explicitly since it is NOT a guaranteed transitive dependency of
#      DESeq2/tximport and the script hard-fails without it.
#    - tximport: imports Salmon quant.sf output into R
#    - DESeq2: differential expression analysis (Bioconductor)
#    - limma: provides removeBatchEffect(), used by the DESeq2 R script
#      to build the batch-adjusted PCA view (visualization only) shown
#      side-by-side with the raw PCA in the Differential Expression
#      workspace whenever a batch column is selected -- previously
#      missing here, which would hard-fail at runtime
#      ("there is no package called 'limma'") for any project using a
#      batch column, since the script does `library(limma)` directly.
RUN R -e "install.packages(c('DBI', 'duckdb', 'tidyr', 'jsonlite', 'BiocManager'), repos='https://cran.r-project.org')" \
    && R -e "BiocManager::install(c('tximport', 'DESeq2', 'limma'), update = FALSE, ask = FALSE)"

# 6. Copy application dependency records
WORKDIR /app
COPY requirements.txt /app/requirements.txt

# 7. Install Python packages inside the container
#    (streamlit, duckdb, plotly, pandas, openpyxl, kaleido — see
#    requirements.txt). `kaleido` is pinned/verified explicitly here
#    (in addition to whatever requirements.txt specifies) since it
#    backs the "Save as PDF" plot export feature in the Differential
#    Expression workspace, and its behavior differs meaningfully by
#    version:
#      - kaleido >= 1.0 no longer bundles Chromium; it downloads a
#        Chrome binary on first use via `kaleido.get_chrome_sync()`.
#        That download requires network access, so it's done here at
#        BUILD time (when the image build has internet access) rather
#        than left to happen at runtime, where the container may be
#        deployed without outbound network access and PDF export would
#        silently fail on its first real use.
#      - kaleido 0.2.1 (the old pin some deployments still use because
#        it bundles Chromium directly, avoiding the download above) is
#        already past Plotly's own deprecation cutoff (Sept 2025) --
#        intentionally NOT used here in favor of a current kaleido
#        release with the Chrome binary pre-fetched at build time.
RUN pip3 install --no-cache-dir -r requirements.txt \
    && pip3 install --no-cache-dir "kaleido>=1.0.0" \
    && python3 -c "import kaleido; kaleido.get_chrome_sync()"

# 8. Copy the rest of the project files into the container
COPY . /app

# Ensure our custom R script remains executable inside Linux environments
RUN chmod +x /app/bin/process_spatial.R
ENV PATH="/app/bin:${PATH}"

# Set the default port exposure window for our user interface server
EXPOSE 8501

# Default command launches the UI if the container is run stand-alone
CMD ["streamlit", "run", "app/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
