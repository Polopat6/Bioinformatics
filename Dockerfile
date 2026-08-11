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
#    - tximport: imports Salmon quant.sf output into R (future step)
#    - DESeq2: differential expression analysis (future step, from
#      Bioconductor)
RUN R -e "install.packages(c('DBI', 'duckdb', 'tidyr', 'BiocManager'), repos='https://r-project.org')" \
    && R -e "BiocManager::install(c('tximport', 'DESeq2'), update = FALSE, ask = FALSE)"

# 6. Copy application dependency records
WORKDIR /app
COPY requirements.txt /app/requirements.txt

# 7. Install Python packages inside the container
#    (streamlit, duckdb, plotly, pandas, openpyxl — see requirements.txt)
RUN pip3 install --no-cache-dir -r requirements.txt

# 8. Copy the rest of the project files into the container
COPY . /app

# Ensure our custom R script remains executable inside Linux environments
RUN chmod +x /app/bin/process_spatial.R
ENV PATH="/app/bin:${PATH}"

# Set the default port exposure window for our user interface server
EXPOSE 8501

# Default command launches the UI if the container is run stand-alone
CMD ["streamlit", "run", "app/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
