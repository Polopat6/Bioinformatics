# Use a stable Ubuntu base image with R pre-installed
FROM rocker/r-ver:4.3.2

# 1. Install system dependencies for spatial genomics and databases
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-pip \
    python3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    procps \
    curl \
    && rm -rf /var/lib/apt/lists/*

# 2. Install R packages inside the container
RUN R -e "install.packages(c('DBI', 'duckdb', 'tidyr'), repos='https://r-project.org')"

# 3. Copy application dependency records
WORKDIR /app
COPY requirements.txt /app/requirements.txt

# 4. Install Python packages inside the container
RUN pip3 install --no-cache-dir -r requirements.txt

# 5. Copy the rest of the project files into the container
COPY . /app

# Ensure our custom R script remains executable inside Linux environments
RUN chmod +x /app/bin/process_spatial.R
ENV PATH="/app/bin:${PATH}"

# Set the default port exposure window for our user interface server
EXPOSE 8501

# Default command launches the UI if the container is run stand-alone
CMD ["streamlit", "run", "app/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
