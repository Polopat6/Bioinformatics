# 🧬 Spatial Transcriptomics Atlas Explorer

A production-grade, reproducible bioinformatics engineering project showcasing cross-language pipelines, embedded relational databases, and interactive data visualization.

## 🚀 Architecture
This project intentionally leverages the unique strengths of both the **R** and **Python** ecosystems, bridging them seamlessly using an embedded analytical database without configuration friction.

[ User Interface ] (Streamlit + Plotly)│ ▲▼ │ (Optimized SQL Queries)[ Nextflow Pipeline ] ──► [ R Statistical Core ] ──► [ Built-in DuckDB ]

- **Workflow Orchestration**: **Nextflow** manages pipeline logic, error handling, and code reproducibility.
- **Spatial Statistics**: **R** performs marker selection and tissue coordinate formatting.
- **Embedded Database Layer**: **DuckDB** bridges R and Python, allowing data frames to be written as a flat file by R and read instantly via SQL in Python.
- **Data Presentation Layer**: **Python (Streamlit + Plotly)** provides an interactive visual dashboard for data exploration.

## ⚙️ How to Run the Project Locally

## ⚙️ How to Run the Project Locally

### 1. Download the 10x Genomics Dataset
Run these commands from the project root to fetch the authentic Adult Mouse Brain (FFPE) spatial layout matrices:
```bash
cd data/
curl -O https://10xgenomics.com
tar -xf Visium_FFPE_Mouse_Brain_spatial.tar.gz
cd ..
```

### 2. Execute the Pipeline
Run the Nextflow script to trigger the data processor and generate the database:
```bash
nextflow run main.nf
```

### 3. Launch the Application
Install the UI layout dependencies and fire up the web interface server:
```bash
python3 -m pip install -r requirements.txt
python3 -m streamlit run app/app.py
```
Your default web browser will automatically open to `http://localhost:8501`.