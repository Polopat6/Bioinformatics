import streamlit as st
import duckdb
import plotly.express as px
import os

st.set_page_config(layout="wide", page_title="Spatial Transcriptomics Atlas")
st.title("🧬 Spatial Transcriptomics Tissue Explorer")
st.markdown("---")

db_path = "spatial_atlas.db"

if not os.path.exists(db_path):
    st.error(f"⚠️ Database file `{db_path}` not found!")
else:
    con = duckdb.connect(database=db_path, read_only=True)
    st.sidebar.header("📊 Visualization Settings")
    
    genes_df = con.execute("SELECT DISTINCT gene FROM gene_expression ORDER BY gene").df()
    gene_list = genes_df['gene'].tolist()
    selected_gene = st.sidebar.selectbox("Select a marker gene to visualize:", gene_list)

    # Cluster configuration query
    cluster_df = con.execute("SELECT spot_id, imagecol, imagerow, seurat_clusters FROM spatial_metadata").df()

    # Expression layout query
    expr_query = f"""
        SELECT m.spot_id, m.imagecol, m.imagerow, e.expression
        FROM spatial_metadata m
        JOIN gene_expression e ON m.spot_id = e.spot_id
        WHERE e.gene = '{selected_gene}'
    """
    expr_df = con.execute(expr_query).df()
    con.close()

    # 📊 NEW SUMMARY METRICS ROWS START HERE 📊
    total_spots = len(cluster_df)
    total_clusters = cluster_df['seurat_clusters'].nunique()
    avg_expression = round(expr_df['expression'].mean(), 2)
    max_expression = round(expr_df['expression'].max(), 2)

    m_col1, m_col2, m_col3, m_col4 = st.columns(4)
    with m_col1:
        st.metric(label="Total Tissue Spots Indexed", value=f"{total_spots:,}")
    with m_col2:
        st.metric(label="Identified Tissue Clusters", value=total_clusters)
    with m_col3:
        st.metric(label=f"Average {selected_gene} Intensity", value=avg_expression)
    with m_col4:
        st.metric(label=f"Peak {selected_gene} Value", value=max_expression)
        
    st.markdown("---")
    # 📊 NEW SUMMARY METRICS ROWS END HERE 📊
