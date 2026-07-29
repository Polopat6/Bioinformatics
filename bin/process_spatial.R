#!/usr/bin/env Rscript

# 1. Load Core Libraries
library(DBI)
library(duckdb)
library(tidyr)

# Helper function to generate mock spatial data if Seurat isn't installed yet
generate_mock_spatial_data <- function() {
  message("Generating mock spatial transcriptomics dataset for development...")
  set.seed(42)
  num_spots <- 200
  
  # Simulate tissue pixel coordinates (like a Visium grid)
  spatial_coords <- data.frame(
    spot_id = paste0("spot_", 1:num_spots),
    imagecol = sample(10:80, num_spots, replace = TRUE),
    imagerow = sample(10:80, num_spots, replace = TRUE),
    seurat_clusters = as.factor(sample(0:4, num_spots, replace = TRUE))
  )
  
  # Simulate marker expression for 3 typical spatial genes
  # GFAP (Astrocytes), MS4A1 (B-Cells), MKI67 (Proliferation marker)
  genes <- c("GFAP", "MS4A1", "MKI67")
  expression_list <- list()
  
  for (gene in genes) {
    expression_list[[gene]] <- data.frame(
      gene = gene,
      spot_id = spatial_coords$spot_id,
      expression = round(runif(num_spots, min = 0, max = 5), 2)
    )
  }
  
  expression_df <- do.call(rbind, expression_list)
  return(list(metadata = spatial_coords, expression = expression_df))
}

# 2. Run Data Processing
data_layers <- generate_mock_spatial_data()

# 3. Connect to Built-in DuckDB File
message("Connecting to DuckDB database file...")
con <- dbConnect(duckdb(), dbdir = "spatial_atlas.db", read_only = FALSE)

# 4. Write Tidy R DataFrames into Database Tables
message("Writing processing layers to relational database...")
dbWriteTable(con, "spatial_metadata", data_layers$metadata, overwrite = TRUE)
dbWriteTable(con, "gene_expression", data_layers$expression, overwrite = TRUE)

# 5. Clean Close
dbDisconnect(con, shutdown = TRUE)
message("Pipeline processing complete. 'spatial_atlas.db' generated successfully!")
