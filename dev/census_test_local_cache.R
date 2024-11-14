# Description:
# This R script uses several libraries to manage, transform, and annotate a dataset of cell metadata.
# It performs the following operations:
# 1. Connects to an in-memory DuckDB database and queries parquet files containing cell metadata.
# 2. Manipulates and annotates the dataset by merging dataset IDs and cell type information.
# 3. Converts the annotated tibble into a data.table for efficient data manipulation.
# 4. Generates unique identifiers for each cell based on dataset ID and cell type, then matches these against existing .h5ad files.
# 5. Outputs the final annotated dataset with matched file identifiers to a parquet file for further analysis or visualization.
# 6. Demonstrates how to load and use metadata from a cache for efficient processing.

library(duckdb)
library(dbplyr)
library(dplyr)
library(tidyr)
library(data.table)

# Connect to DuckDB in memory, query parquet files, and perform initial data transformation
my_consensus = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/consensus_output_new_plus_non_immune_harmonisation.parquet')")
  ) |> 
  mutate(cell_ = paste(cell_, dataset_id, sep="___")) |>  # Concatenate cell ID with dataset ID
  select(cell_, dataset_id, cell_type_immune_consensus = reannotation_consensus)  # Select and rename columns for clarity

# Add a file ID column for database use and convert the modified consensus data into a tibble
my_consensus_tbl <- my_consensus |> 
  mutate(file_id_db = paste(dataset_id, cell_type_immune_consensus, sep = "___")) |> 
  as_tibble()

# Convert the tibble to a data.table for more efficient data manipulation
dt <- as.data.table(my_consensus_tbl)

# Add a new column to data.table using set() for memory efficiency
id_tbl <- set(dt, j = "file_id_cellNexus", value = sapply(dt$file_id_db, digest))

# Define the path to the directory containing .h5ad files
path = "/vast/projects/cellxgene_curated/cellxgene_temp/29_10_2024/original"

# List .h5ad files in the specified directory and create a data.table
files <- list.files(path, pattern = "\\.h5ad$")
files_dt <- data.table(value = files)

# Remove the file extension for matching purposes
files_dt <- files_dt[, value := stringr::str_replace(value, ".h5ad", "")]

# Perform an inner join between the ID table and files data table
result_dt <- merge(
  x = id_tbl, 
  y = files_dt, 
  by.x = "file_id_cellNexus", 
  by.y = "value"
)

# Add the .h5ad extension back to identifiers, convert to tibble and write to parquet
result_tbl <- result_dt |> 
  as_tibble() |> 
  mutate(file_id_cellNexus = paste0(file_id_cellNexus, ".h5ad"))
result_tbl |> 
  arrow::write_parquet("/vast/projects/cellxgene_curated/consensus_output_new_plus_non_immune_harmonisation_plus_file_id_cellNexus_for_cellxgene_temp.parquet")

# Demonstrate loading metadata from cache and retrieving specific columns
get_metadata(cache_directory = "/vast/projects/cellxgene_curated/",
             remote_url = get_database_url("metadata.1.0.0.parquet")) |> 
  select(cell_, file_id_cellNexus) |>
  head(1000) |>
  get_data_container(repository = NULL,
                     cache_directory = "/vast/projects/cellxgene_curated/cellxgene/29_10_2024/", 
                     grouping_column = "file_id_cellNexus")