# Description:
# This script performs data management tasks involving three single-cell datasets:
# cellxgene, fibrosis, and prostate atlas. It reads metadata and dataset-specific information,
# cleans and renames columns, and writes updated data back to disk. The process involves
# connecting to databases in memory, executing SQL queries, and handling data in both
# Parquet and HDF5 formats. Additionally, it sets up directories and tests data conversion
# scripts for compatibility with Anndata structures. This is intended to unify metadata and
# streamline data handling in preparation for analysis.

library(duckdb)
library(dbplyr)
library(dplyr)
library(tidyr)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(tidySingleCellExperiment)
library(stringr)

# read cellxgene
metadata <- tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),  
    sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/Census_final_run/cell_metadata_cell_type_consensus_v1_0_6_mengyuan.parquet')") )

# Function of supporting Read parquet by duckdb, then do something, then write parquet. This avoids converting to tibble
# @example
# con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
# tbl(con, sql(paste0("SELECT * FROM read_parquet('", path_parquet, "')"))) |>
#   mutate(c = "x") |>
#   duckdb_write_parquet(path = path_parquet,
#                        con = con)
duckdb_write_parquet <- function(.tbl_sql, path, con) {
  
  sql_tbl <- 
    .tbl_sql |>
    sql_render()
  
  sql_call <- glue::glue("COPY ({sql_tbl}) TO '{path}' (FORMAT 'parquet')")
  
  res <- dbExecute(con, sql_call)
  
  return(res)
}

# clean duplicated columns in cellxgene
metadata <- metadata |> select(
                   -cell_, 
                   -cell__1,
                   -dataset_id_1,
                   -dataset_id_1_1,
                   -cell__2,
                   -cell__3,
                   -cell__4,
                   -dataset_id_2,
                   -dataset_id_3,
                   -dataset_id_4,
                   -sample_id_1,
                   -sample_id_2,
                   -sample_placeholder,
                   -cell_type_unified_ensemble_1,-cell_,
                   -cell__1,
                   -dataset_id_1,
                   -dataset_id_1_1,
                   -cell__2,
                   -cell__3,
                   -cell_type_1,
                   -dataset_id_2,
                   -observation_joinid_1,
                   -self_reported_ethnicity_1,
                   -donor_id_1,
                   -sample_id_1,
                   -assay_1,
                   -blueprint_first_labels_fine_1,
                   -monaco_first_labels_fine_1,
                   -azimuth_predicted_celltype_l2_1,
                   -cell__4,
                   -dataset_id_3,
                   -empty_droplet_1,
                   -dataset_id_4,
                   -cell_type_unified_ensemble_1,
                   -sample_id_2,
                   -atlas_id_1,
                   -tissue_1,
                   -is_primary_data_1,
                   -cell_type_ontology_term_id_1,
                   -azimuth,
                   -blueprint,
                   -monaco,
                   -matches("^scores")) |> 
  dplyr::rename(cell_annotation_blueprint_singler = blueprint_first_labels_fine,
         cell_annotation_monaco_singler = monaco_first_labels_fine,
         cell_annotation_azimuth_l2 = azimuth_predicted_celltype_l2)

metadata_path = "/vast/scratch/users/shen.m/cellNexus/metadata.1.0.5.parquet"
metadata |>  duckdb_write_parquet(path = metadata_path,
                       con = dbConnect(duckdb::duckdb(), dbdir = ":memory:"))

# Repeat similar steps for the fibrosis atlas
fibrosis <- tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),  
                sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/fibrosis.0.2.3.parquet')") )
fibrosis <- fibrosis |> dplyr::rename(file_id_cellNexus_single_cell = file_id_db) |> 
  mutate(file_id_cellNexus_single_cell = paste0(file_id_cellNexus_single_cell,".h5ad"))
#fibrosis |> as_tibble() |> arrow::write_parquet("/vast/scratch/users/shen.m/cellNexus/fibrosis.1.0.4.parquet")
fibrosis_path = "/vast/scratch/users/shen.m/cellNexus/fibrosis.1.0.4.parquet"
fibrosis |>  duckdb_write_parquet(path = fibrosis_path,
                                  con = con)

# Repeat similar steps for the Prostate atlas
prostate <- tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),  
                sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/ProstateAtlas/prostate.0.1.0.parquet')") )
prostate <- prostate |> 
  mutate(file_id_cellNexus_single_cell = paste0(file_id_db,".h5ad"),
         cell_type_harmonised = NA)
#prostate |> as_tibble() |> arrow::write_parquet("/vast/scratch/users/shen.m/cellNexus/prostate.1.0.4.parquet")
prostate_path = "/vast/scratch/users/shen.m/cellNexus/prostate.1.0.4.parquet"
prostate |>  duckdb_write_parquet(path = prostate_path,
                                  con = con)

# Convert counts in HDF5 SCE to Anndata in atlas and rename gene symbols to ensembl IDs. This is done by Target parallelation:
# script: ~/git_control/CuratedAtlasQueryR/dev/convert_fibrosis_and_prostate_hdf5_to_anndata_targets.R
# store: scratch/cellNexus
original_hdf5_download_path = "/vast/scratch/users/shen.m/cellNexus/original_hdf5/"
cpm_hdf5_download_path = "/vast/scratch/users/shen.m/cellNexus/cpm_hdf5/"
if (!dir.exists(original_hdf5_download_path))  dir.create(original_hdf5_download_path, recursive = TRUE)
if (!dir.exists(cpm_hdf5_download_path))  dir.create(cpm_hdf5_download_path, recursive = TRUE)

# test fibrosis and prostate with anndata 
cache = "/vast/scratch/users/shen.m/cellNexus"
get_metadata(cache_directory = cache,
             get_database_url(databases = NULL)) |> 
  dplyr::filter(file_id_cellNexus %in% c( "0004e421765504041c8a460a83de2d01.h5ad", # this is from cellxgene
                                          "0a54b616d9afd26c9da310e8a504b541.h5ad", # this is from prostate atlas
                                          "12eb5fe25994253c1d320ca590a6e681.h5ad"  # this is from fibrosis atlas
                                          )
               ) |>
  cellNexus:::get_data_container(repository = NULL,
                                 cache_directory = cache,
                                 grouping_column = "file_id_cellNexus",
                                 #assays = "cpm"
                                 )
