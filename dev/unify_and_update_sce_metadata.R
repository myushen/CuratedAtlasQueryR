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
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_0.parquet')") )

# clean duplicated columns in cellxgene
metadata <- metadata |> select(-cell__1,
                   -dataset_id_1,
                   -dataset_id_1_1,
                   -cell__2,
                   -cell__3,
                   -dataset_id_2,
                   -dataset_id_3,
                   -cell_type_consensus_harmonised_1) |> 
  dplyr::rename(cell_annotation_blueprint_singler = blueprint_first_labels_fine,
         cell_annotation_monaco_singler = monaco_first_labels_fine,
         cell_annotation_azimuth_l2 = azimuth_predicted_celltype_l2) |>
  mutate(file_id_cellNexus = paste0(file_id_cellNexus, ".h5ad"))

# This requires a big memory machine
metadata |> as_tibble() |> arrow::write_parquet("/vast/scratch/users/shen.m/cellNexus/metadata.1.0.0.parquet")

# Repeat similar steps for the fibrosis atlas
fibrosis <- tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),  
                sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/fibrosis.0.2.3.parquet')") )
fibrosis <- fibrosis |> dplyr::rename(file_id_cellNexus = file_id_db) |> 
  mutate(file_id_cellNexus = paste0(file_id_cellNexus,".h5ad"))
fibrosis |> as_tibble() |> arrow::write_parquet("/vast/scratch/users/shen.m/cellNexus/fibrosis.1.0.0.parquet")

# Repeat similar steps for the Prostate atlas
prostate <- tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),  
                sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/ProstateAtlas2/prostate.0.1.0.parquet')") )
prostate <- prostate |> 
  mutate(file_id_cellNexus = paste0(file_id_cellNexus,".h5ad"),
         cell_type_harmonised = NA)
prostate |> as_tibble() |> arrow::write_parquet("/vast/scratch/users/shen.m/cellNexus/prostate.1.0.0.parquet")


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
