library(glue)
library(targets)
library(HDF5Array)
library(tibble)
library(dplyr)
library("AnnotationDbi")
library("org.Hs.eg.db")

# Define the path for the anndata directory (run original and cpm separatly)
input_directory = "/vast/scratch/users/shen.m/cellNexus/cpm_hdf5"
output_directory = "/vast/scratch/users/shen.m/cellNexus/cpm_anndata"

# test-pipeline
# input_directory = "/vast/scratch/users/shen.m/Census_rerun/cellNexus/original/"
# output_directory = "/vast/scratch/users/shen.m/Census_rerun/cellNexus/anndata/"

# Create the output directory if it does not exist
if (!dir.exists(output_directory))  dir.create(output_directory, recursive = TRUE)

computing_resources = crew.cluster::crew_controller_slurm(
  slurm_memory_gigabytes_per_cpu = 20, 
  slurm_cpus_per_task = 1,
  workers = 100,
  verbose = TRUE
)
tar_option_set(
  memory = "transient",
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker",
  format = "qs",
  controller = computing_resources
)

list(
  tar_target(
    conversion_tbl,
    tibble(
      hdf5 = list.files(path = glue("{input_directory}/"), full.names = TRUE),
      anndata = glue("{output_directory}/{basename(hdf5)}.h5ad") |> as.character()
    )),
  tar_target(
    sce_hdf5,
    loadHDF5SummarizedExperiment(conversion_tbl$hdf5),
    pattern = map(conversion_tbl),
    iteration = "list"
  ),
  tar_target(
    convert_gene_to_ensembl_hdf5,
    {
      # Convert gene names to ensembl ID
      # Note that gene annotations is constantly changing, 
      # none of org.Hs.eg.db, EnsDb.Hsapiens.v86, biomaRt 100% map gene names to ensembl, choose org.Hs.eg.db here
      ensembl_ids = mapIds(org.Hs.eg.db,
                           keys=sce_hdf5 |> rownames(), 
                           column="ENSEMBL",
                           keytype="SYMBOL",
                           multiVals="first")
      
      result_tibble <- tibble(
        SYMBOL = names(ensembl_ids),
        ENSEMBL = ensembl_ids
      )
      
      # Check for any gene symbols that didn't map and decide how to handle them
      not_missing_genes_indices <- which(!is.na(ensembl_ids))
      sce_hdf5 <- sce_hdf5[not_missing_genes_indices, ]
      # Update row names with Ensembl IDs
      rownames(sce_hdf5) <- ensembl_ids[not_missing_genes_indices]
      
      # For original counts, assay name should be "counts" if it is not 
      # For cpm, assay name should be "cpm" if it is not 
      # Rename assay names to X for consistency reason
      if (assays(sce_hdf5) |> names() == "X") {
        names(assays(sce_hdf5)) <- "counts"
      }
      
      if (assays(sce_hdf5) |> names() == "counts_per_million") {
        names(assays(sce_hdf5)) <- "cpm"
      }
      
      sce_hdf5
    },
    pattern = map(sce_hdf5),
    iteration = "list"
    ),
    tar_target(
      save_hdf5_to_anndata,
      convert_gene_to_ensembl_hdf5 |> zellkonverter::writeH5AD(conversion_tbl$anndata, compression = "gzip"),
      pattern = map(convert_gene_to_ensembl_hdf5, conversion_tbl),
      iteration = "list"
    )
  )


# tar_make(script = "~/git_control/CuratedAtlasQueryR/dev/convert_fibrosis_and_prostate_hdf5_to_anndata_targets.R",
#          store = "~/scratch/cellNexus/")
