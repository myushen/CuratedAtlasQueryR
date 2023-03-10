#' Environment that we use to cache the DuckDB connections
cache = rlang::env(
    metadata_table = rlang::env()
)


#' Gets the Curated Atlas metadata as a data frame.
#' 
#' Downloads a parquet database of the Human Cell Atlas metadata to a local 
#' cache, and then opens it as a data frame. It can then be filtered and 
#' passed into [get_SingleCellExperiment()] 
#' to obtain a [`SingleCellExperiment::SingleCellExperiment-class`]
#'
#' @param remote_url Optional character vector of length 1. An HTTP URL pointing
#'   to the location of the parquet database.
#' @param cache_directory Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   metadata.parquet
#' @param use_cache Optional logical scalar. If `TRUE` (the default), and this
#'  function has been called before with the same parameters, then a cached
#'  reference to the table will be returned. If `FALSE`, a new connection will
#'  be created no matter what.
#' @return A lazy data.frame subclass containing the metadata. You can interact
#'   with this object using most standard dplyr functions. For string matching,
#'   it is recommended that you use `stringr::str_like` to filter character
#'   columns, as `stringr::str_match` will not work.
#' @export
#' @examples
#' library(dplyr)
#' filtered_metadata <- get_metadata() |>
#'     filter(
#'         ethnicity == "African" &
#'             assay %LIKE% "%10x%" &
#'             tissue == "lung parenchyma" &
#'             cell_type %LIKE% "%CD4%"
#'     )
#'
#' @importFrom DBI dbConnect
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl
#' @importFrom httr progress
#' @importFrom cli cli_alert_info
#' 
#' @details 
#' 
#' The metadata was collected from the Bioconductor package `cellxgenedp`. it's vignette `using_cellxgenedp` provides an overview of the columns in the metadata.
#' The data for which the column `organism_name` included "Homo sapiens" was collected collected from `cellxgenedp`.
#' 
#' The columns `dataset_id` and `file_id` link the datasets explorable through `CuratedAtlasQueryR` and `cellxgenedp`to the CELLxGENE portal.
#' 
#'  Our representation, harmonises the metadata at dataset, sample and cell levels, in a unique coherent database table.
#' 
#' Dataset-specific columns (definitions available at cellxgene.cziscience.com)
#' `cell_count`, `collection_id`, `created_at.x`, `created_at.y`, `dataset_deployments`, `dataset_id`, `file_id`, `filename`, `filetype`, `is_primary_data.y`, `is_valid`, `linked_genesets`, `mean_genes_per_cell`, `name`, `published`, `published_at`, `revised_at`, `revision`, `s3_uri`, `schema_version`, `tombstone`, `updated_at.x`, `updated_at.y`, `user_submitted`, `x_normalization`
#' 
#' Sample-specific columns (definitions available at cellxgene.cziscience.com)
#' 
#' `sample_`, `.sample_name`, `age_days`, `assay`, `assay_ontology_term_id`, `development_stage`, `development_stage_ontology_term_id`, `ethnicity`, `ethnicity_ontology_term_id`, `experiment___`, `organism`, `organism_ontology_term_id`, `sample_placeholder`, `sex`, `sex_ontology_term_id`, `tissue`, `tissue_harmonised`, `tissue_ontology_term_id`, `disease`, `disease_ontology_term_id`, `is_primary_data.x`
#' 
#' Cell-specific columns (definitions available at cellxgene.cziscience.com)
#' 
#' `cell_`, `cell_type`, `cell_type_ontology_term_idm`, `cell_type_harmonised`, `confidence_class`, `cell_annotation_azimuth_l2`, `cell_annotation_blueprint_singler` 
#' 
#' Through harmonisation and curation we introduced custom column, not present in the original CELLxGENE metadata
#' 
#' - `tissue_harmonised`: a coarser tissue name for better filtering
#' - `age_days`: the number of days corresponding to the age
#' - `cell_type_harmonised`: the consensus call identity (for immune cells) using the original and three novel annotations using Seurat Azimuth and SingleR
#' - `confidence_class`: an ordinal class of how confident `cell_type_harmonised` is. 1 is complete consensus, 2 is 3 out of four and so on.             
#' - `cell_annotation_azimuth_l2`: Azimuth cell annotation
#' - `cell_annotation_blueprint_singler`: SingleR cell annotation using Blueprint reference
#' - `cell_annotation_blueprint_monaco`: SingleR cell annotation using Monaco reference
#' - `sample_id_db`: Sample subdivision for internal use
#' - `file_id_db`: File subdivision for internal use
#' - `sample_`: Sample ID
#' - `.sample_name`: How samples were defined
#' 
#' 
#' **Possible cache path issues**
#' 
#' If your default R cache path includes non-standard characters (e.g. dash because of your user or organisation name), the following error can manifest
#' 
#' Error in `db_query_fields.DBIConnection()`:
#' ! Can't query fields.
#' Caused by error:
#' ! Parser Error: syntax error at or near "/"
#' LINE 2: FROM /Users/bob/Library/Cach...
#' 
#' The solution is to choose a different cache, for example
#' 
#' get_metadata(cache_directory = path.expand('~'))
#' 
get_metadata <- function(
    remote_url = "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/metadata/metadata.0.2.3.parquet",
    cache_directory = get_default_cache_dir(),
    use_cache = TRUE
) {
    hash = c(remote_url, cache_directory) |> paste0(collapse="") |> cli::hash_sha256()
    cached_connection = cache$metadata_table[[hash]]
    if (!is.null(cached_connection) && isTRUE(use_cache)) {
        cached_connection
    }
    else {
        db_path <- file.path(cache_directory, "metadata.0.2.3.parquet")
        sync_remote_file(
            remote_url,
            db_path,
            progress(type = "down", con = stderr())
        )
        table = duckdb() |>
            dbConnect(drv = _, read_only = TRUE) |>
            tbl(db_path)
        cache$metadata_table[[hash]] = table
        table
    }
}