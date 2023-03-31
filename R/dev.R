# Utility scripts for development purposes, that are not exported to users

#' Upload a file to the Nectar object store
#' @param source A character scalar indicating the local path to the file to
#'   upload
#' @param container A character scalar indicating the name of the container to
#'   upload to
#' @param name An optional character scalar indicating the name the file should
#'   have after being uploaded. Defaults to being the basename of the source
#'   file.
#' @param credential_id The OpenStack application credential ID as a character
#'   scalar. This is optional because you can alternatively source a
#'   `-openrc.sh` file instead of providing it here.
#' @param credential_id The OpenStack application credential secret as a
#'   character scalar
#' @return `NULL`, invisibly
#' @keywords internal
upload_swift <- function(
    source,
    container,
    name = basename(source),
    credential_id = NULL,
    credential_secret = NULL
) {
    # Create the basilisk environment
    swift_env <- basilisk::BasiliskEnvironment(
        envname="swift-nectar-upload",
        pkgname=packageName(),
        packages=c(
          "python-swiftclient==4.2.0",
          "python-keystoneclient==5.1.0",
          "python==3.10.9"
        )
    )
    proc <- basilisk::basiliskStart(swift_env)
    
    # Build the CLI args
    if (!is.null(credential_id) && !is.null(credential_secret)){
        auth <- c(
            "--os-auth-type",
            "v3applicationcredential",
            "--os-application-credential-id",
            credential_id,
            "--os-application-credential-secret",
            credential_secret
        )
    }
    else {
        auth <- character()
    }
    args <- c(
        "-m", 
        "swiftclient.shell",
        "--os-auth-url",
        "https://keystone.rc.nectar.org.au:5000/v3/",
        "--os-project-id",
        "06d6e008e3e642da99d806ba3ea629c5",
        auth,
        "upload",
        container,
        source,
        "--object-name",
        name
    )
    
    # Perform the upload
    system2(reticulate::py_exe(), args=args)
    basilisk::basiliskStop(proc)
    
    invisible(NULL)
}

#' Update the metadata database in nectar using a newly created data frame
#' @param metadata The data frame to upload
#' @param version The version for the new metadata as a character scalar, e.g.
#'   "0.2.3"
#' @inheritDotParams upload_swift
#' @examples
#' \dontrun{
#'  metadata = CuratedAtlasQueryR::get_metadata() |>
#'      head(10) |>
#'      dplyr::collect()
#'  update_database(
#'      metadata, 
#'      "0.2.3", 
#'      credential_id = "ABCDEFGHIJK", 
#'      credential_secret = "ABCD1234EFGH-5678IJK"
#'  )
#'  # Prints "metadata.0.2.3.parquet" if successful
#' }
#' @keywords internal
#' @inherit upload_swift return
update_database <- function(metadata, version, ...){
    # These are optional dev packages
    rlang::check_installed(c("arrow", "glue", "basilisk"))
    
    dir <- tempdir()
    parquet_name <- glue::glue("metadata.{version}.parquet")
    parquet_path <- file.path(dir, parquet_name)
    arrow::write_parquet(metadata, sink=parquet_path)
    
    upload_swift(parquet_path, container="metadata", name=parquet_name, ...)
}

#' Update the unharmonised parquet files
#' @param unharmonised_parquet_dir The path to a directory containing parquet
#'   files, one for each dataset, e.g.
#'   /vast/projects/cellxgene_curated/metadata_non_harmonised_parquet_0.2
#' @inheritDotParams upload_swift
#' @inherit upload_swift return
#' @keywords internal
#' @examples
#' \dontrun{
#' update_unharmonised(
#'     "/vast/projects/cellxgene_curated/metadata_non_harmonised_parquet_0.2", 
#'     credential_id = "ABCDEFGHIJK", 
#'     credential_secret = "ABCD1234EFGH-5678IJK"
#' )
#' }
update_unharmonised <- function(unharmonised_parquet_dir, ...){
    # name="/" forces it have no prefix, ie be at the top level in the bucket
    upload_swift(
        unharmonised_parquet_dir, 
        container="unharmonised_metadata",
        name="/", 
        ...
    )
}

#' Converts a series of HDF5Array-serialized SingleCellExperiments to AnnData
#' @param src A character scalar. The path to a directory containing one or more
#'  directories created by [HDF5Array::saveHDF5SummarizedExperiment()].
#' @param dest A character scalar. The path to a directory in which to save the
#'  created anndata files.
#' @keywords internal
#' @return A character vector of the newly-created anndata files
#' @examples
#' \donttest{
#' dir_to_anndata(
#'     "/vast/projects/cellxgene_curated/splitted_DB2_data_0.2.1",
#'     "/vast/projects/cellxgene_curated/splitted_DB2_anndata_0.2.1"
#' )
#' dir_to_anndata(
#'     "/vast/projects/cellxgene_curated/splitted_DB2_data_scaled_0.2.1",
#'     "/vast/projects/cellxgene_curated/splitted_DB2_anndata_scaled_0.2.1"
#' )
#' }
dir_to_anndata <- function(src, dest){
    dir.create(dest, showWarnings = FALSE)
    # This is a quick utility script to convert the SCE files into AnnData format for use in Pythonlist.files("/vast/projects/RCP/human_cell_atlas/splitted_DB2_data", full.names = FALSE) |>  purrr::walk(function(dir){
    basilisk::basiliskRun(fun = function(sce) {
        list.dirs(src)[-1] |>
            purrr::map_chr(function(sce_dir){
                cli::cli_alert_info("Processing {sce_dir}.")
                prefix <- basename(sce_dir)
                out_path <- glue::glue("{prefix}.h5ad") |>
                    file.path(dest, name=_)
                
                if (file.exists(out_path)) {
                    cli::cli_alert_info("{out_path} already exists. Skipping")
                }
                else {
                    sce <- HDF5Array::loadHDF5SummarizedExperiment(sce_dir)
                    single_column <- length(colnames(sce)) == 1
                    if (single_column){
                        # Hack, so that single-column SCEs will convert 
                        # correctly
                        cli::cli_alert_info(
                            "{sce_dir} has only 1 column. Duplicating column."
                        )
                        sce <- cbind(sce, sce)
                        single_column <- TRUE
                    }
                    ad <- zellkonverter::SCE2AnnData(sce)
                    if (single_column){
                        # Remove the duplicate column
                        sce$X <- sce$X[1]
                    }
                    # TODO: customize chunking here, when anndata supports it
                    # (see https://github.com/scverse/anndata/issues/961)
                    ad$write_h5ad(out_path)
                }
                out_path
            }, .progress = "Converting files")
    }, env = zellkonverter::zellkonverterAnnDataEnv())
}