#' Run microbiome personalized hyper-graph ranking (mPHGR)
#'
#' Call the bundled, unmodified computational logic in Python through reticulate.
#' Input matrices must already be prepared. No preprocessing is performed.
#'
#' @param data_dir Folder with control.txt and numbered test files.
#' @param data_type "1" for longitudinal (PTB-style) normalization, "2" for
#'   stage-aggregated (GC-style) normalization and plotting.
#' @param d Minimum association count at the fixed 0.3 cutoff.
#' @param Q Differential correlation threshold for hyper-graph construction.
#' @param prefix_data File prefix before each numeric ID, e.g. "test2_data".
#' @param stages Named list of stage names and file IDs, for data_type="2".
#'   First stage supplies the normalization baseline.
#' @return Invisibly, a list matching the original Python return value.
#'   Results are written to Score.csv and d.csv in data_dir.
#' @export
mphgr_run <- function(data_dir, data_type, d, Q, prefix_data, stages = NULL) {
  python_folder <- system.file("python", package = "mPHGR")
  core <- reticulate::import_from_path(
    "mPHGR", path = python_folder, convert = TRUE
  )
  stages_py <- if (is.null(stages)) NULL else reticulate::r_to_py(stages)
  previous_dir <- getwd()
  on.exit(setwd(previous_dir), add = TRUE)
  result <- core$run_mphgr(
    data_dir = normalizePath(data_dir, winslash = "/", mustWork = TRUE),
    data_type = as.character(data_type),
    d = as.integer(d),
    Q = as.numeric(Q),
    prefix_data = as.character(prefix_data),
    stages = stages_py
  )
  invisible(result)
}
