
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
