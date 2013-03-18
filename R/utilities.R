checkVerbose <- function(...) {
  verbose <- list(...)$verbose
  if (is.null(verbose)) verbose <- options()$verbose
  verbose
}

dir.exists <- function(path) {
  path <- as.character(path)
  sapply(file.info(path)$isdir, isTRUE)
}

