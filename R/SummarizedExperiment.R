annotateSummarizedExperiment <- function(x, annotations, nuclear=TRUE) {
  stopifnot(inherits(x, "SummarizedExperiment"))
  stopifnot(inherits(annotations, "GenomicRanges"))

  regions <- rowData(x)
  xmeta <- mcols(regions)
  mcols(regions) <- NULL

  anno.cols <- names(mcols(annotations))

  if (nuclear) {
    for (wut in anno.cols) {
      if (wut %in% names(xmeta)) {
        xmeta[[wut]] <- NULL
      }
    }
  }

  values(regions) <- xmeta
  suppressWarnings(regions <- annotateReads(regions, annotations))
  rowData(x) <- regions
  x
}

setMethod("TPM", c("SummarizedExperiment", "numeric"),
function(x, normer, do.group=!missing(group.cols),
         group.cols=c('seqnames', 'strand', 'entrez.id'),
         meta.cols=c(group.cols, "symbol"),
         assay.idx=1, na.rm.keys=TRUE, ...) {
  expt.names <- colnames(x)
  if (is.null(expt.names) || length(expt.names) != ncol(x)) {
    stop("Experiment names are required for the SummarizedExperiment")
  }
  if (!all(expt.names %in% names(normer))) {
    stop("Library sizes are not present for all expt.names")
  }

  meta <- as(rowData(x), 'data.table')
  if (!all(meta.cols %in% names(meta))) {
    stop("Requested meta.cols not in metadata")
  }

  meta <- meta[, meta.cols, with=FALSE]
  dt <- cbind(meta, data.table(assay(x, assay.idx)))

  if (isTRUE(do.group)) {
    if (!all(group.cols %in% names(meta))) {
      stop("group.cols not found in metadata")
    }
    if (na.rm.keys) {
      ## Remove any rows that have NA's in the keys
      keep <- Reduce("&", lapply(group.cols, function(kcol) !is.na(dt[[kcol]])))
      n.axe <- sum(!keep)
      if (n.axe) {
        warning("Removing ", n.axe, " rows due to NA in key column", .immediate=TRUE)
        dt <- dt[keep,]
      }
    }
    setkeyv(dt, group.cols)
    more.cols <- setdiff(meta.cols, group.cols)
    dt <- dt[, {
      cnts <- sapply(expt.names, function(wut) sum(.SD[[wut]]), simplify=FALSE)
      more <- sapply(more.cols, function(wut) .SD[[wut]][1], simplify=FALSE)
      c(more, cnts)
    }, by=group.cols]
  }

  tpmf <- getMethod("TPM", c("numeric", "numeric"))

  for (wut in expt.names) {
    dt[[wut]] <- tpmf(dt[[wut]], normer[[wut]])
  }

  as.data.frame(dt)
})

setMethod("TPM", c("SummarizedExperiment", "missing"),
function(x, normer, assay.idx=1, ...) {
  TPM(x, colSums(assay(x, assay.idx)), assay.idx=assay.idx, ...)
})

