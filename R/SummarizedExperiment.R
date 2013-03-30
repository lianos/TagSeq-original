##' Cluster peaks that are close together, optinaly add reads together.
##'
##' @param x The SummarizedExperiment to do this to
##' @param max.distance The max distance to consider combining peaks together
##' @param peak.pos Where to measure the distance between peaks from? Currently
##' only support end,start,center. edges is to use 'end' to measure dist
##' downstream, and start to measure dist upstream
##' @param combine.counts add the read count from each merged peak togther?
clusterPeaks <- function(x, max.distance=50L,
                         peak.pos=c('end', 'start', 'center', 'edges'),
                         requantify=TRUE) {
  peak.pos <- match.arg(peak.pos)
  if (peak.pos == 'edges') {
    stop("`edges` not supported yet")
  }
  stopifnot(inherits(x, 'SummarizedExperiment'))
  stopifnot(is.numeric(max.distance) && max.distance > 0L)

  xevts <- resize(rowData(x), width=max.distance, fix=peak.pos)
  mcols(xevts) <- NULL
  mask <- reduce(xevts[order(xevts)])

  mm <- as.matrix(findOverlaps(xevts, mask))
  if (any(duplicated(mm[,1]))) {
    stop("There shouldn't be any duplicate queryHits")
  }

  ## Reasign fenceposts from original peak definitions
  xx <- rowData(x)
  mcols(xx) <- NULL
  xx <- as(xx, 'data.table')

  mask.idx <- integer(nrow(mm))
  mask.idx[mm[,1]] <- mm[,2]
  xx[, mask.idx := mask.idx]
  setkeyv(xx, 'mask.idx')

  bounds <- xx[, {
    list(seqnames=seqnames[1L], start=min(start), end=max(end),
         strand=strand[1L], npeaks=.N)
  }, by="mask.idx"]

  out <- as(bounds[, mask.idx := NULL], 'GRanges')
  if (requantify) {
    out <- newSEfromRegions(out, x)
  }
  out
}

newSEfromRegions <- function(regions, original, original.assay.idx=1L) {
  stopifnot(inherits(regions, 'GenomicRanges'))
  stopifnot(inherits(original, 'SummarizedExperiment'))

  ocnts <- assay(original, original.assay.idx)
  mm <- as.matrix(findOverlaps(regions, rowData(original)))

  dt <- data.table(region.id=mm[,1], as.matrix(ocnts), key='region.id')
  ncnts <- dt[, lapply(.SD, sum), by='region.id']

  se <- SummarizedExperiment(as.matrix(ncnts[, region.id := NULL]), rowData=regions,
                             colData=colData(original), exptData=exptData(original))
  names(assays(se)) <- 'counts'
  se
}


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

