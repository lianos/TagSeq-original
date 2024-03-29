## ============================================================================
## Use biosignals
## ============================================================================

##' Smoothes vector \code{x} using kernel density smoothing
##'
##' @author Anshul Kundaje (original MATLAB implementation)
##' @author Steve Lianoglou (ported to R)
##'
##' @param x input vector
##' @param kernel.type name of kernel
##' @param bandwidth kernel half width
##' @param normalize If set to true, kernel is equivalent to smooth averaging,
##' otherwise it is equivalent of a smooth sum
##' @param trim.length \code{logical(1)}. The vector returned by
##' \code{\link{convolve}} might be longer than the one you gave it. If
##' \code{TRUE}, the returned vector will be as long as your original vector.
##'
##' @return a vector of smoothed values for \code{x}
# kdeSmooth <- function(x, kernel, bandwidth=1L, normalize=FALSE,
#                       trim.length=TRUE, rescale=TRUE, upper.quantile=0.98,
#                       convolve.missing=FALSE, kernel.norm=NULL, ...) {
#   if (inherits(x, 'Rle')) {
#     return(kdeSmooth.Rle(x, kernel, bandwidth, normalize=normalize,
#                          trim.length=trim.length, ...))
#   }
#   stopifnot(is.numeric(x))
#   if (!is.null(dim(x))) {
#     ## No convultion of matrices yet
#     stop("x can only be 1D")
#   }
# 
#   if (is.character(kernel)) {
#     kernel <- match.arg(kernel, .convolutionKernels())
#     kernel.norm <- generateKernel(kernel, bandwidth, normalize=TRUE, ...)
#     kernel <- generateKernel(kernel, bandwidth, normalize=normalize, ...)
#   }
# 
#   ## Pad head/tail of vector to sidestep "edge effects"
#   x.padded <- c(rep(x[1], bandwidth), x, rep(tail(x, 1), bandwidth))
#   not.number <- !is.finite(x.padded)
#   has.nan <- any(not.number)
# 
#   convolve.missing <- convolve.missing && !is.null(kernel.norm)
#   if (convolve.missing) {
#     ## Guard against missing values
#     if (has.nan) {
#       missing.values <- x.padded[not.number]
#       x.padded[not.number] <- 0
#     }
#   }
# 
#   xs <- convolve(x.padded, rev(kernel), type='open')
# 
#   if (convolve.missing) {
#     ## Dividing by normalized result "deals with" missing values in x
#     xs <- xs / convolve(as.numeric(!not.number), rev(kernel.norm), type='open')
#   }
# 
#   ## Shift the signal "back" by to account for extra padding + bandwidth
#   xs <- xs[-seq(2*bandwidth)]
#   if (trim.length) {
#     ## The vector you get might be longer than the one you gave.
#     ## xs <- xs[1:min(length(x) + bandwidth, length(xs))]
#     xs <- xs[1:min(length(x), length(xs))]
#   }
# 
#   ## Reset original missing values
#   if (convolve.missing && has.nan) {
#     xs[not.number] <- missing.values
#   }
# 
#   ## TODO: Consider looking for + removing outliers at the tails of vector.
#   ##       Maybe look at the range of `diff`s you get at the edges and remove
#   ##       outliers. Maybe the lenght of tails to look for is a function of
#   ##       bandwidth -- explore further
#   if (rescale) {
#     xs <- (xs / max(xs, na.rm=TRUE)) * max(x, na.rm=TRUE)
#   }
# 
#   xs
# }

# kdeSmooth.Rle <- function(x, kernel, bandwidth, normalize=TRUE,
#                           trim.length=TRUE, rescale=TRUE, upper.quantile=0.98,
#                           ...) {
#   stopifnot(inherits(x, 'Rle'))
#   if (is.character(kernel)) {
#     kernel <- match.arg(kernel, .convolutionKernels())
#     kernel <- generateKernel(kernel, bandwidth, normalize=normalize, ...)
#   }
# 
#   regions <- slice(x, 1, rangesOnly=TRUE)
#   v <- Views(x, regions)
# 
#   ## Calculate new values in runs
#   vs <- viewApply(v, function(vv) {
#     vector <- as.numeric(vv)
#     ## conv <- convolve(padded, kval, type='open')
#     conv <- kdeSmooth(vector, kernel, bandwidth, rescale=rescale,
#                       trim.length=TRUE)
#     conv
#   }, simplify=FALSE)
# 
#   .as.Rle(vs, regions, length(x))
# }

##' Generates the kernel window function
##'
##' @author Anshul Kundaje (original MATLAB implementation)
##' @author Steve Lianoglou (ported to R)
##'
##' @param kernel.type Name of kernel to use
##' @param bandwidth kernel half width
##' @param normalize \code{logical(1)}: normalizes kernel if (\code{TRUE}).
##'
##' @return a vector with kernel values
# generateKernel <- function(kernel.type='normal', bandwidth=1L, normalize=TRUE, ...) {
#   kernel.type <- match.arg(kernel.type, .convolutionKernels())
#   if (bandwidth <= 0) {
#     stop("kernel bandwidth must be > 0")
#   }
#   win.len <- 2L * ceiling(bandwidth) + 1L
#   kval <- seq(-1, 1, length.out=win.len)
# 
#   args <- list(...)
#   mu <- if (!is.null(args$mu)) args$mu else 0
#   sd <- if (!is.null(args$sd)) args$sd else 1
# 
#   if (kernel.type %in% c('normal', 'gaussian')) {
#     kval <- dnorm(seq(-3, 3, length.out=win.len), mu, sd)
#   } else if (kernel.type %in% c('unform', 'rectangular')) {
#     kval <- rep(1, win.len)
#   } else if (kernel.type == 'triangular') {
#     kval <- 1 - abs(kval)
#   } else if (kernel.type == 'epanechnikov') {
#     kval <- (3/4) * (1 - kval^2)
#   } else if (kernel.type %in% c('quartic', 'biweight')) {
#     kval <- (15 / 16) * (1 - kval^2)^2
#   } else if (kernel.type %in% c('triweight', 'tricube')) {
#     kval <- (35/32) * (1 - kval^2)^3
#   } else if (kernel.type == 'cosine') {
#     kval <- (pi / 4) * cos(pi * kval / 2)
#   } else if (kernel.type == 'laplacian') {
#     kval <- c(1, -2, 1)
#     return(kval)
#   } else if (kernel.type == 'lofg') {
#     ## http://www.ce.rit.edu/~cockburn/courses/cvf05/Hw1_corrected.pdf
#     ## \frac{d}{dx}G = - \frac{2}{\sqrt{2\pi}\sigma^3} \times
#     ## \exp (- \frac{x^2}{\sigma^2} )
#     kval <- (2 / (sqrt(2*pi) * sd^3)) * kval * exp(-1 * kval^2 / sd^2)
#   } else if (kernel.type == 'tukey') {
#     stop("wut's a tukey?")
#   } else {
#     stop("wut kernel is this: ", kernel.type)
#   }
# 
#   ## if (normalize) {
#   ##   kval <- kval / sum(kval)
#   ## } else {
#   ##   kval <- kval / max(kval)
#   ## }
# 
#   kval
# }

##' @nord
##' The types of kernel we know how to use
# .convolutionKernels <- function() {
#   c('uniform', 'rectangular', 'gaussian', 'normal', 'triangular',
#     'epanechnikov', 'quartic', 'biweight', 'tricube', 'triweight',
#     'cosine', 'tukey', 'laplacian', 'lofg')
# }

# .fixEdgeEffects <- function(x, bandwidth, is.normed=FALSE) {
#   if (is.normed) {
#     xn <- x
#   } else {
#     xn <- xn / max(xn)
#   }
# 
#   edge.size <- ceiling(length(x) * 0.1)
#   ## head
#   ## idx.head <- 1:edge.size
#   ## diff.head <- diff(x[idx.head])
#   ## fix.head <- max(diff.head) > 1.5 * sd(diff.head)
#   ## if (any(fix.head)) {
#   ##   x[fix.head] <-
#   ## }
# 
#   ## idx.tail <- (length(x) - edge.size + 1L):length(x)
#   ## diff.tail <- diff(x[idx.tail])
#   ## fix.tail <- max(diff.tail) > 1.5 * sd(diff.tail)
#   ## if (any(fix.tail)) {
# 
#   ## }
# }

# .as.Rle <- function(new.vals, val.ranges, len) {
#   if (is.unsorted(start(val.ranges)) || is.unsorted(end(val.ranges))) {
#     stop("Illegal ranges")
#   }
#   o <- findOverlaps(val.ranges, ignoreSelf=TRUE, ignoreRedundant=TRUE,
#                     type="any")
#   if (length(o) > 0L) {
#     stop("Overlapping ranges not handled yet")
#   }
# 
#   prev.end <- 0L
#   runs <- lapply(1:length(new.vals), function(idx) {
#     ## Builds 2 column matrix. Col 1 is runValue, col2 is runLength
#     val <- new.vals[[idx]]
#     .range <- val.ranges[idx]
# 
#     ans <- cbind(val, rep(1, length(val)))
#     dist <- start(.range) - prev.end -1L
#     if (dist > 0) {
#       ans <- rbind(c(0, dist), ans)
#     }
#     prev.end <<- end(.range)
#     ans
#   })
# 
#   end.point <- end(val.ranges[length(val.ranges)])
#   if (end.point < len) {
#     runs[[length(runs) + 1L]] <- c(0, len - end.point + 1)
#   }
#   runs <- do.call(rbind, runs)
# 
#   smoothed <- Rle(runs[,1], runs[,2])
#   smoothed
# }

################################################################################
## At one point in time, I was thinking about laplacian of gaussians
## Using Laplacian of Gaussians to call peaks:
## http://homepages.inf.ed.ac.uk/rbf/HIPR2/log.htm
## http://www.cse.shirazu.ac.ir/~eghbali/cee1/log.html
## For derviation(?):
## http://fourier.eng.hmc.edu/e161/lectures/gradient/node8.html
##
##
##   The LoG operator calculates the second spatial derivative of an image.
##   This means that in areas where the image has a constant intensity
##   (i.e. where the intensity gradient is zero), the LoG response will be zero.
##   In the vicinity of a change in intensity, however, the LoG response will be
##   positive on the darker side, and negative on the lighter side. This means
##   that at a reasonably sharp edge between two regions of uniform but
##   different intensities, the LoG response will be:
##     * zero at a long distance from the edge,
##     * positive just to one side of the edge,
##     * negative just to the other side of the edge,
##     * zero at some point in between, on the edge itself.

##' assumes smooth/continuous function that doesn't "bounce off" 0
# .cross0 <- function(x, as.fence.post=FALSE) {
#   pos <- x > 0
#   neg <- x < 0
#   pos2neg <- head(pos, -1) & tail(neg, -1)
#   neg2pos <- head(neg, -1) & tail(pos, -1)
# 
#   pos2neg.idx <- which(pos2neg)
#   neg2pos.idx <- which(neg2pos)
# 
#   if (as.fence.post) {
# 
#   } else{
# 
#   }
# 
#   ret
# }

# findLoGPeaks <- function(LoG.values) {
#   pos <- LoG.values > 0
#   neg <- LoG.values < 0
# 
#   pos2neg <- head(pos, -1) & tail(neg, -1)
#   neg2pos <- head(neg, -1) & tail(pos, -1)
# 
#   pos2neg.idx <- which(pos2neg)
#   neg2pos.idx <- which(neg2pos)
# 
#   if (neg2pos.idx[1] < pos2neg.idx[1]) {
#     ## first "peak" is a trough(?)
#     neg2pos.idx <- neg2pos.idx[-1]
#   }
# 
#   if (length(neg2pos.idx) < length(pos2neg.idx)) {
#     ## This shouldn't happen
#     warning("Trimming `starts` in LoG", immediate.=TRUE)
#     pos2neg.idx <- pos2neg.idx[1:neg2pos.idx]
#   }
# 
#   IRanges(pos2neg.idx, neg2pos.idx[1:length(pos2neg.idx)])
# }
