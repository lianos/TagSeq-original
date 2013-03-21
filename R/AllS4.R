##' @exportClass TagSeqExperiment
setClass("TagSeqExperiment", contains="SummarizedExperiment")

##' @exportClass SAGEseqExperiment
setClass("SAGEseqExperiment", contains="TagSeqExperiment")

##' @exportClass CompressedReads
setClass("CompressedReads", contains="GRanges")

## An (annotated) compressed container for the tags (per chromosome).
## The values() DataFrame will store:
##  * strand (factor)
##  * chr (Rle<character>)
##  * is.unique (Rle<logical>)
##  * tag.type (tagType factor)
##  * symbol (character)
##  * tag.type.anti (character) -- the annotation on the opposite strand
##  * symbol.anti (character) -- the annotated gene on the opposite strand
##  * tag.type.extra -- add 1,2,3 for splitting of utrs?
## setClass("AnnotatedTags", contains="TagSequences",
##          representation(.tags='GRanges'),
##          prototype(.tags=GRanges()))

## -----------------------------------------------------------------------------
## Methods


