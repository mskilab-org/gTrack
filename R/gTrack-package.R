#' @import rtracklayer
#' @import Rsamtools
#' @importFrom bamUtils read.bam splice.cigar
#' @importFrom data.table data.table rbindlist := setkeyv
#' @importFrom GenomeInfoDb Seqinfo seqinfo keepSeqlevels seqlevels seqlengths seqlevels<- seqlengths<- genome<- seqnames
#' @importFrom GenomicRanges GRanges values ranges width strand values<- strand<- seqnames coverage ranges<- reduce seqinfo
#' @importFrom gUtils grl.unlist si2gr grbind gr.string gr.fix grl.pivot gr.findoverlaps gr.flatten gr.chr gr.match gr.sub
#' @importFrom methods setClass setGeneric setMethod setRefClass
#' @importFrom RCurl url.exists
NULL