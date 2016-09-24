#' @name karyogram
#' @title karyogram
#' @description
#'
#' Returns gTrack displaying band pattern for hg18 or hg19 karyotype
#'
#' @param hg19 logical scalar flag, if T returns gTrack for hg19
#' @param bands logical scalar, if T returns gTrack with colored giemsa bands
#' @param arms logical scalar, if T and bands F, returns chromosome arms with different colors and centromeres and telomeres marked, if arms is F and bands F returns chromosomes, each with a different color
#' @param tel.width numeric scalar, specifies telomere width in bases (only relevant if arms = T, bands = F)
#' @param ... Additional arguments sent to the \code{gTrack} constructor
#' @export
#' @author Marcin Imielinski
karyogram = function(hg19 = TRUE, bands = TRUE, arms = TRUE, tel.width = 2e6, ... )
{
  #    ucsc.bands = get_ucsc_bands(hg19 = hg19)

  if (hg19)
    ucsc.bands = readRDS(system.file("extdata", "ucsc.bands.hg19.rds", package = 'gTrack'))
  else
    ucsc.bands = readRDS(system.file("extdata", "ucsc.bands.hg18.rds", package = 'gTrack'))

  if (bands)
  {
    si = si2gr(ucsc.bands);
    si = suppressWarnings(si[order(as.numeric(as.character(seqnames(si))))])
    si = Seqinfo(as.character(seqnames(si)), width(si)+1)
    ucsc.bands = sort(gr.fix(ucsc.bands, si, drop = T))

    ucsc.bands = split(ucsc.bands, seqnames(ucsc.bands))
    td = gTrack(list(ucsc.bands), colormaps = list(stain = c('gneg' = 'white', 'gpos25' = 'gray25', 'gpos50' = 'gray50', 'gpos75'= 'gray75', 'gpos100' = 'black', 'acen' = 'red', 'gvar' = 'pink', 'stalk' = 'blue')), border = 'black', ...)
    formatting(td)$stack.gap[1] = 0

  }
  else
  {
    tmp = c(RColorBrewer::brewer.pal(11, 'BrBG'), RColorBrewer::brewer.pal(11, 'PiYG'), RColorBrewer::brewer.pal(11, 'RdYlGn'))
    col.adjust = 10;

    col.karyo = data.frame(
      qtel = rep('black', length(tmp)),
      qarm = lighten(tmp, col.adjust),
      centro = rep('gray', length(tmp)),
      parm = lighten(tmp, -col.adjust),
      ptel = rep('black', length(tmp)), stringsAsFactors = F)

    if (arms) ## draw arms with slightly different hues of the same color and black ranges for telomere / centromere
    {
      tmp.tel = aggregate(formula = end ~ seqnames, data = GenomicRanges::as.data.frame(ucsc.bands), FUN = max)
      tmp.tel = structure(tmp.tel[,2], names = tmp.tel[,1])+1
      telomeres = c(GRanges(names(tmp.tel), IRanges::IRanges(start = rep(1, length(tmp.tel)), end = rep(tel.width, length(tmp.tel))),
                            seqlengths = GenomeInfoDb::seqlengths(seqinfo(ucsc.bands))),
                    GRanges(names(tmp.tel), IRanges::IRanges(tmp.tel-tel.width+1, tmp.tel), seqlengths = GenomeInfoDb::seqlengths(seqinfo(ucsc.bands))))
      values(telomeres)$lwd.border = 1;
      values(telomeres)$border = 'black';
      centromeres = reduce(ucsc.bands[values(ucsc.bands)$stain=='acen'])
      values(centromeres)$lwd.border = 1;
      values(centromeres)$border = 'black';
      #            arms = reduce(setdiff(ucsc.bands, centromeres))
      arms = IRanges::gaps(c(telomeres, centromeres))
      arms = arms[which(GenomicRanges::strand(arms)=='*')]
      values(arms)$lwd.border = 1;
      values(arms)$border = 'black';
      karyotype = sort(c(arms, centromeres, telomeres))
      values(karyotype)$region = paste(as.character(seqnames(karyotype)), c('ptel', 'p', ' cen', 'q', 'qtel'), sep = '')
      colmap = structure(as.vector(t(as.matrix(col.karyo)))[1:length(karyotype)], names = values(karyotype)$region)
      GenomicRanges::strand(karyotype) = '+'

      #            karyotype = split(karyotype, seqnames(karyotype))
      si = si2gr(karyotype);
      si = suppressWarnings(si[order(as.numeric(as.character(seqnames(si))))])
      si = Seqinfo(as.character(seqnames(si)), width(si)+1)
      karyotype = gr.fix(karyotype, si, drop = T)
      names(karyotype) = NULL;
      td = gTrack(list(karyotype), colormaps = list(region= colmap), border = 'black',
                  stack.gap = 0, xaxis.interval = 1e7, ...)
    }
    else
    {
      karyotype = sort(reduce(ucsc.bands))
      values(karyotype)$region = as.character(seqnames(karyotype))
      colmap = structure(col.karyo$qrm[1:length(karyotype)], names = as.character(seqnames(karyotype)))
      GenomicRanges::strand(karyotype) = '+'
      si = si2gr(karyotype);
      si = suppressWarnings(si[order(as.numeric(as.character(seqnames(si))))])
      si = Seqinfo(as.character(seqnames(si)), width(si)+1)
      karyotype = gr.fix(karyotype, si, drop = T)
      names(karyotype) = NULL;
      td = gTrack(list(karyotype), colormaps = list(region = colmap),
                  border = 'black', stack.gap = 0, xaxis.interval = 1e7, ...)
    }

  }

  formatting(td)$angle = 10
  formatting(td)$ywid = 2

  return(td)
}
