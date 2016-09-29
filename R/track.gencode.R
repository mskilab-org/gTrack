#' @name track.gencode
#' @title Constructor a \code{gTrack} from GENCODE transcripts
#'
#' @description
#' Returns gTrack object representing GENCODE transcripts and their components (utr, cds etc) with assigned colors.
#' Usually built from cached data objects but can also be built from provided GRangesList
#'
#' @param rg (optional) GRangesList representing transcript models imported from GENCODE gff3 file using rtracklayer import
#' @param genes (optional) character vector specifying genes to limit gTrack object to
#' @param gene.collapse scalar logical specifying whether to collapse genes by transcript (or used stored version of transcripts)
#' @param grep character vector for which to grep genes to positively select
#' @param grepe character vector for which to grep genes which to exclude
#' @param bg.col scalar character representing background color for genespan
#' @param cds.col scalar character representing background color for CDS
#' @param cds.utr scalar character representing background color for UTR
#' @param st.col scalar character representing color of CDS start
#' @param en.col scalar character representing color of CDS end
#' @param cached logical scalar whether to use "cached" version provided with package
#' @param gr.srt.label scalar numeric specifying angle on exon label
#' @param gr.cex.label scalar numeric > 0 specifying character expansion on exon label
#' @param labels.suppress.gr scalar logical specifying whether to suppress exon label plotting
#' @param stack.gap stack.gap argument to gTrack
#' @param ... additional arguments passed down to gTrack
#'
#' @export
#' @author Marcin Imielinski
track.gencode = function(gencode = NULL,
                         gene.collapse = TRUE,
                         genes = NULL,
                         grep = NULL,
                         grepe = NULL, ## which to exclude
                         bg.col = alpha('blue', 0.1), cds.col = alpha('blue', 0.6), utr.col = alpha('purple', 0.4),
                         st.col = 'green',
                         en.col = 'red',
                         grl.labelfield, ## Don't touch these
                         gr.labelfield,
                         col,
                         cached = T, ## if true will use cached version
                         cached.dir = Sys.getenv('GENCODE_DIR'),
                         cached.path = paste(cached.dir, "gencode.composite.rds", sep = '/'),  ## location of cached copy
                         cached.path.collapsed = paste(cached.dir, "gencode.composite.collapsed.rds", sep = '/'),  ## location of cached copy
                         gr.srt.label = 0,
                         gr.cex.label = 0.3,
                         cex.label = 0.5,
                         #labels.suppress.gr = T,
                         drop.rp11 = TRUE,
                         stack.gap = 1e6,
                         ...)
{

  if (nchar(cached.dir)==0)
  {
    cached.path = system.file("extdata", "gencode.composite.rds", package = 'gTrack')  ## location of cached copy
    cached.path.collapsed = system.file("extdata", "gencode.composite.collapsed.rds", package = 'gTrack')
  }


  if (!gene.collapse)
  {
    if (file.exists(cached.path))
      cat(sprintf('Pulling gencode annotations from %s\n', cached.path))
    else
      cat(sprintf('Caching gencode annotations to %s\n', cached.path))
  }
  else
  {
    if (file.exists(cached.path.collapsed))
      cat(sprintf('Pulling gencode annotations from %s\n', cached.path.collapsed))
    else
      cat(sprintf('Caching gencode annotations to %s\n', cached.path.collapsed))
  }


  if (!cached | (!gene.collapse  & !file.exists(cached.path)) | (gene.collapse  & !file.exists(cached.path.collapsed)))  ## if no composite refgene copy, then make from scratch
  {
    cat('recreating composite gencode object\n')
    if (is.null(gencode))
      gencode = read_gencode()

    cat('loaded gencode rds\n')

    tx = gencode[gencode$type =='transcript']
    genes = gencode[gencode$type =='gene']
    exons = gencode[gencode$type == 'exon']
    utr = gencode[gencode$type == 'UTR']
    ## ut = unlist(utr$tag)
    ## utix = rep(1:length(utr), sapply(utr$tag, length))
    ## utr5 = utr[unique(utix[grep('5_UTR',ut)])]
    ## utr3 = utr[unique(utix[grep('3_UTR',ut)])]
    ## utr5$type = 'UTR5'
    ## utr3$type = 'UTR3'
    startcodon = gencode[gencode$type == 'start_codon']
    stopcodon = gencode[gencode$type == 'stop_codon']
    OUT.COLS = c('gene_name', 'transcript_name', 'transcript_id', 'type', 'exon_number', 'type')
    tmp = c(genes, tx, exons, utr, startcodon, stopcodon)[, OUT.COLS]

    cat('extracted intervals\n')

    ## compute tx ord of intervals
    ord.ix = order(tmp$transcript_id, match(tmp$type, c('gene', 'transcript', 'exon', 'UTR', 'start_codon','stop_codon')))
    tmp.rle = rle(tmp$transcript_id[ord.ix])
    tmp$tx.ord[ord.ix] = unlist(lapply(tmp.rle$lengths, function(x) 1:x))
    tmp = tmp[order(match(tmp$type, c('gene', 'transcript', 'exon', 'UTR', 'start_codon','stop_codon')))]

    cat('reordered intervals\n')

    tmp.g = tmp[tmp$type != 'transcript']
    gencode.composite = split(tmp.g, tmp.g$gene_name)
    ix = sapply(split(1:length(tmp.g), tmp.g$gene_name), function(x) x[1])
    values(gencode.composite)$id = tmp.g$gene_name[ix]
    values(gencode.composite)$gene_sym = tmp.g$gene_name[ix]
    values(gencode.composite)$type = tmp.g$gene_type[ix]
    values(gencode.composite)$status = tmp.g$gene_status[ix]
    cat('saving gene collapsed track\n')
    saveRDS(gencode.composite, cached.path.collapsed)

    tmp.t = tmp[tmp$type != 'gene']
    gencode.composite = split(tmp.t, tmp.t$transcript_id)
    ix = sapply(split(1:length(tmp.t), tmp.t$transcript_id), function(x) x[1])
    values(gencode.composite)$gene_sym = tmp.t$gene_name[ix]
    values(gencode.composite)$id = paste(tmp.t$gene_name[ix], tmp.t$transcript_name[ix], sep = '-')
    values(gencode.composite)$type = tmp.t$transcript_type[ix]
    values(gencode.composite)$status = tmp.t$transcript_status[ix]
    cat('saving transcript track\n')
    saveRDS(gencode.composite, cached.path)

    genes = NULL

    cat(sprintf('cached composite tracks at %s and %s\n', cached.path, cached.path.collapsed))
  }

  if (gene.collapse)
    gencode.composite = readRDS(cached.path.collapsed)
  else
    gencode.composite = readRDS(cached.path)


  if (drop.rp11)
    gencode.composite = gencode.composite[!grepl('^RP11', values(gencode.composite)$gene_sym)]

  if (!is.null(genes))
    gencode.composite = gencode.composite[values(gencode.composite)$gene_sym %in% genes]

  if (!is.null(grep))
  {
    ix = rep(FALSE, length(gencode.composite))
    for (g in grep)
      ix = ix | grepl(g, values(gencode.composite)$gene_sym)

    gencode.composite = gencode.composite[ix]
  }

  if (!is.null(grepe))
  {
    ix = rep(TRUE, length(gencode.composite))
    for (g in grepe)
      ix = ix & !grepl(g, values(gencode.composite)$gene_sym)

    gencode.composite = gencode.composite[ix]
  }


  cmap = list(type = c(gene = bg.col, transcript = bg.col, exon = cds.col, start_codon = st.col, stop_codon = en.col, UTR = utr.col))

  ## set the label params
  lab <- glabel(grl.labelfield = 'id', gr.labelfield = 'exon_number', gr.angle = gr.srt.label,
                gr.cex = cex.label, grl.cex = cex.label)

  return(suppressWarnings(gTrack(gencode.composite, col = NA, label = lab,
                                 stack.gap = stack.gap, colormaps = cmap, ...)))
}
