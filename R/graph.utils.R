#' @name barbs
#' @title barbs
#' @description
#'
#' (improvement over pencils)
#' used to represent directional ranges.  Shape consists of two horizontally reflected parallelograms
#' with a specified angle.  The middle segment of the parallelogram represents the provided range.
#'
#' Parameterized by x0, x1 (range values), y0, y1 (thickness), angle, where angle is between -90 and 90
#'
#' eg angle of 0 = rectangle
#'    angle of 45 = rightward pointing barb
#'    angle of -45 = leftward pointing barb
#'
#' Note: the shape will NOT be necessarily contained inside the box outlined by x0, x1, y0, y1
#' @author Marcin Imielinski
#' @keywords internal
barbs = function(x0, y0, x1, y1, angle = 0, ...)
{

  y.thick = (y1-y0)/2;
  y.midpoint = rowMeans(cbind(y0, y1))
  angle = pmax(-80, pmin(80, angle))/180*pi;

  # to make the angle visible irrespective of plot coordinates
  # we need to compute the x.offset in viewer coordinates (ie inches)
  x.wid = diff(par('usr')[1:2])
  y.wid = diff(par('usr')[3:4])
  x.offset = tan(angle)*y.thick/y.wid*par('pin')[2]*x.wid/par('pin')[1]

  x.vertices = rbind(x0, x0-x.offset, x1-x.offset, x1, x1-x.offset, x0-x.offset, x0, NA)
  y.vertices = rbind(y.midpoint, y1, y1, y.midpoint, y0, y0, y.midpoint, NA);

  if (ncol(x.vertices)==1)
    x.vertices = rep(x.vertices, ncol(y.vertices))
  else if (ncol(y.vertices)==1)
    y.vertices = rep(y.vertices, ncol(x.vertices))

  ## convert lwd.border 0 to NA. JEREMIAH
  other.args <- list(...)
  if (!any(other.args$lwd > 0)) {
    other.args$border = NA
    other.args[['lwd']] <- NULL
  }

  do.call('polygon', c(list(x=as.numeric(x.vertices), y=as.numeric(y.vertices)), other.args[names(other.args) %in% c("lwd", "border", "col")]))
}

#' @name bernsteinp
#' @title bernsteinp
#' @description
#' bernsteinp
#'
#' computes matrix of bernstein polynomials of degree n at m points t in the
#' interval 0, 1
#'
#' out[i,j] = choose(n, j)*t[i]^j*(1-t[i])^(n-j)
#' where t = seq(0, 1, length.out = m) and j in 0 .. n
#'
#' this matrix can be multiplied by a n x 2 matrix of control points to yield
#' an m x 2 matrix of m equally spaced points along the curve
#' @author Marcin Imielinski
#' @keywords internal
bernsteinp = function(n, m)
{
  m = round(m);
  if (m<0)
    return(NULL)
  t = seq(0, 1, length.out = m)
  out = matrix(sapply(0:(n-1), function(i) choose(n-1,i)*t^i*(1-t)^(n-1-i)), nrow = m);
  return(out)
}

#' @name affine.map
#' @title affine.map
#' @description
#'
#'
#' affinely maps 1D points in vector x from interval xlim to interval ylim,
#' ie takes points that lie in
#' interval xlim and mapping onto interval ylim using linear / affine map defined by:
#' (x0,y0) = c(xlim(1), ylim(1)),
#' (x1,y1) = c(xlim(2), ylim(2))
#' (using two point formula for line)
#' useful for plotting.
#'
#' if cap.max or cap.min == T then values outside of the range will be capped at min or max
#' @rdname affine-map-methods
#' @author Marcin Imielinski
#' @keywords internal
affine.map = function(x, ylim = c(0,1), xlim = c(min(x), max(x)), cap = F, cap.min = cap, cap.max = cap, clip = T, clip.min = clip, clip.max = clip)
{
  #  xlim[2] = max(xlim);
  #  ylim[2] = max(ylim);

  if (xlim[2]==xlim[1])
    y = rep(mean(ylim), length(x))
  else
    y = (ylim[2]-ylim[1]) / (xlim[2]-xlim[1])*(x-xlim[1]) + ylim[1]

  if (cap.min)
    y[x<min(xlim)] = ylim[which.min(xlim)]
  else if (clip.min)
    y[x<min(xlim)] = NA;

  if (cap.max)
    y[x>max(xlim)] = ylim[which.max(xlim)]
  else if (clip.max)
    y[x>max(xlim)] = NA;

  return(y)
}


#' @name alpha
#' @title alpha
#' @description
#' Give transparency value to colors
#'
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @author Marcin Imielinski
#' @param col RGB color
#' @keywords internal
alpha = function(col, alpha)
{
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}

#' Blends colors
#'
#' @param cols colors to blend
#' @keywords internal
blend = function(cols)
{
  col.rgb = rowMeans(col2rgb(cols))
  out = rgb(red = col.rgb['red']/255, green = col.rgb['green']/255, blue = col.rgb['blue']/255)
  return(out)
}


#' @name col.scale
#' @title col.scale
#' @description
#'
#' Assign rgb colors to numeric data
#'
#' Assigns rgb colors to numeric data values in vector "x".. maps scalar values
#' in val.range (default c(0,1)) to a linear color scale of between col.min (default white)
#' and col.max (default black), each which are length 3 vectors or characters.  RGB values are scaled between 0 and 1.
#' Values below and above val.min and val.max are mapped to col.max and col.max respectively
#'
#' @author Marcin Imielinski
#' @keywords internal
col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white',
                     invert = F # if T flips rgb.min and rgb.max
)
{

  ## NOTE fix
  error = NULL

  if (!is.numeric(col.min))
    if (is.character(col.min))
      col.min = col2rgb(col.min)/255
    else
      error('Color should be either length 3 vector or character')

    if (!is.numeric(col.max))
      if (is.character(col.max))
        col.max = col2rgb(col.max)/255
      else
        error('Color should be either length 3 vector or character')

      col.min = as.numeric(col.min);
      col.max = as.numeric(col.max);

      x = (pmax(val.range[1], pmin(val.range[2], x))-val.range[1])/diff(val.range);
      col.min = pmax(0, pmin(1, col.min))
      col.max = pmax(0, pmin(1, col.max))

      if (invert)
      {
        tmp = col.max
        col.max = col.min
        col.min = tmp
      }

      nna = !is.na(x);

      out = rep(na.col, length(x))
      out[nna] = rgb((col.max[1]-col.min[1])*x[nna] + col.min[1],
                     (col.max[2]-col.min[2])*x[nna] + col.min[2],
                     (col.max[3]-col.min[3])*x[nna] + col.min[3])

      return(out)
}

#' @name brewer.master
#' @title brewer.master
#' @description
#' Make brewer colors using entire palette
#'
#' Makes a lot of brewer colors using entire brewer palette
#' @param scalar positive integer of number of colors to return
#' @param palette any brewer.pal palette to begin with (default 'Accent')
#'
#' @keywords internal
#' @author Marcin Imielinski
brewer.master = function(n, palette = 'Accent')
{
  palettes = list(
    sequential = c('Blues'=9,'BuGn'=9, 'BuPu'=9, 'GnBu'=9, 'Greens'=9, 'Greys'=9, 'Oranges'=9, 'OrRd'=9, 'PuBu'=9, 'PuBuGn'=9, 'PuRd'=9, 'Purples'=9, 'RdPu'=9, 'Reds'=9, 'YlGn'=9, 'YlGnBu'=9, 'YlOrBr'=9, 'YlOrRd'=9),
    diverging = c('BrBG'=11, 'PiYG'=11, 'PRGn'=11, 'PuOr'=11, 'RdBu'=11, 'RdGy'=11, 'RdYlBu'=11, 'RdYlGn'=11, 'Spectral'=11),
    qualitative = c('Accent'=8, 'Dark2'=8, 'Paired'=12, 'Pastel1'=8, 'Pastel2'=8, 'Set1'=9, 'Set2'=8, 'Set3'=12)
  );

  palettes = unlist(palettes);
  names(palettes) = gsub('\\w+\\.', '', names(palettes))

  if (palette %in% names(palettes))
    i = match(palette, names(palettes))
  else
    i = ((max(c(1, suppressWarnings(as.integer(palette))), na.rm = T)-1) %% length(palettes))+1

  col = c();
  col.remain = n;

  while (col.remain > 0)
  {
    if (col.remain > palettes[i])
    {
      next.n = palettes[i]
      col.remain = col.remain-next.n;
    }
    else
    {
      next.n = col.remain
      col.remain = 0;
    }

    col = c(col, RColorBrewer::brewer.pal(max(next.n, 3), names(palettes[i])))
    i = ((i) %% length(palettes))+1
  }

  col = col[1:n]
  return(col)
}


#' @name lighten
#' @title lighten
#' @description
#' lighten
#'
#' lightens / darkens colors by brighness factor f (in -255 .. 255) that will make lighter if > 0 and darker < 0
#' @author Marcin Imielinski
#' @keywords internal
lighten = function(col, f)
{
  M = col2rgb(col)
  return(apply(matrix(pmax(0, pmin(255, M + f*matrix(rep(1, length(M)), nrow = nrow(M)))), ncol = length(col))/255, 2, function(x) rgb(x[1], x[2], x[3])))
}


#' @name plot.blank
#' @title plot.blank
#' @description
#' Make a blank plot
#'
#' Shortcut for making blank plot with no axes
#' @author Marcin Imielinski
#' @keywords internal
plot.blank = function(xlim = c(0, 1), ylim = c(0,1), xlab = "", ylab = "", axes = F, bg.col = "white", ...)
{
  par(bg = bg.col)
  plot(0, type = "n", axes = axes, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  #    par(usr = c(xlim, ylim))
}

#' Function to plot color bar
#' http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
#' @keywords internal
color.bar <- function(lut, nticks=11, ticks=seq(min, max, len=nticks), title='',
                      xpos, ypos, width, height, text=c("min", 'max'), cex=1) {
  scale = (length(lut)-1)/(height)
  #scale = 1
  #dev.new(width=1.75, height=5)
  #plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  #axis(2, ticks, las=1)
  rect(xpos, ypos, xpos+width, ypos + height, lwd=2)
  #text(xpos+width, ypos, text[1])
  #text(xpos+width, ypos+height, text[2])
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + ypos
    rect(xpos, y, width, y+1/scale, col=lut[i], border=NA)
  }

  ## add the text
  for (i in 1:length(text)) {
    y = ypos + (i-1)/(length(text)-1) * height
    segments(xpos+width, y, xpos+width*1.3, y)
    text(xpos+width*2.5, y, text[i], cex=cex)
  }
}


#' utility function to load data files into variables
#' @keywords internal
readData = function(..., environment = NULL) {
  my.env  = new.env()
  data(..., envir = my.env);
  return(as.list(my.env))
}

#' @keywords internal
gr.flatten = function(gr, gap = 0)
{
  if (length(gr) == 0)
    return(data.frame())
  else if (length(gr) == 1)
    return(data.frame(start = 1, end = width(gr) + 1)) ## added +1 for [1,2] is width 2
  else
  {
    starts = as.numeric(cumsum(c(1, width(gr[1:(length(gr)-1)])+gap)))
    #ends = as.numeric(starts+width(gr)-1) ## [1,2] is width 1
    ends = as.numeric(starts+width(gr)) ## [1,2] is width 2
    return(cbind(data.frame(start=starts, end=ends), as.data.frame(mcols(gr))))
  }
}


#' gr.flatmap
#'
#' Takes \code{GRanges} pile and maps onto a flattened coordinate system defined by windows (GRanges object)
#' a provided "gap" (in sequence units).  If squeeze == TRUE then will additionally squeeze ranges into xlim.
#'
#' output is list with two fields corresponding to data frames:
#' $grl.segs = data frame of input gr's "lifted" onto new flattened coordinate space (NOTE: nrow of this not necessarily equal to length(gr))
#' $window.segs = the coordinates of input windows in the new flattened (and squeezed) space
#' @param gr \code{GRanges} pile to flatten
#' @param windows \code{GRanges} pile of windows defining the coordinate system
#' @param gap TODO
#' @param strand.agnostic TODO
#' @param squeeze TODO
#' @param xlim TODO
#' @return TODO
#' @name gr.flatmap
#' @keywords internal
gr.flatmap = function(gr, windows, gap = 0, strand.agnostic = TRUE, squeeze = FALSE, xlim = c(0, 1))
{
  if (strand.agnostic)
    GenomicRanges::strand(windows) = "*"

  ## now flatten "window" coordinates, so we first map gr to windows
  ## (replicating some gr if necessary)
  #    h = findOverlaps(gr, windows)

  h = gr.findoverlaps(gr, windows);

  window.segs = gr.flatten(windows, gap = gap)

  grl.segs = BiocGenerics::as.data.frame(gr);
  grl.segs = grl.segs[values(h)$query.id, ];
  grl.segs$query.id = values(h)$query.id;
  grl.segs$window = values(h)$subject.id
  grl.segs$start = start(h);
  grl.segs$end = end(h);
  grl.segs$pos1 = pmax(window.segs[values(h)$subject.id, ]$start,
                       window.segs[values(h)$subject.id, ]$start + grl.segs$start - start(windows)[values(h)$subject.id])
  grl.segs$pos2 = pmin(window.segs[values(h)$subject.id, ]$end,
                       window.segs[values(h)$subject.id, ]$start + grl.segs$end - start(windows)[values(h)$subject.id])
  grl.segs$chr = grl.segs$seqnames

  if (squeeze)
  {
    min.win = min(window.segs$start)
    max.win = max(window.segs$end)
    grl.segs$pos1 = affine.map(grl.segs$pos1, xlim = c(min.win, max.win), ylim = xlim)
    grl.segs$pos2 = affine.map(grl.segs$pos2, xlim = c(min.win, max.win), ylim = xlim)
    window.segs$start = affine.map(window.segs$start, xlim = c(min.win, max.win), ylim = xlim)
    window.segs$end = affine.map(window.segs$end, xlim = c(min.win, max.win), ylim = xlim)
  }

  return(list(grl.segs = grl.segs, window.segs = window.segs))
}

#' @keywords internal
gr.stripstrand = function(gr)
{
  GenomicRanges::strand(gr) = "*"
  return(gr)
}

#' @keywords internal
gr.in = function(query, subject, ...)
{
  tmp = gr.findoverlaps(query, subject, ...)
  out = rep(FALSE, length(query))
  out[tmp$query.id] = TRUE

  return(out)
}

#' \code{stats::aggregate}, but returns vector
#'
#' @description
#' Same as \code{stats::aggregate} except returns named vector
#' with names as first column of output and values as second
#'
#' Note: there is no need to ever use aggregate or vaggregate, just switch to data.table
#'
#' @param ... arguments to aggregate
#' @return named vector indexed by levels of "by"
#' @author Marcin Imielinski
#' @keywords internal
vaggregate = function(...)
{
  out = aggregate(...);
  return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
}

# clip.seg
#
# (data frame seg function)
#
# Clips a pile of segments / ranges to a given region (a single seg) .. ie removing all segments outside of that region
# and then trimming the segments overlapping "region" to the region boundaries
#
clip.seg = function(x, region)
{
  x = standardize_segs(x);
  region = standardize_segs(region);
  if (nrow(region)>1)
    stop('Can only region to a single range')

  ## if either x or region do not specify chrom then we ignore chrom
  if (is.null(x$chr))
    region$chr = NULL;

  if (is.null(region$chr))
    x$chr = NULL;

  ix = which(seg.on.seg(x, region));
  x = x[ix, ]
  if (nrow(x)>0)
  {
    x$pos1 = pmax(region$pos1, x$pos1)
    x$pos2 = pmin(region$pos2, x$pos2)
  }
  return(x)
}

#' standardize_segs
#'
#' (data frame seg function)
#'
#' Takes and returns segs data frame standardized to a single format (ie $chr, $pos1, $pos2)
#'
#' if chr = TRUE will ensure "chr" prefix is added to chromossome(if does not exist)
#' @keywords internal
standardize_segs = function(seg, chr = FALSE)
{
  #if (inherits(seg, 'IRangesList'))
  #  seg = irl2gr(seg);

  if (is(seg, 'matrix'))
    seg = as.data.frame(seg, stringsAsFactors = FALSE)

  if (inherits(seg, 'RangedData') | inherits(seg, 'GRanges') | inherits(seg, 'IRanges'))
  {
    val = as.data.frame(values(seg));
    values(seg) = NULL;
    seg = as.data.frame(seg, row.names = NULL);  ## returns compressed iranges list
    seg$seqnames = as.character(seg$seqnames)
  }
  else
    val = NULL;

  field.aliases = list(
    ID = c('id', 'patient', 'Sample'),
    chr = c('seqnames', 'chrom', 'Chromosome', "contig", "seqnames", "seqname", "space", 'chr', 'Seqnames'),
    pos1 = c('start', 'loc.start', 'begin', 'Start', 'start', 'Start.bp', 'Start_position', 'pos', 'pos1', 'left', 's1'),
    pos2 =  c('end', 'loc.end', 'End', 'end', "stop", 'End.bp', 'End_position', 'pos2', 'right', 'e1'),
    strand = c('strand', 'str', 'strand', 'Strand', 'Str')
  );

  if (is.null(val))
    val = seg[, setdiff(names(seg), unlist(field.aliases))]

  seg = seg[, intersect(names(seg), unlist(field.aliases))]

  for (field in setdiff(names(field.aliases), names(seg)))
    if (!(field %in% names(seg)))
      names(seg)[names(seg) %in% field.aliases[[field]]] = field;

  if (chr)
    if (!is.null(seg$chr))
      if (!grepl('chr', seg$chr[1]))
        seg$chr = paste('chr', seg$chr, sep = "");

  if (is.null(seg$pos2))
    seg$pos2 = seg$pos1;

  missing.fields = setdiff(names(field.aliases), c(names(seg), c('chr', 'ID', 'strand')));

  if (length(missing.fields)>0)
    warning(sprintf('seg file format problem, missing an alias for the following fields:\n\t%s',
                    paste(sapply(missing.fields, function(x) paste(x, '(can also be', paste(field.aliases[[x]], collapse = ', '), ')')), collapse = "\n\t")));

  if (!is.null(val))
    seg = cbind(seg, val)

  return(seg)
}

# seg.on.seg
#
# (data frame version of GenomicRanges %in%)
#
# Given two collections of segs segs1 and segs2 (each a data fram of $chr, $start, $end) returns a logical vector of length nrow(segs1)
# that is TRUE in position i if segment i in segs1 overlaps one or more segments in segs2
#
# also works on data frames with the following segment nomenclatures, however data frame should only contain one of these:
# $chrom, $loc.start, $loc.end, $pos1, $pos2
#
# if $chr or $chrom not included then will only consider start and end
#
# e.g. given segments in data frame segs and snps in data frame map snp.in.seg(snps, segs) will flag snps that are in segment
#' @keywords internal
seg.on.seg = function(segs1, segs2, pad = 0, within = F # within = T -> output only T if segs1 range is within at least one segs2 range
)
{
  # convert col names to a standard format
  segs1 = standardize_segs(segs1);
  segs2 = standardize_segs(segs2);

  out = vector(mode = "logical", length = nrow(segs1))

  # if neither of the sets of segs have chr specified then we act as if all segs are on a single chromosome
  if (is.null(segs1$chr) & is.null(segs2$chr))
    segs1$chr = segs2$chr = 1;

  for (chr in unique(intersect(segs1$chr, segs2$chr)))
  {
    chrind1 = which(segs1$chr == chr)
    chrind2 = which(segs2$chr == chr)
    for (i in chrind2)
    {
      region = segs2[i,, drop = FALSE]
      if (within)
        out[chrind1] = out[chrind1] | (segs1$pos1[chrind1]>=(region$pos1-pad) & segs1$pos2[chrind1]<=(region$pos2+pad))   # seg1 is inside region
      else
        out[chrind1] = out[chrind1] | !(segs1$pos2[chrind1]<(region$pos1-pad) | segs1$pos1[chrind1]>(region$pos2+pad))
    }
  }

  return(out)
}

#' Deduplicate character vectors
#
#' Relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#' @param x Character vector to deduplicate.
#' @param suffix User defined suffix
#' @return dedupclicated character vector
#' @keywords internal
dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = unique(x[dup])
  udup.ix = lapply(udup, function(y) which(x==y));
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)
}

get_seqinfo <- function(.Object, seqinfo) {

  if (is.na(seqinfo) & is.list(.Object@data)) {

    slen = c()
    for (i in 1:length(.Object@data))
    {
      x = .Object@data[[i]]



      if (is(x, 'GRanges') || is(x, 'GRangesList'))
      {
        slen = GenomeInfoDb::seqlengths(x)

        if (any(is.na(slen)))
        {
          if (is(x, 'GRanges'))
            slen = GenomeInfoDb::seqlengths(gUtils::gr.fix(x))
          else if (is(x, 'GRangesList'))
            slen = GenomeInfoDb::seqlengths(gUtils::gr.fix(unlist(x)))
        }
      }
      else if (is(x, 'ffTrack'))
      {
        if (is.null(slen))
          slen = GenomeInfoDb::seqlengths(x)
        else
          slen[seqlevels(x)] = pmax(slen[seqlevels(x)], GenomeInfoDb::seqlengths(x), na.rm = TRUE)
        .Object[i, 'yaxis'] = TRUE
      }
      else if (inherits(x, 'RleList'))
      {
        .Object@data[[i]] = as(x, 'RleList')
        .Object[i, 'yaxis'] = TRUE
        slen = structure(sapply(x, length), names = names(x))
      }
      else if (is.character(x))
      {
        if (!file.exists(x))
          stop('External file does not exist.  Note the file must be a supported UCSC file format or R .rds file with the associated file extension .bw, .wig, .bed., .gff, .2bit, .bedgraph, .rds.')

        if (grepl('(\\.bw)|(\\.bigwig)', x, ignore.case = TRUE))
        {
          f = rtracklayer::BigWigFile(normalizePath(x))
          slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
          .Object[i, 'yaxis'] = TRUE
          if (is.null(slen))
            stop('External file must be a valid and existing .bigwig file')
        }
        else if (grepl('\\.wig', x, ignore.case = TRUE))
        {
          f = rtracklayer::WIGFile(normalizePath(x))
          slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
          .Object[i, 'yaxis'] = TRUE
          if (is.null(slen))
            stop('External file must be a valid and existing .wig file')
        }
        else if (grepl('\\.bed', x, ignore.case = TRUE)) ## only option is to load
        {
          f = rtracklayer::BEDFile(normalizePath(x))
          tmp.out = tryCatch(rtracklayer::import(f, asRangedData = FALSE), error = function(x) NULL)
          if (is.null(tmp.out))
            stop('External file must be a valid and existing .bed file')
          else
            .Object@data[[i]] = tmp.out

          slen = GenomeInfoDb::seqlengths(tmp.out)
          #                            slen = tryCatch(seqlengths(f), error = function(x) NULL)
          #                            if (is.null(slen))
          #                              stop('External file must be a valid and existing .bed file')
        }
        else if (grepl('\\.gff', x, ignore.case = TRUE))
        {
          f = rtracklayer::GFFFile(normalizePath(x))
          slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
          if (is.null(slen))
            stop('External file must be a valid and existing .gff file')
        }
        else if (grepl('\\.2bit', x, ignore.case = TRUE)) {
          f = rtracklayer::TwoBitFile(normalizePath(x))
          slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
          if (is.null(slen))
            stop('External file must be a valid and existing .2bit file')
        }
        else if (grepl('\\.bedgraph', x, ignore.case = TRUE)) ## only option is to load
        {
          f = rtracklayer::BedGraphFile(normalizePath(x))
          tmp.out = tryCatch(rtracklayer::import(f, asRangedData = FALSE), error = function(x) NULL)
          if (is.null(tmp.out))
            stop('External file must be a valid and existing .bedgraph file')
          else
            .Object@data[[i]] = tmp.out

          slen = GenomeInfoDb::seqlengths(tmp.out)
        }
        else if (grepl('\\.rds', x, ignore.case = TRUE)) ## assume this is GRanges or GRangesList
        {
          tmp.out = tryCatch(readRDS(x), error = function(x) NULL)

          if (is.null(tmp.out))
            stop('External file must be a valid and existing .rds file containing GRanges or GRangesList')
          else if (!inherits(tmp.out, 'GRanges') | !inherits(tmp.out, 'GRangesList'))
            stop('External file must be a valid and existing .rds file containing GRanges or GRangesList')
          else
            .Object@data[[i]] = tmp.out

          slen = GenomeInfoDb::seqlengths(tmp.out)
        }
      }
      else
        if (is.null(slen))
          slen= structure(sapply(x, length), names = names(x))
        else
          slen[names(x)] = pmax(slen[x], sapply(x, length), na.rm = TRUE)
    }

    if (any(.Object@formatting$chr.sub))
      names(slen) = gsub('chr', '', names(slen))

    if (any(duplicated(names(slen))))
      slen = data.table::data.table(len = slen, nm = names(slen))[, list(len = max(len)), by = nm][, structure(len, names = nm)]

    .Object@seqinfo = Seqinfo(seqnames = names(slen), seqlengths = slen);

  } else {
    if (is.null(seqinfo))
      .Object@seqinfo = Seqinfo()
    else if (is.list(seqinfo) & all(sapply(seqinfo, inherits, 'Seqinfo')) & length(seqinfo) == length(.Object))
      .Object@seqinfo = seqinfo
    else
    {
      if (!inherits(seqinfo, 'Seqinfo'))
        seqinfo = GenomeInfoDb::seqinfo(seqinfo);

      .Object@seqinfo = seqinfo;
      .Object@data = lapply(.Object@data, gUtils::gr.fix, seqinfo)
    }

  }
  return(.Object)
}

#' convert input data into a list of length 'len' of type FUN
#' @keywords internal
listify <- function(x, FUN, len = 1) {
  if (is.null(x))
    return(rep(list(FUN()), len))
  if (!is.list(x))
    return(list(x))
  return(x)
}



