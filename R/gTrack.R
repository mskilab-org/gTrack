#' @importFrom data.table data.table
#'
#' @section Slots:
#' \describe{
#'   \item{data}{length(.Object) length list containing genomic data (e.g. GRanges, GrangesLists, RleLists, path to ucsc file or ffTrack file on disk)}
#'   \item{seqinfo}{Seqinfo object}
#'   \item{colormap}{length(.Object) length named list of named vectors whose entry i specifies the colormap for the meta data field of object entry i (specified by the name) and maps unique value of that data field to colors (specified by the named vector)}
#'   \item{edges}{list of data.frames of length length(.Object) which has columns $from, $to, and optional fields $col, $lwd, and $lty to specify splined edges joining data items in the corresponding track}
#'   \item{formatting}{\code{data.frame} holding all of the formatting options}
#'   \item{mdata}{\code{matrix} holding interaction data between genomic loci (e.g. hiC)}
#'   \item{vars}{List of variants associated with track}
#' }
#'
#' @name gTrack-class
#' @rdname gTrack-class
#' @description
#' S4 class for \code{gTrack}
#'
#' Class \code{gTrack} defines a subsettable object that wraps formatting information around
#' genomic data (GRanges, GRangesList, RleLists, UCSC formats, ffTrack) and can be displayed
#' in a "genome browser" style plot.
#'
#' @exportClass gTrack
#' @author Marcin Imielinski
#' @importFrom methods setClass setGeneric setMethod setRefClass
#' @importFrom gUtils grl.unlist si2gr grbind gr.string gr.fix grl.pivot gr.findoverlaps gr.flatten gr.chr gr.match
setClass('gTrack', representation(data = 'list', mdata= 'list', seqinfo = 'Seqinfo', formatting = 'data.frame', colormap = 'list', edges = 'list', vars = 'list'))

#setClass('trackData', contains = "gTrack") ## for legacy, backwards compatibility with old trackData class


setMethod('initialize', 'gTrack', function(.Object, data, mdata, edges, vars, colormaps, seqinfo, y.field, yaxis, format, ...) ## only place NON formatting fields here. The rest passed with ...
    {
        .Object@data <- listify(data, GRanges)
  # if (is.null(data))
  #   .Object@data = list(GRanges())
  # else if (!is.list(data))
  #   .Object@data = list(data)
  # else
  #   .Object@data = data;

  if (any(ix <- !(sapply(.Object@data, class) %in% c('GRanges', 'GRangesList', 'character', 'RleList', 'ffTrack'))))
    stop('check input: gTrack objects can only be defined around GRanges, GRangesLists, RleLists, ffTrack, file paths to .rds files of the latter object types, or file paths to  UCSC format files')

  .Object@mdata <- listify(mdata, matrix, length(.Object@data))
  # if (is.null(mdata))
  #   .Object@mdata = rep(list(matrix()), length(.Object@data))
  # else if (is.matrix(mdata) | inherits(mdata, 'Matrix'))
  #   .Object@mdata = list(mdata)
  # else if (is.list(mdata))
  #   .Object@mdata = mdata
  # else
  #   stop('optional arg mdata must be a square matrix of same dimensions as data GRanges')

  for (i in 1:length(.Object))
    if (!is.matrix(.Object@mdata[[i]]) & !inherits(.Object@mdata[[i]], 'Matrix'))
      stop('optional arg mdata be either empty matrix or a square matrix of same dimensions as data GRanges')
  else if (!identical(.Object@mdata[[i]], matrix()))
  {
    if (!identical(dim(.Object@mdata[[i]]), rep(length(.Object@data[[i]]), 2)))
      stop('mdata for each entry must be a square matrix of same dimensions as data GRanges')
    .Object@mdata[[i]] = (t(.Object@mdata[[i]]) + .Object@mdata[[i]])/2
  }

#  .Object@edges <- listify(edges, data.frame, length(.Object@data)) ## listify not working here i.e. for concatenating objects with @edges field
  if (is.null(edges))
    .Object@edges = rep(list(data.frame()), length(.Object@data))
  else if (!is.list(edges) | inherits(edges, 'data.frame'))
    .Object@edges = list(edges)
  else
      .Object@edges = edges

  .Object@vars <- listify(vars, list, length(.Object@data))
  # if (is.null(vars))
  #   .Object@vars = rep(list(list()), length(.Object@data))
  # else if (!is.list(edgevars) | inherits(vars, 'data.frame'))
  #   .Object@vars = list(vars)
  # else
  #   .Object@vars = vars

  null.ix = sapply(.Object@data, is.null)
  .Object@data[null.ix] = GRanges()

  edges.df = sapply(.Object@edges, is.data.frame)
  edges.mat = sapply(.Object@edges, function(x) is.array(x) | inherits(x, 'Matrix'))

  if (any(edges.df))
      .Object@edges[edges.df] = lapply(.Object@edges[edges.df], function(x)
    {
      if (nrow(x)>0)
      {
        if (!all(c('from', 'to') %in% names(x)))
          stop('edges data frame missing $to and $from columns')

        if (data.table::is.data.table(x))
            x = as.data.frame(x)
        x
      }
      else
        data.frame()
  })
 
  if (any(edges.mat))
    .Object@edges[edges.mat] = lapply(which(edges.mat), function(x)
    {
      y = .Object@edges[[x]]
      if (nrow(y) != length(.Object@data[[x]]) & ncol(y) != length(.Object@data[[x]]))
        warning('dimensions of matrix do not match length of track data item.  Adjacency matrix of edges should be square and have as many rows as there are input ranges in the corresponding gTrack item')
      tmp.edges = which(y>0, arr.ind = TRUE)
      return(data.frame(from = tmp.edges[,1], to = tmp.edges[,2], lwd = y[tmp.edges]))
    })
  
  .Object@colormap <- listify(colormaps, list)
  # if (is.null(colormaps))
  #   .Object@colormap = list(list())
  # else if (!is.list(colormaps))
  #   .Object@colormap = list(colormaps)
  # else
  #   .Object@colormap = colormaps;

  if (length(.Object@colormap)==1)
    .Object@colormap = rep(.Object@colormap, length(.Object@data))

  if (is.null(names(.Object@colormap)))
    names(.Object@colormap) = NA

  ## convert options to formatting string
  if (!is.data.frame(format) && is.na(format))
    .Object@formatting = as.data.frame(list(...), stringsAsFactors=FALSE)
  else if (is.data.frame(format) && nrow(format) == length(.Object))
    .Object@formatting <- format #formatting(.Object) <- format
  else
    stop("Expecting formatting (as supplied directly to constructor) to be data.frame of nrow = length(gTrack)")

  # if (length(setdiff(c(length(height), length(col), length(lift), length(format)), c(1, length(.Object@data))))>0)
  #   stop('Length of formatting attributes (height, lift, format, col) have to be either 1 or the gTrack object length');

  ## figure out seqinfo from provided tracks, make sure that all the seqinfos are compatible
  ## if only Rle's provided

  ## JEREMIAH 3/20/16 t.name does not need default
  # if (is.na(name))
  #   t.name = names(.Object@data)
  # else
  #   t.name = name
  # if (is.null(t.name))
  #   if (!all(is.na(y.field)))
  #     t.name = y.field
  #
  # if (is.null(t.name))
  #   t.name = rep(NA, length(.Object@data))


  # .Object@formatting = data.frame(name = t.name, height = height, ygap = ygap, stack.gap = stack.gap,
  #                                 lift = lift, split = split, angle = angle, format = format, col = col,
  #                                 lwd.border = lwd.border, ypad = ypad,  ywid = ywid, border = border,
  #                                 hadj.label = hadj.label, gr.colorfield = gr.colorfield, smooth = smooth,
  #                                 round = round, vadj.label = vadj.label, y.field = y.field,
  #                                 circles = circles, bars = bars, y0.bar = y0.bar, lines = lines,
  #                                 source.file.chrsub = source.file.chrsub, yaxis.cex = yaxis.cex,
  #                                 max.ranges = max.ranges, yaxis = yaxis,
  #                                 yaxis.pretty = yaxis.pretty, y.quantile = y.quantile, triangle = triangle,
  #                                 is.null = is.null, stringsAsFactors = FALSE)

  # .Object@formatting$y0 = NA
  # .Object@formatting$y1 = NA

  # other.args = list(...)
  #
  # avail.args <- c("")
  # if (any(!names(other.args) %in% avail.args))
  #
  #   if (length(other.args)>0)
  #     for (i in 1:length(other.args))
  #       .Object@formatting[, names(other.args)[i]] = other.args[[i]]
  #
  # if (!is.null(.Object@formatting$label))
  #   .Object@formatting$label = NULL

  ## JEREMIAH DO THIS LATER, no need to track now. 3/20/2016. Confusing to bc of formatting issue
  ## turn off the colors if a colormap is supplied
  # cix = sapply(.Object@colormap, length)>0 | !is.na(.Object@formatting$gr.colorfield)
  # if (any(cix))
  #   .Object@formatting$col[cix] = NA
  # ## set the default color
  # if (any(nix <- !cix & is.na(.Object@formatting$col)))
  #   .Object@formatting$col[nix] = alpha('black', 0.5)
  # ## set the default border
  # if (any(nix <- is.na(.Object@formatting$border)))
  #   .Object@formatting$border[nix] = .Object@formatting$col[nix]

  ## set the seqinfo if not provided
  .Object <- get_seqinfo(.Object, seqinfo)

  ## replicate if length y.field > 1 and length dat is 1
  #if (length(y.field)>1 & length(.Object)==1)
  #  .Object@data = rep(.Object@data, length(y.field))

  ## if single data input but multiple yfields, replicate

  if (length(y.field) > 1 && length(.Object@data) == 1) {
    new_object = .Object
    new_object$y.field = y.field[1]
    new_object$yaxis = yaxis[1]
    for (i in 2:length(y.field)) {
      o = .Object
      o$y.field = y.field[i]
      o$yaxis = yaxis[i]
      new_object = c(new_object, o)
    }
    .Object = new_object
  } else if (!all(is.na(y.field))) {
    formatting(.Object)$y.field <- y.field
    formatting(.Object)$yaxis <- yaxis
  }
  else {
    formatting(.Object)$y.field <- NA
    formatting(.Object)$yaxis <- TRUE
}
  ##
  if (any(!is.na(y.field))) {
      ix <- nchar(formatting(.Object)$name) == 0
      formatting(.Object)$name[ix] <- y.field[ix]
  }
  
  return(.Object)  
})




#' Construct a new \code{gTrack}
#'
#' Arguments described as "formatting" are vectors.  They are replicated (if necessary) to match the length of the object.
#' @param data Instance of one of the following objects: (1) GRanges (2) GRangesList (3) ffTrack or
#' (4) character representing Bed, Wig, BigWig, or .rds GRanges file.  It can also be a list of the above (if specifying multiplying
#'  tracks.
#' @param mdata Square matrix representing values at connections between genomic elements (e.g. chr8:1,000-2,000 to chr3:5,000-6,000). The number of
#' rows (and thus columns) must be equal to the length of the \code{data} element. The value of the (i,j) element of this matrix represents a value
#' for the intersection of \code{data[i]} and \code{data[j]}. Note that all connections are assumed to be symmetric, and only the upper triangle
#' will be read. Displays a \code{warning} if the matrix is non-symmetric
#' @param y.field vector or scalar character referencing meta data field of GRanges or GRangesList to use for "y axis" coordinate when plotting
#' numeric tracks, (note: for a RleList or ffTrack this field is automatically set to "score"), default is NA for non-numeric tracks
#' @param name vector or scalar character specifying name of this track, which will be displayed on label to the left of the track
#' @param height vector or scalar numeric specifying height of track(s) (in relative units)
#' @param gr.labelfield vector or scalar character specifying which GRanges meta data field to use for GRanges label (default "label") (formatting)
#' @param grl.labelfield vector or scalar character specifying which GRanges meta data field to use for GRangesList label (default "label") (formatting)
#' @param legend.maxitems Scalar positive integer specifying what is the maximum number of items to include in legend [Inf]
#' @param label.suppress vector or scalar logical flag specifying whether to suppress all GRanges / GRangesList label drawing  (formatting)
#' @param label.suppress.gr vector or scalar logical flag specifying whether to suppress GRanges label drawing  (formatting)
#' @param label.suppress.grl vector or scalar logical flag specifying whether to suppress GRangesList label drawing  (formatting)
#' @param ygap vector or scalar numeric specifying gap between tracks
#' @param stack.gap vector or scalar.numeric specifying x gap between stacking non-numeric GRanges or GrangesLists items in track(s)
#' @param cex.label size of top level label (e.g. GRAngesList item)
#' @param cex.label.gr size of bottom level label (i.e. GRanges in a GrangesListitem, e.g. exons in a gene)
#' @param gr.srt.label rotation of GRanges label
#' @param ywid vector or scalar numeric specifying the y-extent of individual ranges (in local plot coordinates)
#' @param col vector or scalar character specifying static color for track(s), if NA then color is specified by colormaps() property or gr.colorfield or col meta data field of GRanges / GRangesList data object
#' @param border vector or scalar character specifying static border for polygons in track(s), if NA then $border is determine duing gr.colorfield / colormap or meta field $border of GRanges / GRangesList
#' @param max.ranges vector or scalar numeric specifying what is the max number of ranges to plot in a window (formatting)
#' @param angle vector of scalar numeric specifying angle of polygons that represent signed range'
#' items (only relevant if used within chainedTracks object)
#' @param split vector or scalar logical flag specifying whether to split when lifting (only relevant if used wihtin chainedTracks object)
#' @param colormaps length(.Object) length named list of named vectors whose entry i maps uniques value of a data field to colors.  The data.field is specified by the name of the list entry, and the unique values / colors are specified by the named vector.
#' @param edges Data frame of columns $from, $to, and optional fields $col, $lwd, and $lty, specifying edges linking data items.
#' Also can be a list of the above if specifying multiple tracks (and must be compatible in length with data arg)
#' @param xaxis.suffix vector or scalar numeric specifying the suffix that will be used in describin x axis coordinates (TODO: move to display) (formatting)
#' @param xaxis.unit vector or scalar numeric specifying the unit that will be used in describing x axis coordinates (TODO: move to display)  (formatting)
#' @param xaxis.round vector or scalar non-neg integer specifying number of decimals to round xaxis coordinate labels (formatting)
#' @param xaxis.nticks vector or scalar positive integer specifying how many xaxis ticks to optimally draw (formatting)
#' @param xaxis.label.angle vector or scalar numeric between 0 and 360 specifying angle with which to draw xaxis coordinate labels (formatting)
#' @param xaxis.newline vector or scalar logical specifying whether to draw a newline in the xaxis coordinate labels (formatting)
#' @param xaxis.width Logical scalar specifying whether to add window width to xaxis window labels [TRUE]
#' @param lwd.border vector or scalar integer specifying the thickness of the polygon borders (formatting)
#' @param cex.label vector or scalar numeric specifying the expansion factor of the range labels (formatting)
#' @param hadj.label vector or scalar numeric specifying the horizontal adjustment of the range labels (formatting)
#' @param vadj.label vector or scalar numeric specifying the vertical adjustment of the range labels (formatting)
#' @param ypad vector or scalar numeric between 0 and 1 specifying how much whitespace padding to add within panel (formatting)
#' @param circles vector or scalar logical specifying whether to scatter plot range data (formatting)
#' @param lines vector or scalar logical specifying whether to line plot range data (formatting)
#' @param bars vector or scalar logical specifying whether to bar plot range data (formatting)
#' @param y0.bar vector or scalar numeric specifying where to draw the lower boundary of a bar in a bar plot (only applicable if bars == T) (formatting)
#' @param source.file.chrsub vector or scalar logical specifying whether or not sub "chr" out of any external files (e.g. UCSC style files) (formatting)
#' @param y.grid.col vector or scalar character specifying color of "gridlines" used to specify numeric track data (formatting)
#' @param y.grid.cex vector or scalar non-neg numeric specifying character expansion for y tick / y grid labels (formatting)
#' @param y.grid.lty vector or scalar positive integer specifying line style of y grid lines for numeric tracks (formatting)
#' @param y.grid.lwd vector or scalar positive integer specifying thickness of y grid lines for numeric tracks (formatting)
#' @param y.grid.labx vector or scalar positive integer specifying fraction of xlim left of plot to place y axis labels (formatting)
#' @param yaxis vector or scalar logical specifying whether to print yaxis (formatting)
#' @param yaxis.pretty vector or scalar positive integer specifying how many ticks to optimally draw on yaxis (formatting)
#' @param draw.var vector or scalar logical specifying whether to draw.var for GRanges / GRangesList  specifying reads
#'  (GRanges must contain $cigar +/- $MD field) (formatting)
#' @param draw.paths vector or scalar logical specifying whether to interpret GRangesLists as "paths" and connect them with a
#'  a set of spline curves.  (formatting)
#' @param path.col vector or scalar character specifying color of path (only applicable for tracks in which draw.paths = T) (formatting)
#' @param path.col.arrow vector or scalar character specifying color of arrow of path (only applicable for tracks in which draw.paths = T) (formatting)
#' @param path.cex.arrow vector or scalar numeric > 0 specifying expansion factor of arrow of path (only applicable for tracks
#' in which draw.paths = T)  (formatting)
#' @param path.stack.y.gap vector or scalar numeric > 0 specifying y stack gap of paths (only applicable for tracks in which draw.paths = T)  (formatting)
#' @param path.stack.x.gap vector or scalar numeric > 0 specifying x stack gap for paths  (only applicable for tracks in which draw.paths = T)  (formatting)
#' @param path.cex.v vector or scalar numeric > 0 specifying vertical bulge of sline in paths (only applicable for tracks in which draw.paths = T)  (formatting)
#' @param path.cex.h vector or scalar numeric > 0 specifying horizontal bulge of spline in paths (formatting)
#' (only applicable for tracks in which draw.paths = T)  (formatting)
#' @param draw.backbone vector or scalar logical specifying whether to draw "backbone" connecting different items in a GRangesList item (formatting)
#' @param xaxis.cex.tick Scalar numeric specifying expansion factor for axis tick labels  (formatting)
#' @param xaxis.ticklen Scalar numeric specifying lengths for axis ticks  (formatting)
#' @param gr.cex.label vector or scalar numeric specifying GRanges label character expansion (default is cex.label)  (formatting)
#' @param gr.srt.label vector or scalar numeric between 0 and 180 specifying rotation of GRanges labels  (formatting)
#' @param sep.lty Scalar integer specifying line style for window separators [2] (dashed line)
#' @param sep.lwd Scalar numeric specifying line thickness for window separators [1]
#' @param sep.bg.col Color (supplied by name or hex code) for the background of the windows [gray95]
#' @param sep.draw Logical allowing separators to be turned off [FALSE]
#' @param cmap.min minimum saturating data value of color map for triangle plot
#' @param cmap.max maximum saturating data value of color map for triangle plot
#' @param labels.suppress whether to suppress labels (both GRangesList and GRanges)
#' @param labels.suppress.gr whether to suppress GRanges labels
#' @param labels.suppress.gr whether to suppress GRangesList labels
#' @param ... additional formatting arguments to gTrack
#'
#' @name gTrack
#' @rdname gTrack-class
#' @export
#' @author Marcin Imielinski
gTrack = function(data = NULL, ##
                  y.field = NA, ## character field specifying what numeric values field of the underlying gr or grl should be used as y coordinate (otherwise ranges are stacked)
                  mdata = NULL, ## this is matrix of values for triangle plot
                  edges = NULL, ## this is a list of data.frames of length numtracks or a data.frame, data.frame has fields $from, $to, and optional fields $col, $lwd, $lty to specify edges joining ranges and their display, can also be an adjacency matrix of dim n x n where n =  length(GRanges) or length(GRangesList), the values of the matrix will specify $lwd of the resulting lines
                  vars = NULL, ## GRangesList of same length as data (if length is 1 ) or list of GRangesList representing variants (output of draw.var)
                  colormaps = NULL, # (named) list same length as data
                  height = 10,
                  ygap = 2,
                  stack.gap = 0,
                  cex.label = 1,
                  gr.cex.label = cex.label *0.8,
                  gr.srt.label = 0,
                  col = NA,
                  border = NA,
                  angle = 15,
                  name = "",
                  gr.colorfield = NA,
                  y.quantile = 0.01, ## if y0 or y1 is not specified then will draw between y.quantile and 1-y.quantile of data
                  y.cap = T, ## whether to cap values at y0, y1 (only relevant if y.field specified)
                  lwd.border = 1,
                  hadj.label = 1,
                  vadj.label = 0.5,
                  smooth = NA, ## smooth with running mean with this window
                  round = NA, ## round the output of running mean to this precision
                  ywid = NA,
                  ypad = 0,
                  seqinfo = NA,
                  circles = FALSE,
                  lines = FALSE,
                  bars = FALSE,
                  draw.paths = FALSE,
                  draw.var = FALSE, 
                  triangle = !is.null(mdata),
                  max.ranges = 5e4, ## parameter to limit max number of ranges to draw on canvas, will downsample to this amount
                  source.file.chrsub = T, ## if source file has chr for seqnames this will sub it out
                  y0.bar = NA,
                  yaxis = !is.na(y.field), # logical whether to print yaxis
                  yaxis.pretty = 5, # how many ticks to optimally draw on yaxis
                  yaxis.cex = 1, ## size of text at yaxis
                  chr.sub = TRUE, ## remove 'chr' from slen if drawing from UCSC style format
                  edgevars = NA,
                  gr.labelfield = NA,
                  grl.labelfield = NA,
                  xaxis.prefix = "",
                  xaxis.unit = 1,
                  xaxis.suffix = "",
                  xaxis.round = 3,
                  xaxis.cex.label = 1,
                  xaxis.newline = TRUE,
                  xaxis.chronly = FALSE,
                  xaxis.width= TRUE,
                  xaxis.interval = 'auto',
                  xaxis.label.angle = 0,
                  xaxis.ticklen = 1,
                  xaxis.cex.tick = 1,
                  sep.lty = 2,
                  sep.lwd = 1,
                  sep.bg.col = 'gray95',
                  sep.draw = TRUE,
                  y0 = NA,
                  y1 = NA,
                  m.sep.lwd = 1,
                  m.bg.col = 'white',
                  cmap.min = NA,
                  cmap.max = NA,
    labels.suppress = FALSE,
    labels.suppress.grl = labels.suppress,
    labels.suppress.gr = labels.suppress,
                  bg.col = 'white', ## background color of whole thing
    formatting = NA) {

    ## make border batch color, but be a little less transparent
    if (is.na(col) & is.na(gr.colorfield) & is.null(colormaps))
        col = alpha('black', 0.5)
    
    if (is.na(border) & is.na(gr.colorfield) & is.null(colormaps))
        {
            rgb = col2rgb(col, alpha = TRUE)            
            border = rgb(rgb['red', ]/255, rgb['green', ]/255, rgb['blue', ]/255, alpha = 0.9)
        }


    
    ## TODO: FIX THIS USING formals() and some eval / do.call syntax or something similar 
    new('gTrack', data = data, y.field = y.field, mdata = mdata, name = name, format = formatting,
      edges = edges, vars = vars, draw.paths = draw.paths, colormaps = colormaps, height = height, ygap = ygap,
      stack.gap = stack.gap, col = col, border = border, angle = angle, draw.var = draw.var,
      gr.colorfield = gr.colorfield, y.quantile = y.quantile,
      cex.label = cex.label, gr.cex.label.gr = gr.cex.label, gr.srt.label = gr.srt.label,
      y.cap = y.cap, lwd.border = lwd.border, hadj.label = hadj.label, vadj.label = vadj.label, smooth = smooth,
      round = round, ywid = ywid, ypad = ypad, seqinfo = seqinfo, circles = circles, lines = lines,
      bars = bars, triangle = triangle, max.ranges = max.ranges, source.file.chrsub = source.file.chrsub,
      y0.bar = y0.bar, yaxis = yaxis, yaxis.pretty = yaxis.pretty, yaxis.cex = yaxis.cex,
      chr.sub = chr.sub, edgevars = edgevars, gr.labelfield = gr.labelfield,
      grl.labelfield = grl.labelfield, xaxis.prefix = xaxis.prefix, xaxis.unit = xaxis.unit,
      xaxis.suffix = xaxis.suffix, xaxis.round = xaxis.round, xaxis.interval = xaxis.interval,
      xaxis.cex.label = xaxis.cex.label, xaxis.newline = xaxis.newline,
        xaxis.chronly = xaxis.chronly, xaxis.width = xaxis.width,
        labels.suppress = labels.suppress, labels.suppress.gr = labels.suppress.gr, labels.suppress.grl = labels.suppress.grl,
      xaxis.label.angle = xaxis.label.angle, xaxis.ticklen = xaxis.ticklen,
      xaxis.cex.tick = xaxis.cex.tick, sep.lty = sep.lty, sep.lwd = sep.lwd, sep.bg.col = sep.bg.col,
      sep.draw = sep.draw, y0 = y0, y1 = y1, m.sep.lwd = m.sep.lwd, m.bg.col = m.bg.col,
      cmap.min = cmap.min, cmap.max = cmap.max, bg.col = bg.col)
}


setValidity('gTrack', function(object)
{
  problems = c();

  ALLOWED.DATA.CLASSES = c('GRangesList', 'GRanges', 'RleList', 'SimpleRleList', 'character', 'ffTrack');
  ## deprecated a while ago... 3/20/16


  ##ALLOWED.FORMATS = list(RleList = c('scatter', 'line', 'bar', 'ranges'), SimpleRleList = c('scatter', 'line', 'bar', 'ranges'), GRanges = c('ranges'), GRangesList = c('ranges'), character = c('ranges'), ffTrack = 'ranges');

  if (length(object@data)>0)
    if (!all(sapply(object@data, function(x) class(x) %in% ALLOWED.DATA.CLASSES)))
      problems = c(problems, (paste('One or more of the provided data objects in data slot are not of the allowed data classes', paste(ALLOWED.DATA.CLASSES, collapse = ", "))))
    else
    {
      # for (i in 1:nrow(object@formatting))
      #   if (!(object@formatting$format[i] %in% ALLOWED.FORMATS[[class(object@data[[i]])]]))
      #     problems = c(problems, sprintf('Format specified for data track %s (%s) is not allowed.  Available formats for this track are: %s.', i, object@formatting$format[i], paste(ALLOWED.FORMATS[[class(object@data[[i]])]], collapse = ", ")))
      if (length(object@data) != nrow(object@formatting))
        problems = c(problems, 'Length of data is not equal to the number of formatting entries')
    }
    else
      if (nrow(object@formatting) != 0)
        problems = c(problems, 'Null trackdata object has incompatible fields');


      if (!all(sapply(object@edges, is.data.frame)))
        problems = c(problems, 'Some trackdata edges attributes are not data.frames')
      else if (any(!sapply(object@edges, function(x) if (nrow(x)>0) all(c('from', 'to') %in% colnames(x)) else T)))
        problems = c(problems, 'Some nonempty trackdata edges attributes are missing $to and $from fields')
      else if (any(!sapply(1:length(object@data), function(x)
        if (nrow(object@edges[[x]])>0)
          all(object@edges[[x]]$from <= length(object@data[[x]])) & all(object@edges[[x]]$to <= length(object@data[[x]]))
        else T)))
        problems = c(problems, 'Some nonempty trackdata edges $to and $from fields are out of bounds (ie exceed the length of the data field of the corresponding gTrack item')

      if (!is.null(formatting(object)$y.field) && !is.na(formatting(object)$y.field))
      {
        nix = !is.na(object$y.field) & sapply(dat(object), inherits, 'GRanges')
        if (any(nix))
          tmp = sapply(which(nix), function(x) object$y.field[x] %in% names(values(dat(object)[[x]])))
        else
          tmp = TRUE

        if (any(!tmp))
        {
          problems = c(problems, paste('These y.fields are not found in their respective GRanges object:',
                                       paste(object$y.field[nix][!tmp], collapse = ', ')))
        }
      }


      if (length(object@data) != length(object@colormap))
        problems = c(problems, 'Length of object is not the same length as colormap')

      if (any(ix <- sapply(object@data, class) == 'character'))
      {
        for (iix in which(ix))
        {
          x = object@data[[iix]]

          if (grepl('(\\.bw$)|(\\.bigwig$)', x, ignore.case = T))
          {
            f = BigWigFile(normalizePath(x))
            slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
            if (is.null(slen))
              problems = c(problems, sprintf('External file %s is not a valid and existing .bw / .bigwig file', x))
          }
          else if (grepl('\\.wig', x, ignore.case = T))
          {
            f = WIGFile(normalizePath(x))
            slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
            if (is.null(slen))
              problems = c(problems, sprintf('External file %s is not a valid and existing .wig file', x))

          }
          else if (grepl('\\.bed', x, ignore.case = T))
          {
            f = BEDFile(normalizePath(x))
            #                        slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
            #                        if (is.null(slen))
            #                          problems = c(problems, sprintf('External file %s is not a valid and existing .bed file', x))
          }
          else if (grepl('\\.gff', x, ignore.case = T))
          {
            f = GFFFile(normalizePath(x))
            slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
            if (is.null(slen))
              problems = c(problems, sprintf('External file %s is not a valid and existing .gff file', x))
          }
          else if (grepl('\\.2bit', x, ignore.case = T))
          {
            f = TwoBitFile(normalizePath(x))
            slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
            if (is.null(slen))
              problems = c(problems, sprintf('External file %s is not a valid and existing .2bit file', x))

          }
          else if (grepl('\\.bedgraph', x, ignore.case = T))
          {
            f = BEDGraphFile(normalizePath(x))
            #                       slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
            #                       if (is.null(slen))
            #                         problems = c(problems, sprintf('External file %s is not a valid and existing .bedgraph file', x))
          }
          else
            problems = c(problems, sprintf('External file %s does not map to a supported UCSC format or .rds. Supported files must have one of the following extensions: .bw, .bed. .bedgraph, .2bit, .wig, .gff, .rds.', x))
        }
      }

      if (length(problems)==0)
        TRUE
      else
        problems
}
);

#suppressWarnings(removeMethod('[', 'gTrack')) ## takes care of stupid R 2.15 bug


#' @name [
#' @title [
#' @description
#'
#' Subsetting tracks of gTrack gt
#' gt[1]
#' gt[ix] # where ix is an integer vector
#'
#' @param x \code{gTrack} object
#' @param i Integer specifying which \code{gTrack} to grab
#' @docType methods
#' @export
#' @aliases [,gTrack,ANY,ANY,ANY-method
#' @rdname sub-methods
#' @author Marcin Imielinski
setMethod('[', 'gTrack', function(x, i)
{
  if (is.logical(i))
    i = which(i)
  x@data = x@data[i]
  x@seqinfo = x@seqinfo
  #            x@height = x@height[i]
  x@formatting = x@formatting[i, ]
  x@colormap = x@colormap[i]
  x@edges = x@edges[i]

  if (.hasSlot(x, 'mdata'))
    x@mdata = x@mdata[i]
  else
    x@mdata = rep(list(matrix()), length(x))

  return(x)
})

#' @name mdata
#' @title mdata
#' @description
#'
#' Accessing submatrices of mdata object for length 1 gTrack gt
#'
#' Usage: (igr, jgr are GRanges objects corresponding to slices of matrix to be accessed from gt)
#' mdata(gt, igr, jgr)
#' @param x \code{gTrack} object containing a matrix in the \code{mdata} slot
#' @param igr \code{GRegion} object specifying positions to subset matrix with
#' @param jgr \code{GRegion} object specifying positions to subset matrix with
#' @author Jeremiah Wala
#' @docType methods
#' @rdname mdata-methods
#' @export
setGeneric('mdata', function(x, igr=NULL, jgr=NULL) standardGeneric('mdata'))


#' @rdname mdata-methods
#' @aliases mdata,gTrack,ANY,ANY,ANY-method
setMethod('mdata', signature=c("gTrack", "ANY", "ANY"), function(x, igr = NULL, jgr = igr)
{
  if (is.null(igr) & is.null(jgr))
    return(x@mdata)

  if (is.null(igr))
    igr = si2gr(x)

  if (is.null(jgr))
    jgr = si2gr(x)

  if (is(igr, 'Rle'))
    igr = as.character(igr)

  if (is(jgr, 'Rle'))
    jgr = as.character(jgr)

  if (is.character(igr))
    igr = si2gr(x)[igr]

  if (is.character(jgr))
    jgr = si2gr(x)[jgr]

  if (length(igr)==0 | length(jgr)==0)
    return(NULL)

  out = lapply(1:length(x), function(y)
  {
    if (is.null(x@mdata[[y]]))
      return(NULL)
    i = gUtils::gr.in(x@data[[y]], igr)
    j = gUtils::gr.in(x@data[[y]], jgr)
    rown = dedup(gr.string(x@data[[y]][i]))
    coln = dedup(gr.string(x@data[[y]][j]))
    tmp = (x@mdata[[y]][i, j] + t(x@mdata[[y]][j, i]))/2
    rownames(tmp) = rown
    colnames(tmp) = coln
    return(tmp)
  })
  if (length(x)==1)
    out = out[[1]]
  return(out)
})

#' @name $
#' @title $
#' @description
#'
#' Accessing columns of gTrack formatting data.frame
#'
#' @param x \code{gTrack} object
#' @param name Name of the \code{formatting} field to view
#' @docType methods
#' @rdname cash-methods
#' @aliases $,gTrack-method
#' @export
#' @author Marcin Imielinski
setMethod('$', 'gTrack', function(x, name)
{
  return(x@formatting[, name])
})

#' @name $<-
#' @title $<-
#' @description
#'
#' Setting formats of gTrack object - ie modifying the formatting(gt) data frame after
#' an object has already been instantiated
#'
#' Usage:
#' gt$y.field = 'score'
#' gt$gr.colorfield[1] = 'readtype'
#'
#' @param x \code{gTrack} object to alter \code{formatting} field of
#' @param name \code{formatting} field to alter
#' @param value New value
#' @docType methods
#' @rdname cash-set-methods
#' @aliases $<-,gTrack-method
#' @export
#' @author Marcin Imielinski
setMethod('$<-', 'gTrack', function(x, name, value)
{
  x@formatting[, name] = value
  return(x)
})

#' @name length
#' @title length
#' @description
#'
#' Getting length of gTrack object gt
#' Usage:
#' length(gt)
#' @param x \code{gTrack} object
#' @docType methods
#' @rdname length-methods
#' @aliases length,gTrack-method
#' @export
#' @author Marcin Imielinski
setMethod('length', 'gTrack', function(x)
{
  return(length(x@data))
})

#' @name reduce
#' @title reduce
#' @description
#'
#' Computing the GRanges footprint of the gTrack object on the genome
#' usage:
#' reduce(gt) # outputs a GRanges
#' @param x \code{gTrack} object to retrieve reduced \code{GRanges} from
#' @param ... additional arguments to GRanges reduce function
#' @return \code{GRanges} with the minimal footprint of the \code{gTrack} data
#' @importFrom GenomicRanges reduce
#' @importFrom gUtils gr.sub
#' @docType methods
#' @rdname reduce-methods
#' @aliases reduce,gTrack-method
#' @export
#' @author Marcin Imielinski
setMethod('reduce', 'gTrack', function(x, ... )
{
  if (length(x)==0)
    return(GRanges(seqlengths = GenomeInfoDb::seqlengths(x)))

  .dat2gr = function(y)
  {
    if (is(y, 'Rle'))
    {
      y[is.na(y)] = 0;
      out = as(y!=0, 'GRanges');
      return(out[out$score])
    }
    else if (is(y, 'GRangesList'))
      gr.stripstrand(unlist(y))
    else
      gr.stripstrand(y)
  }

  if (length(x)==1)
    return(reduce(.dat2gr(x@data[[1]]), ...))
  else
    return(reduce(do.call('grbind', lapply(x@data, .dat2gr)), ... ))
})

#uppressWarnings(removeMethod('seqinfo', 'gTrack')) ## takes care of stupid R 2.15 bug


#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' returns Seqinfo of gTrack object gt
#' Usage:
#' seqinfo(gt)
#' @docType methods
#' @param x \code{gTrack} object
#' @return \code{seqinfo}
#' @importFrom GenomicRanges seqinfo
#' @export
#' @author Marcin Imielinski
setMethod("seqinfo", signature(x = "gTrack"), function(x)
{
  return(x@seqinfo)
})

#' @name seqinfo<-
#' @title seqinfo,-
#' @description
#'
#' set seqinfo property of gTrack
#'
#' @export
#' @param .Object \code{gTrack} object
#' @param value new \code{Seqinfo} object
#' @docType methods
#' @rdname seqinfo-set-methods
#' @author Marcin Imielinski
setGeneric('seqinfo<-', function(.Object, value) standardGeneric('seqinfo<-'))


#' @rdname seqinfo-set-methods
#' @aliases seqinfo,gTrack-method
setReplaceMethod('seqinfo', 'gTrack', function(.Object, value)
{
  .Object@seqinfo = value;
  .Object@data = lapply(dat(.Object), gr.fix, value)
  validObject(.Object)
  return(.Object)
});

#' @name c
#' @title c
#' @description
#'
#' Concatenate gTrack objects gt1, gt2, gt3
#' c(gt1, gt2, gt3) # returns a gTrack object with the component tracks "stacked"
#'
#' @param x Initial \code{gTrack} object
#' @param ... Any number of \code{gTrack} objects
#' @param recursive If recursive = TRUE, the function recursively descends through lists (and pairlists) combining all their elements into a vector [FALSE]
#' @docType methods
#' @rdname c-methods
#' @aliases c,gTrack-method
#' @export
#' @author Marcin Imielinski
setMethod('c', 'gTrack', function(x, ..., recursive = FALSE)
{
  args = list(x, ...)

  if (any(ix <- sapply(args, inherits, 'trackData')))
    args[ix] = lapply(args[ix], function(y) {class(y) = 'gTrack'; return(y)})

  if (any(!sapply(args, inherits, 'gTrack')))
    stop('some objects in gTrack concatenation are not gTrack')

  has.mdata = sapply(args, .hasSlot, 'mdata')

  if (any(!has.mdata))
    for (j in which(!has.mdata))
      args[[j]]@mdata = rep(list(matrix()), length(args[[j]]))

  args = args[!sapply(args, is.null)]

  out <- gTrack(data = do.call('c', lapply(1:length(args), function(y) args[[y]]@data)),
                #                             seqinfo = do.call('c', lapply(1:length(args), function(y) args[[y]]@seqinfo)),
                colormaps = do.call('c', lapply(args, function(y) y@colormap)),
                edges = do.call('c', lapply(args, function(y) y@edges)),
                mdata = do.call('c', lapply(1:length(args), function(y) args[[y]]@mdata)),
                format = do.call('rrbind', lapply(args, formatting)),
                y.field = do.call('c', lapply(args, function(y) formatting(y)$y.field)),
                yaxis = do.call('c', lapply(args, function(y) formatting(y)$yaxis)))
  #out@mdata = do.call('c', lapply(1:length(args), function(y) args[[y]]@mdata))

##  formatting(out) <- do.call('rrbind', lapply(args, formatting))

  if ('is.null' %in% names(formatting(out)))
  {
    is.n = out$is.null
    is.n[is.na(is.n)] = F
    out = out[!is.n]
}
  return(out)
})

#' @name formatting
#' @title formatting
#' @description
#'
#' Get data frame specifying formatting of gTrack object gt
#' usage:
#' formatting(gt)
#' gt # will just display the formatting data.frame
#'
#' If you want to access particular fields of the formatting data.frame, just use the "$" accessor like you would for a data.frame
#'
#' @param .Object \code{gTrack} object to extracting the formatting data.frame from
#' @docType methods
#' @rdname formatting-methods
#' @export
#' @author Marcin Imielinski
setGeneric('formatting', function(.Object) standardGeneric('formatting'))

#' @name xaxis
#' @title Retrieves the xaxis parameters
#' @description
#'
#' Return the portion of the gTrack @format field responsible for
#' formatting the x-axis
#' @param .Object \code{gTrack} object to retrieve xaxis parameters from
#' @return \code{data.frame} which is a subset of \code{formatting} showing xaxis params
#' @author Jeremiah Wala
#' @docType methods
#' @rdname xaxis-methods
#' @export
setGeneric('xaxis', function(.Object) standardGeneric('xaxis'))

#' @rdname xaxis-methods
#' @aliases xaxis,gTrack-method
setMethod('xaxis', 'gTrack', function(.Object) {
  xaxis.fields <- c("xaxis.prefix", "xaxis.suffix", "xaxis.unit",
                    "xaxis.round", "xaxis.interval", "xaxis.pos",
                    "xaxis.nticks", "xaxis.pos.label",
                    "xaxis.cex.label","xaxis.newline", "xaxis.chronly",
                    "xaxis.width", "xaxis.label.angle","xaxis.ticklen")
  return(.Object@formatting[,which(colnames(.Object@formatting) %in% xaxis.fields)])
})

#' @name sep
#' @title Retrieves the seperator graphical parameters
#' @description
#'
#' Return the portion of the gTrack @format field responsible for
#' formatting the windows and their separators
#' @author Jeremiah Wala
#' @param .Object \code{gTrack} object to extract separator params from
#' @return \code{data.frame} with the seperator paramters
#' @docType methods
#' @rdname sep-methods
#' @export
setGeneric('sep', function(.Object) standardGeneric('sep'))

#' @rdname sep-methods
#' @aliases sep,gTrack-method
setMethod('sep', 'gTrack', function(.Object) {
  sep.fields <- c("sep.lty", "sep.lwd", "sep.bg.col", "sep.draw")
  return(.Object@formatting[,which(colnames(.Object@formatting) %in% sep.fields)])

})

#' @rdname formatting-methods
#' @aliases formatting,gTrack-method
setMethod('formatting', 'gTrack', function(.Object)
{
  return(.Object@formatting)
})

#' @name edgs
#' @title edgs
#' @description
#'
#' Get edges list associated with gTrack object, is is a list of data.frames
#' each which contains a data.frame specifying formatted edges, as pairs of indices into the
#' underlying GRanges object $from, $to and additional formats $col $lwd $lty
#'
#' usage:
#' edgs(gt)
#' @param .Object \code{gTrack} object to extract the edge list from
#' @docType methods
#' @rdname edgs-methods
#' @export
setGeneric('edgs', function(.Object) standardGeneric('edgs'))


#' @rdname edgs-methods
#' @aliases edgs,gTrack-method
setMethod('edgs', 'gTrack', function(.Object)
{
  #             .Object = list(...)[[1]]
  return(.Object@edges)
})


setGeneric('vars', function(.Object) standardGeneric('vars'))

#' @name vars
#' @title vars
#' @description
#'
#' Get variants associated with track
#'
#' @keywords internal
#' @author Marcin Imielinski
setMethod('vars', 'gTrack', function(.Object)
{
  #             .Object = list(...)[[1]]
  return(.Object@vars)
})

#' @name edgs<-
#' @title edgs<-
#' @description
#' Set edges data.frame associated with associated with gTrack object, the data.frame
#' that is being used to replace must have fields $from, $to, and can have optional fields $lwd
#' $lty, $col specifyign color and line type.
#'
#' usage:
#' edgs(gt)[[1]] <- new.edges.
#'
#' @param .Object \code{gTrack} object on which to set new edges
#' @param value New value of the edges
#' @docType methods
#' @rdname edgs-set-methods
#' @export
#' @author Marcin Imielinski
setGeneric('edgs<-', function(.Object, value) standardGeneric('edgs<-'))

#' @rdname edgs-set-methods
#' @aliases edgs<-,gTrack-method
setMethod('edgs<-', 'gTrack', function(.Object, value)
{
  if (!all(sapply(value, is.data.frame)))
    stop('edges attribute must be a list of data.frames')

  non.empty = sapply(value, nrow)>0

  if (any(non.empty) )
    if (!all(sapply(value[non.empty], function(x) all(c('from', 'to') %in% colnames(x)))))
      stop('data.frames of edges field of gTrack must be either empty or have $from and $to fields specified')

  .Object@edges = value
  validObject(.Object)
  return(.Object)
})

#' @name clear
#' @title clear
#' @description
#'
#' Clear data from gTrack object (makes into an empty container with the same seqinfo and display features)
#' usage:
#' clear(gTrack)
#'
#' @param .Object \code{gTrack} object to clear
#' @author Marcin Imielinski
#' @docType methods
#' @rdname clear-methods
#' @export
setGeneric('clear', function(.Object) standardGeneric('clear'))

#' @rdname clear-methods
#' @aliases clear,gTrack-method
setMethod('clear', 'gTrack', function(.Object)
{
  .Object = .Object[1]
  .Object@data[[1]] =  .Object@data[[1]][NULL]
  validObject(.Object)
  return(.Object);
})


#' @name dat
#' @title dat
#' @description
#'
#' Extract list of data contained inside a gTrack object gt.  Each list item will contain
#' either a GRanges, RleList, or path to a file.
#'
#' usage
#' dat(gt)
#'
#' @param .Object \code{gTrack} to retrive data from
#' @docType methods
#' @rdname dat-methods
#' @export
#' @author Marcin Imielinski
setGeneric('dat', function(.Object) standardGeneric('dat'))

#' @rdname dat-methods
#' @aliases dat,gTrack-method
setMethod('dat', 'gTrack', function(.Object)
{
  return(.Object@data)
})

#' @name formatting<-
#' @title formatting<-
#' @description
#'
#' Set formatting of gTrack object  gt
#'
#' usage:
#' #gt is a gTrack object
#' formatting(td)$height <- 2
#'
#' #can also just use $ directly, like for a data frame, which is more convenient
#' td$height <- 2
#' @param .Object \code{gTrack} to set the formatting field for
#' @param value \code{data.frame} with the new formatting information
#' @docType methods
#' @rdname formatting-set-methods
#' @export
setGeneric('formatting<-', function(.Object, value) standardGeneric('formatting<-'))

#' @rdname formatting-set-methods
#' @aliases formatting<-,gTrack-method
setReplaceMethod('formatting', 'gTrack', function(.Object, value)
{
  REQUIRED.COLUMNS = c('height', 'col', 'ygap', 'y.field');
  if (nrow(value) != length(.Object))
    stop('Replacement data frame has %s rows and the object has %s', nrow(value), length(value))

  if (length(setdiff(REQUIRED.COLUMNS, colnames(value))))
    stop('Replacement data frame is missing some columns')
  .Object@formatting = value;

  validObject(.Object)
  return(.Object)
});


#' @name colormap
#' @title colormap
#' @description
#' Access colormap of gTrack object, this is a named list of named character vectors that specifies
#' the field of the underlying GRanges object that will be used to map a set of values
#' to a set of colors.
#' usage:
#' colormap(gt)
#'
#' @param .Object \code{gTrack} object to retrieve colormap from
#' @docType methods
#' @rdname colormap-methods
#' @export
#' @author Marcin Imielinski
setGeneric('colormap', function(.Object) standardGeneric('colormap'))

#' @rdname colormap-methods
#' @aliases colormap,gTrack-method
setMethod('colormap', 'gTrack', function(.Object)
{
  return(.Object@colormap)
})


#' @name colormap<-
#' @title colormap<-
#' @description
#' Set colormap of gTrack object, this is a named list of named character vectors that specifies
#' the field of the underlying GRanges object that will be used to map a set of values
#' to a set of colors.
#'
#' usage:
#' colormap(gt)[1] = list(tumortype = c(lung = 'red', pancreatic = 'blue', colon = 'purple'))
#' @param .Object \code{gTrack} to set the colormap for
#' @param value New \code{colormap}
#' @docType methods
#' @rdname colormap-set-methods
#' @export
setGeneric('colormap<-', function(.Object, value) standardGeneric('colormap<-'))

#' @rdname colormap-set-methods
#' @aliases colormap<-,gTrack-method
setReplaceMethod('colormap', 'gTrack', function(.Object, value)
{
  .Object@colormap = value;
  validObject(.Object)
  return(.Object)
});
#' @name show
#' @title show
#' @description Display a \code{gTrack} object
#' @docType methods
#' @param object \code{gTrack} to display
#' @rdname show-methods
#' @aliases show,gTrack-method
#' @export
#' @author Marcin Imielinski
setMethod('show', 'gTrack', function(object)
{
  cat(sprintf('gTrack object with %s tracks with formatting:\n', length(object)))
  print(formatting(object))
})

### utility function to allow "softer" matching of seqinfo's of genomes in chains
.identical.seqinfo = function(a, b)
{
  df.a = data.frame(seqnames = GenomeInfoDb::seqnames(a), seqlengths = GenomeInfoDb::seqlengths(a), stringsAsFactors = F);
  df.b = data.frame(seqnames = GenomeInfoDb::seqnames(b), seqlengths = GenomeInfoDb::seqlengths(b), stringsAsFactors = F);

  df.a = df.a[order(df.a$seqnames), ]
  df.b = df.b[order(df.b$seqnames), ]

  rownames(df.a) = NULL
  rownames(df.b) = NULL

  return(identical(df.a, df.b))
}

# @importFrom graphics plot
if (!isGeneric("plot"))
   setGeneric("plot", function(x, ...) standardGeneric("plot"))

#' @name plot
#' @title plot
#' @description
#' Plot gTrack object in multi-track genome browser view.  gTracks will appear as stacked interval
#' plots, line / bar / scatter plots, node-edge graphs, or triangle plots depending on the formatting
#' settings and additional data (eg mdata) provided.  gTracks can be drawn across several non-contiguous
#' "windows" of the genome
#'
#' Additional argument "links" takes a GRangesList of signed interval pairs (ie each item is length 2)
#' representing genomic junctions as input. The junctions will be pointing right for + intervals and
#' left for - intervals.  The meta data of the links GRangesList ($lwd $col $lty) can be set
#' to modulate the color, witwidth, and line style of the junction links .
#'
#' usage:
#' display(gt, win) # where win is a GRanges object
#' display(gt) # this will show the entire span of seqinfo(gt)
#' display(gt, '3:1e6-2e6') ## can use UCSC style strings
#' display(gt, GRanges(20, IRanges(20e6, 21e6)))
#' display(gt, '2', links = ra) # here, the entire chromosome 2 is shown, with rearrangements on top
#'
#'
#' @param windows GRanges specifying windows to view (can also be GRangesList), default is whole genome
#' @param links GRangesList of signed locus pairs specifying links to draw above the plot,
#' optional GRangesList values meta data specify formatting of individual links:
#' $label (text label), $col (line color), $lwd (line weight), $lty (line style),
#' $arrow (arrow flag), $col.arrow (arrow color), $v (vertical spline bulge), $w (horizontal spline bulge)
#' @param gap scalar numeric specifying gap between windows (only relevant if windows has length>1). Units of the gap are in genome coordinates
#' @param max.ranges scalar numeric > 0 specifying max number of ranges to draw in a window (via sampling).  If specified, overrides gTrack max.ranges formatting feature.
#' @param ..., additional last-minute formatting changes to the gtrack can be entered here (eg col = 'blue')
#'
#' @docType methods
#' @rdname plot-methods
#' @author Marcin Imielinski, Jeremiah Wala
#' @aliases plot,gTrack,ANY-method
#' @export
setMethod('plot', c("gTrack","ANY"),
          #signature(x = "gTrack", y = "ANY"),
          function(x,  ##pplot  (for easy search)
                   y,
                   windows = si2gr(seqinfo(x)), ## windows to plot can be Granges or GRangesList
                   links = NULL, ## GRangesList of pairs of signed locations,
                   gap = NULL,  ## spacing betwen windows (in bp)
                   y.heights = NULL, # should be scalar or length(windows) if windows is a GRangesList
                   y.gaps = NULL, # relative heights of gaps below and above stacks (xaxes will be drawn here)
                   cex.xlabel = 1,
                   cex.ylabel = 1,
                   max.ranges = NA, # parameter for max ranges to draw on canvas in each track (overrides formatting)
                   links.feat = NULL, # links features override for links (must be nrow 1 or length(links) data frame
                   verbose=FALSE,
                   legend.params = list(),
                   ... ## additional args to draw.grl OR last minute formatting changes to gTrack object
                   )
{
  if (!missing(y))
    windows = y

  .Object = x
  if (!missing(y))
    windows = y

  
  ## make sure we have min legend data
  if (!"xpos" %in% names(legend.params))
    legend.params$xpos = 0
  if (!"ypos" %in% names(legend.params))
    legend.params$ypos = 1
  if (!"plot" %in% names(legend.params))
    legend.params$plot = TRUE

  win.gap = gap ## recasting some variable names
  new.plot = TRUE
  window.segs = list();
  dotdot.args = list(...);

  ## parse the windows into GRanges
  windows = format_windows(windows, .Object)

    ## if totally empty, plot blank and leave
  if(!length(windows)) {
    plot.blank(bg.col = bg.col)
    return()
  }

  ## make sure gTrack has all fields that are expected later
  .Object <- prep_defaults_for_plotting(.Object)
  
  if (is.null(formatting(.Object)$legend))
      formatting(.Object)$legend = TRUE
  else (any(is.na(formatting(.Object)$legend)))
       formatting(.Object)$legend[is.na(formatting(.Object)$legend)] = TRUE


  if (is.null(formatting(.Object)$legend.title))
      formatting(.Object)$legend.title = NA
  
  has.colormap = sapply(colormap(.Object), length)>0
  has.colorfield = !is.na(formatting(.Object)$gr.colorfield)
  formatting(.Object)$legend = formatting(.Object)$legend == TRUE & (has.colormap | has.colorfield)
  numlegends = sum(formatting(.Object)$legend)
  which.legend = which(formatting(.Object)$legend)
  .Object$legend.title = ifelse(!is.na(.Object$legend.title), .Object$legend.title,
      ifelse(has.colormap, names(colormap(.Object))[1:length(.Object)],
             ifelse(has.colorfield, .Object$gr.colorfield,
                    ifelse(!is.na(.Object$name), .Object$name, ''))))
  
    ## add last minute formatting changes to gTrack
  if (length(dotdot.args)>0)
    for (f in intersect(names(dotdot.args), names(formatting(.Object))))
      formatting(.Object)[, f] = dotdot.args[[f]]

  ## set the window gap
  if (is.null(win.gap))
    win.gap = sapply(windows, function(x) {wx = width(x); min(sum(as.numeric(wx))/length(wx)/10, sum(as.numeric(wx))/20)})

  ## get the height of the stacks
  if (is.null(y.heights) | length(y.heights) != length(windows))
    ##y.heights = rep(1, length(windows)) ## old from when we had windows as GRangesList
    y.heights <- 1

  ## set the gaps between the gTracks
  if (is.null(y.gaps) )
    y.gaps = y.heights*0.8
  else if (length(y.gaps) != length(windows))
    y.gaps = rep(y.gaps[1], length(windows))

  ## ensure that we don't plot too much
  if (!is.na(max.ranges))
    formatting(.Object)$max.ranges = pmin(max.ranges, formatting(.Object)$max.ranges, na.rm = TRUE)

  oth.ix = 1:(length(windows)-1);
  top.gaps = 0.5*y.gaps
  bottom.gaps = 0.5*y.gaps
  if (length(windows)==1)
    oth.ix = c()
  ylim.stacks = data.frame(start = c(bottom.gaps[1], bottom.gaps[1] + cumsum(y.heights[oth.ix] + top.gaps[oth.ix] + bottom.gaps[oth.ix+1])),
                           end = cumsum(y.heights + top.gaps + bottom.gaps) - top.gaps)

  oth.ix = 1:(length(.Object)-1);

  if (length(.Object)==1)
    oth.ix = c()

  tmp.top.gaps = 0.5 * formatting(.Object)$ygap
  tmp.bottom.gaps = 0.5 * formatting(.Object)$ygap
  tmp.ylim.subplot = data.frame(start = c(tmp.bottom.gaps[1], tmp.bottom.gaps[1] +

                                            cumsum(formatting(.Object)$height[oth.ix] + tmp.top.gaps[oth.ix] + tmp.bottom.gaps[oth.ix+1])),
                                end = cumsum(formatting(.Object)$height + tmp.top.gaps + tmp.bottom.gaps) - tmp.top.gaps)

  ylim = c(0, max(ylim.stacks$end)+top.gaps[length(top.gaps)])
  ylim.parent=ylim
  window.ylims = data.frame(start = rep(NA, length(windows)), end = NA);

  new.axis = TRUE;
  this.windows = windows
  #end(this.windows) <- end(this.windows) + 1 ## +1 added
  #this.windows = gUtils::streduce(windows[[i]]) ##gr.stripstrand(GenomicRanges::trim(windows[[i]]))
  ##if (!inherits(this.windows, 'GRanges'))
  ##  this.windows = gUtils::si2gr(this.windows)
  i=1
  this.ylim.subplot = tmp.ylim.subplot;
  this.ylim.subplot$start = affine.map(pmin(1, this.ylim.subplot$start), ylim = unlist(ylim.stacks[i, c('start', 'end')]), xlim = c(0, 1))
  this.ylim.subplot$end = affine.map(pmin(1, this.ylim.subplot$end), ylim = unlist(ylim.stacks[i, c('start', 'end')]), xlim = c(0, 1))
  this.tmp.bottom.gap = tmp.bottom.gaps[1]*(ylim.stacks$end[i]-ylim.stacks$start[i])

  this.xaxis.pos = this.ylim.subplot$start[1]-bottom.gaps[i]*0-this.tmp.bottom.gap
  this.xaxis.pos.label = this.ylim.subplot$start[1]-5*bottom.gaps[i]/6-this.tmp.bottom.gap
  ylim.stacks[i, 'xaxis.pos'] = this.xaxis.pos

  ## loop through the gTracks
  for (j in 1:length(.Object))
  {
    par(xpd = NA);
    cmap = colormap(.Object)[[j]];
    cfield = names(colormap(.Object))[j]

    if (is.na(cfield))
      cfield = formatting(.Object)$gr.colorfield[j]

    if (length(cmap)==0)
      cmap = NA

    ## get the data into GRanges or GRangesList format
    tt <- extract_data_from_tmp_dat(.Object, j, this.windows)
    .Object = tt$o
    tmp.dat = tt$t

    ## flag to tell us whether data is pre-filtered to window (ie in fftrack or rlelist)
    pre.filtered = FALSE;
    if (.Object@formatting$triangle[j])
      pre.filtered = TRUE

    ## subsample if we need to for enforcing max.ranges
    if (!is.na(formatting(.Object)$max.ranges[j]) && formatting(.Object)$max.ranges[j] > 0) {
      tt <- enforce_max_ranges(.Object, pre.filtered, j, tmp.dat, this.windows)
      tmp.dat = tt$t
      pre.filtered = tt$p
    }

    ## adjust y0 .bar
    if (is.null((formatting(.Object)$y0.bar[j])) || is.na((formatting(.Object)$y0.bar[j])))
      formatting(.Object)$y0.bar[j] = 0

    ## smooth the y.field data
    if (!is.na(formatting(.Object)$y.field[j]) && is(tmp.dat, 'GRanges') && !is.na(formatting(.Object)$smooth[j]))
      tmp.dat <- smooth_yfield(.Object, j, tmp.dat)

    ## fix y limits and apply log transform if needed
    if (!is.na(formatting(.Object)$y.field[j]) && (is.na(formatting(.Object)$y0[j]) || is.na(formatting(.Object)$y1[j])))
      .Object <- format_yfield_limits(.Object, j, tmp.dat, pre.filtered, this.windows)

    # if (formatting(.Object[j])$format != 'ranges')
    #   stop("violated assumption. need to fix")

    all.args = as.list(formatting(.Object[j]))
    all.args = c(dotdot.args, all.args[!names(all.args) %in% names(dotdot.args)])

    #                     all.args = list(
    #   col = formatting(.Object)$col[j],
    #   ywid = formatting(.Object)$ywid[j],
    #   border = formatting(.Object)$border[j],
    #   lwd.border = formatting(.Object)$lwd.border[j],
    #   adj.label = c(formatting(.Object)$hadj.label[j],
    #     formatting(.Object)$vadj.label[j]),
    #   gr.adj.label = c(0.5,
    #     formatting(.Object)$vadj.label[j]),
    #   angle = formatting(.Object)$angle[j],
    #   y.pad = formatting(.Object)$ypad[j],
    #   circles = formatting(.Object)$circles[j],
    #   lines = formatting(.Object)$lines[j],
    #   bars = formatting(.Object)$bars[j],
    #   y.grid.cex = formatting(.Object)$yaxis.cex[j],
    #   edges = edgs(.Object)[[j]],
    #   y0.bar = formatting(.Object)$y0.bar[j],
    #   stack.gap = formatting(.Object)$stack.gap[j])
    # na.fields = names(formatting(.Object))[sapply(1:ncol(formatting(.Object)), function(field) is.na(formatting(.Object)[j, field]))]
    # other.fields = setdiff(names(formatting(.Object)), c('name', 'height', 'ygap', 'stack.gap', 'lift', 'split', 'angle', 'format', 'lwd.border', 'source.file', 'source.file.chrsub', 'ypad', 'ywid', 'border', 'col', 'hadj.label', 'vadj.label', 'y.field', 'round', 'cex.ylabel', 'y.quantile', 'max.ranges', 'yaxis', 'yaxis.cex', 'is.null', 'yaxis.pretty', names(all.args))) ## remove na fields and anything else that might mess up draw.grl
    #
    # other.formats = structure(names = other.fields,
    #   lapply(other.fields, function(x) formatting(.Object)[j, x]))
    # all.args[names(other.formats)] = other.formats;
    #
    # all.args = all.args[setdiff(names(all.args), setdiff(na.fields, c('col')))]

    this.y.field = formatting(.Object)$y.field[j]
    this.y.grid = NA;

    if (is.na(this.y.field) | !(this.y.field %in% names(values(tmp.dat))))
      this.y = this.ylim.subplot[j, ]
    else
    {
      if (is.null(formatting(.Object)$log))
        formatting(.Object)$log = NA

      if (!is.na(formatting(.Object)$log[j]))
        if (formatting(.Object)$log[j])
        {
          if (!is.null(tmp.dat$ywid))
            tmp.dat$ywid = log10(tmp.dat$ywid)
          values(tmp.dat)[, this.y.field] = log10(values(tmp.dat)[, this.y.field])
          formatting(.Object)[j, 'y0'] = log10(formatting(.Object)[j, 'y0'])
          formatting(.Object)[j, 'y1'] = log10(formatting(.Object)[j, 'y1'])

        }
      range.y = NULL;
      if (all(c('y0', 'y1') %in% names(formatting(.Object))))
      {
        if (!is.na(formatting(.Object)[j, 'y0']) & !is.na(formatting(.Object)[j, 'y1']))
          range.y = c(formatting(.Object)[j, 'y0'], formatting(.Object)[j, 'y1'])
        else if (!is.na(formatting(.Object)[j, 'y0']) & is.na(formatting(.Object)[j, 'y1']))
          range.y = c(formatting(.Object)[j, 'y0'], max(setdiff(values(tmp.dat)[, this.y.field], c(Inf, -Inf)), na.rm = T))
        else if (is.na(formatting(.Object)[j, 'y0']) & !is.na(formatting(.Object)[j, 'y1']))
          range.y = c(min(setdiff(values(tmp.dat)[, this.y.field], c(Inf, -Inf)), na.rm = T), formatting(.Object)[j, 'y1'])
      }

      if (!is.null(tmp.dat$ywid)) ## remove any weird infinite ywids
        if (any(ix <- is.infinite(tmp.dat$ywid)))
          tmp.dat$ywid[ix] = NA

      if (is.null(range.y)) ## if y range is empty then pull from data
        range.y = range(setdiff(values(tmp.dat)[, this.y.field], c(Inf, -Inf)), na.rm = T);
      #                              range.y = range(setdiff(values(dat(.Object)[[j]])[, this.y.field], c(Inf, -Inf)), na.rm = T);

      ## however if there is a single data value, then we need to scale appropriately
      if (diff(range.y)==0)
      {
        if (any(ix <- !is.na(tmp.dat$ywid))) ## use ywid
          range.y = range.y + 5*max(tmp.dat$ywid[ix])*c(-1, 1)
        else # otherwise use some arbitrary proportion around the value
          range.y = range.y + abs(range.y)*0.2*c(-1, 1)
      }
      
      this.y.ticks = pretty(range.y, formatting(.Object)$yaxis.pretty[j])

      if (is.null(formatting(.Object)$y.cap))
          formatting(.Object)$y.cap = NA
      
      if (!is.na(formatting(.Object)$y.cap[j])) ## cap values from top and bottom
        this.y = affine.map(values(tmp.dat)[, this.y.field], ylim = unlist(this.ylim.subplot[j, ]), xlim = range(this.y.ticks), cap = formatting(.Object)$y.cap[j])
      else
        this.y = affine.map(values(tmp.dat)[, this.y.field], ylim = unlist(this.ylim.subplot[j, ]), xlim = range(this.y.ticks), cap = TRUE)
      #                            this.y = affine.map(values(dat(.Object)[[j]])[, this.y.field], ylim = unlist(this.ylim.subplot[j, ]), xlim = range(this.y.ticks))

      ## if need, bump the range to include ybar base
      # if (formatting(.Object)$y0.bar[j] < min(unlist(this.ylim.subplot[j, ])))
      #   this.ylim.subplot[j,'start'] <- formatting(.Object)$y0.bar[j]

      all.args$y0.bar = affine.map(formatting(.Object)$y0.bar[j], ylim = unlist(this.ylim.subplot[j, ]), xlim = range(this.y.ticks))
      if (formatting(.Object)$yaxis[j])
      {
        # make pretty grid in range.y
        this.y.grid = structure(affine.map(this.y.ticks, ylim = unlist(this.ylim.subplot[j, ]), xlim = range(this.y.ticks)), names = this.y.ticks)

        if (!is.na(formatting(.Object)$log[j]))
          if (formatting(.Object)$log[j])
            names(this.y.grid) = signif(10^this.y.ticks)
      }

      ## ## fix ywids if necessary
      ## if (!is.null(values(tmp.dat)$ywid))
      ##   tmp.dat$ywid = tmp.dat$ywid * (this.ylim.subplot[j, 2]-this.ylim.subplot[j, 1])/diff(range.y)
    }

    if (.Object[j]$bars && is.na(all.args$y0.bar))
        all.args$y0.bar = this.ylim.subplot[j, 1]

    if (.Object[j]$chr.sub)
        tmp.windows = gr.sub(windows, 'chr', '')
    else
        tmp.windows = this.windows

    ## fix legend params
    this.legend.params = legend.params
    if (!formatting(.Object)$legend[j])
        this.legend.params$plot = FALSE
    else
        {            
            this.legend.params$xpos = NA
            if (!is.null(formatting(.Object)$legend.xpos))
                this.legend.params$xpos = is.formatting(.Object)$legend.xpos[j]
            if (!is.null(formatting(.Object)$legend.ypos))
                this.legend.params$ypos = formatting(.Object)$legend.ypos[j]
            
            jj = match(j, which.legend)
            if (is.na(this.legend.params$xpos))
                this.legend.params$xpos = seq(0, 1, length.out = numlegends)[jj]

            this.legend.params$xpos = seq(0, 1, length.out = numlegends)[jj]
            this.legend.params$xjust = c(0, 0.5, 1)[as.integer(cut(this.legend.params$xpos, c(-0.01, 0.2, 0.8, 1)))]
            this.legend.params$title = .Object$legend.title[j]
        }         

    main.args <- list(grl=tmp.dat,y = this.y, ylim = ylim,
                      xaxis.pos = this.xaxis.pos,xaxis.pos.label = this.xaxis.pos.label,
                      win.gap = win.gap[i],windows = tmp.windows,
                      new.plot = new.plot, new.axis = new.axis,
                      gr.colorfield = cfield,gr.colormap = cmap,
                      y.grid = this.y.grid, verbose=verbose,
                      ylim.parent=ylim.parent,mdata=.Object@mdata[[j]],
                      leg.params = this.legend.params,
                      adj.label = c(formatting(.Object)$hadj.label[j],
                                    formatting(.Object)$vadj.label[j]),
                      gr.adj.label = c(0.5,
                                       formatting(.Object)$vadj.label[j]),
                      y.pad = formatting(.Object)$ypad[j],
                      y.grid.cex = formatting(.Object)$yaxis.cex[j],
                      edges = edgs(.Object)[[j]])
    all.args <- c(main.args, all.args[!names(all.args) %in% names(main.args)])
    
    # main.args = c(main.args, all.args[setdiff(names(all.args), names(main.args))]);
    #
    # other.args = dotdot.args
    # other.args = other.args[setdiff(names(other.args), names(main.args))]
    ## make empty plot

    if (new.plot)
    {
      blank.main.args <- all.args
      ###blank.main.args = main.args;
      ###blank.main.args[[1]] = blank.main.args[[1]][c()]
      #blank.main.args$y = list(start = min(this.ylim.subplot$start), end = max(this.ylim.subplot$end))

      # if (any(is.na(.Object@formatting$triangle)) & any(.Object@formatting$triangle))
      #   blank.main.args$y = list(start = min(this.ylim.subplot$start[is.na(.Object@formatting$triangle)]), end = max(this.ylim.subplot$end[is.na(.Object@formatting$triangle)])) ## JEREMIAH
      # else
      blank.main.args$grl <- GRanges()
      blank.main.args$y = list(start=min(this.ylim.subplot$start), end = max(this.ylim.subplot$end))
      blank.main.args$triangle=FALSE

      do.call('draw.grl', blank.main.args)
      ##do.call('draw.grl', c(blank.main.args, other.args))
      all.args$new.plot = FALSE
      all.args$new.axis = FALSE
    }

    new.plot = FALSE
    new.axis = FALSE
    ##arrg <- c(main.args, other.args)
    ##form <- as.list(formatting(.Object[j]))
    ##arrg <- c(arrg, form[!names(form) %in% names(arrg)]) ## add in remaining args

    window.segs = list();
    ##window.segs[[i]] = do.call('draw.grl', c(main.args, other.args))
    if (formatting(.Object[j])$triangle)
    {
      window.segs[[i]] <- do.call('draw.triangle', all.args[names(all.args) %in% c("grl","y","mdata","ylim.parent","windows","win.gap","sigma",
                                                                                   "cmap.min","cmap.max", "m.sep.lwd","m.bg.col","leg.params",
                                                                                   "islog","gr.colormap")])
    } else {
      window.segs[[i]] <- do.call('draw.grl', all.args)
    }

    this.tname = formatting(.Object[j])$name

    if (!is.na(this.tname))
    {
      this.cex.ylabel = ifelse(!is.null(formatting(.Object[j])$cex.ylabel), formatting(.Object[j])$cex.ylabel, cex.ylabel)
      text(par('usr')[2], mean(unlist(this.ylim.subplot[j, c('start', 'end')])),
           this.tname, srt = -90, adj = c(0.5, 1), cex = this.cex.ylabel)
    }

  }

  if (is.null(links))
    links = GenomicRanges::GRangesList()

  if (length(links)>0) # draw rearrangement links
  {
    # first map rearrangements to various windows>
    win.u = this.windows
    win.u$grl.ix = 1  ##holdover from grangeslist windows
    ##win.u = gr.stripstrand(grl.unlist(windows))
    window.segs.u = do.call(rbind, window.segs)
    window.segs.u$width = window.segs.u$end - window.segs.u$start + 1
    window.segs.xlim = do.call('rbind', lapply(window.segs, function(x) data.frame(start = min(x$start), end = max(x$end))))

    links.u = grl.unlist(links)

    if (any(table(links.u$grl.ix)!=2))
      stop('Links should be GRangesList of range pairs.')

    links.p = grl.pivot(links)

    ## flip strands to conform to connectors convention (- connection to left side and + is connection to right side)
    ## with ra specification convention (- refers to segment to left of breakpoint and + refers to segment to right)
    GenomicRanges::strand(links.p[[1]]) = c("-" = "+", "+" = "-")[as.character(GenomicRanges::strand(links.p[[1]]))]
    GenomicRanges::strand(links.p[[2]]) = c("-" = "+", "+" = "-")[as.character(GenomicRanges::strand(links.p[[2]]))]

    ## find overlaps with windows and calculate their window specific coordinates
    l1 = gr.findoverlaps(links.p[[1]], win.u)
    values(l1) = cbind(as.data.frame(values(l1)), as.data.frame(values(links)[l1$query.id, , drop = FALSE]))
    GenomicRanges::strand(l1) = GenomicRanges::strand(links.p[[1]])[l1$query.id]
    l1$stack.id = win.u$grl.ix[l1$subject.id]
    l1$y.pos = ylim.stacks$xaxis.pos[l1$stack.id]
    l1$x.pos = mapply(function(x,y,z,a) (y-z)*a + x, x = window.segs.u[l1$subject.id,]$start, y = start(l1),
                      z = start(win.u[l1$subject.id]), a = window.segs.u$width[l1$subject.id] / width(win.u)[l1$subject.id])

    l2 = gr.findoverlaps(links.p[[2]], win.u)
    values(l2) = cbind(as.data.frame(values(l2)), as.data.frame(values(links)[l2$query.id, , drop = FALSE]))
    GenomicRanges::strand(l2) = GenomicRanges::strand(links.p[[2]])[l2$query.id]
    l2$stack.id = win.u$grl.ix[l2$subject.id]
    l2$y.pos = ylim.stacks$xaxis.pos[l2$stack.id]
    l2$x.pos = mapply(function(x,y,z,a) (y-z)*a + x, x = window.segs.u[l2$subject.id,]$start, y = start(l2),
                      z = start(win.u[l2$subject.id]), a = window.segs.u$width[l2$subject.id] / width(win.u)[l2$subject.id])

    
    .fix.l = function(ll)
    {
      if (!is.null(links.feat))
        for (col in names(links.feat))
          values(ll)[, col] = links.feat[ll$query.id, col]

        # set up connectors
        if (is.null(ll$v))
          ll$v = y.gaps[ll$stack.id]/4
        else
          ll$v = y.gaps[ll$stack.id]*ll$v/2

        if (is.null(ll$h))
          ll$h = (window.segs.xlim$end[ll$stack.id] - window.segs.xlim$start[ll$stack.id])/20
        else
          ll$h = (window.segs.xlim$end[ll$stack.id] - window.segs.xlim$start[ll$stack.id])*ll$h

        if (is.null(ll$arrow))
          ll$arrow = TRUE

        if (is.null(ll$cex.arrow))
          ll$cex.arrow = 1

        if (is.null(ll$lwd))
          ll$lwd = 1


        if (is.null(ll$lty))
          ll$lty = 3


        if (is.null(ll$col))
          ll$col = 'red'

        if (is.null(ll$col.arrow))
          ll$col.arrow = ll$col

        return(ll)
    }

    if (length(l1)>0)
      l1 = .fix.l(l1)

    if (length(l2)>0)
      l2 = .fix.l(l2)


    ## now pair up / merge l1 and l2 using query.id as primary key
    if (length(l1)>0 & length(l2)>0)
    {
      pairs = merge(data.frame(l1.id = 1:length(l1), query.id = l1$query.id), data.frame(l2.id = 1:length(l2), query.id = l2$query.id))[, c('l1.id', 'l2.id')]
      l1.paired = GenomicRanges::as.data.frame(l1)[pairs[,1], ]
      l2.paired = GenomicRanges::as.data.frame(l2)[pairs[,2], ]
    }
    else
    {
      l2.paired = l1.paired = data.frame()
      pairs = data.frame()
    }


    ## some l1 and l2 will be unpaired
    l.unpaired = GRanges(seqlengths = GenomeInfoDb::seqlengths(links));
    p1 = p2 = c();
    if (nrow(pairs)>0)
    {
      p1 = pairs[,1]
      p2 = pairs[,2]
    }

    if (length(l1)>0)
      l.unpaired = c(l.unpaired, l1[setdiff(1:length(l1), p1)])

    if (length(l2)>0)
      l.unpaired = c(l.unpaired, l2[setdiff(1:length(l2), p2)])

    if (length(l.unpaired)>0)
    {
      l.unpaired$v = l.unpaired$v/4
      l.unpaired$h = l.unpaired$h/2
      l.unpaired$col = alpha(l.unpaired$col, 0.5)
      ## unpaired will get a "bridge to nowhere" in the proper direction (eg down)
      l.unpaired$y.pos = ylim.stacks$end[l.unpaired$stack.id]
      #                    l.unpaired$y.pos2 = l.unpaired$y.pos + top.gaps[l.unpaired$stack.id]
      l.unpaired$y.pos2 = l.unpaired$y.pos + l.unpaired$v

      connectors(l.unpaired$x.pos, l.unpaired$y.pos, as.character(GenomicRanges::strand(l.unpaired)),
                 l.unpaired$x.pos, l.unpaired$y.pos2,
                 as.character(GenomicRanges::strand(l.unpaired)),
                 v = abs(l.unpaired$v), h = l.unpaired$h, type = 'S',
                 f.arrow = l.unpaired$arrow, b.arrow = l.unpaired$arrow,
                 cex.arrow = 0.2*l.unpaired$cex.arrow,
                 col.arrow = l.unpaired$col.arrow,
                 lwd = l.unpaired$lwd, lty = l.unpaired$lty, col = l.unpaired$col)

      if (!is.null(l.unpaired$label))
      {
        cex = 0.5;
        if (!is.null(l.unpaired$cex.label))
          cex = l.unpaired$cex.label

        #                        text(l.unpaired$x.pos, l.unpaired$y.pos2, l.unpaired$label, adj = c(0.5, 0.5), cex = cex)
        l.unpaired$text.y.pos =l.unpaired$y.pos - (ylim.stacks$end[l.unpaired$stack.id]-ylim.stacks$start[l.unpaired$stack.id])/100
        text(l.unpaired$x.pos, l.unpaired$text.y.pos, l.unpaired$label, adj = c(0.5, 1), cex = cex)
      }
    }

    ## now draw connectors for paired links
    ## pairs on the same level will get a "U" link,
    ## pairs on different levels will get "S" links

    ## fix y distances so that connections go from topmost part of bottom most connection to the bottom most part of
    ## top connection (ie the xaxis position)

    if (nrow(l1.paired)>0)
    {
      l1.paired$bottom = l1.paired$y.pos < l2.paired$y.pos

      if (any(l1.paired$bottom))
        l1.paired$y.pos[l1.paired$bottom] = ylim.stacks$end[l1.paired$stack.id[l1.paired$bottom]]

      if (any(!l1.paired$bottom))
        l2.paired$y.pos[!l1.paired$bottom] = ylim.stacks$end[l2.paired$stack.id[!l1.paired$bottom]]

      # also make all c type connections top connections with positive v
      ctype = ifelse(l1.paired$stack.id == l2.paired$stack.id, 'U', 'S')
      l1.paired$y.pos[ctype == 'U'] = ylim.stacks$end[l1.paired$stack.id[ctype == 'U']]
      l2.paired$y.pos[ctype == 'U'] = ylim.stacks$end[l2.paired$stack.id[ctype == 'U']]
      l1.paired$v[ctype == 'U'] = abs(l1.paired$v[ctype == 'U'])

      win.width = diff(par('usr')[1:2])
      l1.paired$v = l1.paired$v * affine.map(abs(l2.paired$x.pos - l1.paired$x.pos)/diff(par('usr')[1:2]), xlim = c(0, 1), ylim = c(0.5, 1)) ## make longer links taller
      connectors(l1.paired$x.pos, l1.paired$y.pos, l1.paired$strand, l2.paired$x.pos, l2.paired$y.pos, l2.paired$strand,
                 v = l1.paired$v, h = l1.paired$h, type = ctype,
                 f.arrow = l1.paired$arrow, b.arrow = l1.paired$arrow, cex.arrow = 0.2*l1.paired$cex.arrow, col.arrow = l1.paired$col.arrow,
                 lwd = l1.paired$lwd, lty = l1.paired$lty, col = l1.paired$col)

      if (!is.null(l1.paired$label))
      {
        cex = 0.5;
        if (!is.null(l1.paired$cex.label))
          cex = l1.paired$cex.label

        ##                         text(l1.paired$x.pos, l1.paired$y.pos+l1.paired$v/2, l1.paired$label, adj = c(0.5, 0.5), cex = cex)
        ##                         text(l2.paired$x.pos, l2.paired$y.pos+l1.paired$v/2, l2.paired$label, adj = c(0.5, 0.5), cex = cex)

        l1.paired$text.y.pos = l1.paired$y.pos - (ylim.stacks$end[l1.paired$stack.id]-ylim.stacks$start[l1.paired$stack.id])/100
        l2.paired$text.y.pos = l2.paired$y.pos - (ylim.stacks$end[l2.paired$stack.id]-ylim.stacks$start[l2.paired$stack.id])/100

        text(l1.paired$x.pos, l1.paired$text.y.pos, l1.paired$label, adj = c(0.5, 1), cex = cex)
        text(l2.paired$x.pos, l2.paired$text.y.pos, l2.paired$label, adj = c(0.5, 1), cex = cex)
      }
    }

  }
})

###############
# Tracked data helper functions
#
###############

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


# @name karyogram
# @title karyogram
# @description
#
# Returns gTrack object representing refGene transcripts and their components (utr, cds etc) with assigned colors.
# Usually built from cahced data objects but can also be built from provided GRangesList
#
# @param rg (optional) GRangesList representing transcript models obtained from refgene, with
# GrangesList meta data fields chr, s1, s2, e1, e2, str, gene_sym, Uniprot,
# @param genes (optional) character vector specifying genes to limit gTrack object to
# @param gene.collapse scalar logical specifying whether to collapse genes by transcript (or used stored version of transcripts)
# @param bg.col scalar character representing background color for genespan
# @param cds.col scalar character representing background color for CDS
# @param cds.utr scalar character representing background color for UTR
# @param st.col scalar character representing color of CDS start
# @param en.col scalar character representing color of CDS end
# @param genespan logical scalar whether to include genespan range around entire gne
# @param utr logical scalar whether to include range specifying UTR
# @param cds logical scalar whether to include range specifying CDS
# @param cached logical scalar whether to use "cached" version provided with package
# @param gr.srt.label scalar numeric specifying angle on exon label
# @param gr.cex.label scalar numeric > 0 specifying character expansion on exon label
# @param labels.suppress.gr scalar logical specifying whether to suppress exon label plotting
# @param stack.gap stack.gap argument to gTrack
# @param ... additional arguments passed down to gTrack
#
# @export
# @importFrom IRanges IRanges
# @importFrom GenomicRanges GRanges
# @author Marcin Imielinski
# track.refgene = function(rg = NULL,
#   gene.collapse = T,
#   genes = NULL,
#   bg.col = alpha('blue', 0.1), cds.col = alpha('blue', 0.6), utr.col = alpha('purple', 0.4),
#   st.col = 'green',
#   en.col = 'red',
#   genespan = T, ## flags determining whether to include ranges for these features
#   utr = T, ## show utr?
#   cds = T, ## show cds start / end?
#
#   grl.labelfield, ## Don't touch these
#   gr.labelfield,
#   col,
#   cached = T, ## if true will use cached version
#   cached.path = system.file("extdata", "refgene.composite.rds", package = 'gTrack'),  ## location of cached copy
#   cached.path.collapsed = system.file("extdata", "refgene.composite.collapsed.rds", package = 'gTrack'),  ## location of cached collapsed copy
#   gr.srt.label = 0,
#   gr.cex.label = 0.8,
#   labels.suppress.gr = T,
#   stack.gap = 1e6,
#   ...)
#   {
#
#     if (!cached | (!gene.collapse & !file.exists(cached.path)) | (gene.collapse  & !file.exists(cached.path.collapsed)))  ## if no composite refgene copy, then make from scratch
#       {
#         if (is.null(rg))
#           rg = read_refGene(grl = T)
#
#         values(rg)$label = values(rg)$gene_sym
#         ix.pos = values(rg)$str == '+'
# #         OUT.COLS = c('rg.id', 'type', 'exon_id', 'exon_frame', 'border', 'col', 'type')
#         OUT.COLS = c('rg.id', 'type', 'exon_id', 'exon_frame', 'type', 'is.exon')
#
#         rg.exons = gUtils::grl.unlist(rg)
#         rg.exons$rg.id = rg.exons$grl.ix
#         rg.exons$border = rg.exons$col = cds.col
#         rg.exons$is.exon = T;
#         rg.exons$type = 'exon'
#
#         if (genespan)
#           rg.genespan = GRanges(values(rg)$chr, IRanges(values(rg)$s1, values(rg)$e1), strand = values(rg)$str, seqlengths = hg_seqlengths(), rg.id = 1:length(rg), col = bg.col, type = 'genespan')
#         else
#           rg.genespan = NULL
#
#         if (utr)
#           {
#             type1 = ifelse(ix.pos, 'utr.5', 'utr.3')
#             type2 = ifelse(ix.pos, 'utr.3', 'utr.5')
#             tmp.utr = c(
#               GRanges(values(rg)$chr, IRanges(values(rg)$s1, pmax(values(rg)$s1, values(rg)$s2)), strand = values(rg)$str, seqlengths = hg_seqlengths(), rg.id = 1:length(rg), col = utr.col, border = utr.col, type = type1),
#               GRanges(values(rg)$chr, IRanges(values(rg)$e1, pmax(values(rg)$e1, values(rg)$e2)), strand = values(rg)$str, seqlengths = hg_seqlengths(), rg.id = 1:length(rg), col = utr.col, border = utr.col, type = type2))
#
#             ##utr.ix = which(rg.exons %over% tmp.utr)
#             utr.ix = which(gUtils::gr.in(rg.exons, tmp.utr))
#             utr.exons = rg.exons[utr.ix]
#             utr.exons$gr.ix = utr.ix
#             ix.check = merge(data.frame(i = 1:length(utr.exons), key = utr.exons$grl.ix), data.frame(j = 1:length(tmp.utr), key = tmp.utr$rg.id))
#             rg.utr.exons = pintersect(utr.exons[ix.check$i, ], tmp.utr[ix.check$j, ], resolve.empty = 'start')
#             non.empty = width(rg.utr.exons)!=0
#             rg.utr.exons = rg.utr.exons[non.empty]
#             values(rg.utr.exons) = values(utr.exons)[ix.check$i[non.empty], ]
#             rg.utr.exons$type = tmp.utr$type[ix.check$j[non.empty]]
#             rg.utr.exons$rg.id = rg.utr.exons$grl.ix
#             IRanges::ranges(rg.exons[rg.utr.exons$gr.ix]) = IRanges::ranges(psetdiff(rg.exons[rg.utr.exons$gr.ix], rg.utr.exons)) ## trim rg.exons by rg.utr.exons
#             ## remove width 0 rg.exons
#             rg.exons = rg.exons[width(rg.exons)>0]
#           }
#         else
#           rg.utr.exons = NULL
#
#         if (cds)
#           {
#             sten.col1 = ifelse(ix.pos, st.col, en.col)
#             sten.col2 = ifelse(ix.pos, en.col, st.col)
#             type1 = ifelse(ix.pos, 'cds.start', 'cds.end')
#             type2 = ifelse(ix.pos, 'cds.end', 'cds.start')
#             rg.sten = c(
#               GRanges(values(rg)$chr, IRanges(values(rg)$s2, pmax(values(rg)$s2, values(rg)$s2)), strand = values(rg)$str, seqlengths = hg_seqlengths(), rg.id = 1:length(rg), col = sten.col1, border = sten.col1, type = type1),
#               GRanges(values(rg)$chr, IRanges(values(rg)$e2, pmax(values(rg)$e2, values(rg)$e2)), strand = values(rg)$str, seqlengths = hg_seqlengths(), rg.id = 1:length(rg), col = sten.col2, border = sten.col2, type = type2))
#           }
#         else
#           rg.sten = NULL
#
#         tmp = rbind(rg.genespan, rg.exons, rg.utr.exons, rg.sten)[, OUT.COLS]
#         tmp$label = NULL;
#
#         ## compute tx ord of intervals
#         ord.ix = order(tmp$rg.id, match(tmp$type, c('genespan', 'utr.5', 'cds.start', 'exon', 'cds.end', 'utr.3')))
#         tmp.rle = rle(tmp$rg.id[ord.ix])
#         tmp$tx.ord[ord.ix] = unlist(lapply(tmp.rle$lengths, function(x) 1:x))
#
#         if (gene.collapse)
#           {
#             tmp = tmp[order(match(tmp$type, c('genespan', 'exon', 'utr.5', 'utr.3', 'cds.start', 'cds.end')))]
#             rg.composite = split(tmp, values(rg)[tmp$rg.id, ]$gene_sym)
#             values(rg.composite) = values(rg)[match(names(rg.composite), values(rg)$gene_sym), c('gene_sym', 'Uniprot')]
#             saveRDS(rg.composite, cached.path.collapsed)
#           }
#         else
#           {
#             rg.composite = split(tmp, tmp$rg.id)
#             values(rg.composite) = values(rg)[as.numeric(names(rg.composite)), ]
#             saveRDS(rg.composite, cached.path)
#           }
#       }
#     else if (gene.collapse)
#         {
#             rg.composite = readRDS(cached.path.collapsed)
#         }
#     else
#         {
#             rg.composite = readRDS(cached.path)
#         }
#
#
#     if (!is.null(genes))
#       rg.composite = rg.composite[values(rg.composite)$gene_sym %in% genes]
#
#     cmap = list(type = c(genespan = bg.col, exon = cds.col, cds.start = st.col, cds.end = en.col, utr.5 = utr.col, utr.3 = utr.col))
#
#     return(gTrack(rg.composite, col = NA, grl.labelfield = 'gene_sym', gr.labelfield = 'exon_id',
#                      gr.srt.label = gr.srt.label, gr.cex.label = gr.cex.label, labels.suppress.gr = labels.suppress.gr, stack.gap = stack.gap, colormaps = cmap, ...))
#   }


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
  labels.suppress.gr = T,
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

    return(suppressWarnings(gTrack(gencode.composite, col = NA, grl.labelfield = 'id', gr.labelfield = 'exon_number',
                     gr.srt.label = gr.srt.label, cex.label = cex.label,
                     gr.cex.label = gr.cex.label,
                     labels.suppress.gr = labels.suppress.gr,
                     stack.gap = stack.gap, colormaps = cmap, ...)))
  }


# @name track.splice
# @title track.splice
#
# Given set of exons and rna bam (eg from tophat) determines junction and exon read density and returns a gTrack object
# of splicing graph
#
# @param ex GRanges of candidate exons
# @param bam path to indexed RNA seq bam
# @param verbose
# @import Rsamtools
# @export
# @author Marcin Imielinski
# track.splice = function(ex = NULL, region = NULL, bam, verbose = TRUE,
#     infer.exons = FALSE,
#     min.reads = 0, ## only relevant if ex is null or infer.exons = TRUE, here exons are inferred from "N" intervals
#     min.exon.width = 10,
#     max.exon.width = 1000 ## only relevant if ex is null
#     )
#     {
#         if (!is.null(ex))
#             {
#                 ex = gr.stripstrand(ex)
#                 ex.ov = gr.findoverlaps(ex, ex)
#             }
#
#         if (is.null(region))
#             reads = read.bam(bam, intervals = ex, pairs.grl = FALSE, verbose = verbose)
#         else
#             reads = read.bam(bam, intervals = region, pairs.grl = FALSE, verbose = verbose)
#
#         if (verbose)
#             cat(length(reads), 'reads\n')
#
#         if (length(reads)==0)
#             {
#                 if (!is.null(ex))
#                     {
#                         ex$expr = 0
#                         ex$log.expr = log(ex$expr)
#                         return(gTrack(ex, y.field = 'log.expr', name = 'Log Read Density', height = 30))
#                     }
#                 else
#                     return(gTrack(region[c()]))
#             }
#
#         sp.reads = splice.cigar(reads, return.grl = FALSE)
#
#         if (verbose)
#             cat(length(sp.reads), 'spliced fragments\n')
#
#         if (is.null(ex))
#             infer.exons = TRUE
#
#         if (infer.exons)
#             {
#                 sp.reads = sp.reads[seqnames(sp.reads)==seqnames(region) & ranges(sp.reads) %over% ranges(region), ]
#                 tmp = grdt(sp.reads)
#                 m.reads = tmp[type != 'N']
#                 n.reads = tmp[type == 'N']
#                 ustart = c(min(start(sp.reads)), n.reads[, length(seqnames), by = end][V1>=min.reads, end]+1)
#                 uend = c(n.reads[, length(seqnames), by = start][V1>=min.reads, start]-1, max(end(sp.reads))) ## 1 after N ends are starts of exons
#
#                 ustart.maxwidth = m.reads[, max(end-start), keyby = start]
#                 uend.maxwidth = m.reads[, max(end-start), keyby = end]
#
#                 new.ex = data.table(seqnames = n.reads$seqnames[1], start = rep(ustart, length(uend)), end = rep(uend, each = length(ustart)))[ (end-start) <= max.exon.width & (end-start) >= min.exon.width, ]
#
#                 ## we will get exon overload unless we cull a bit
#                 ## to be parsimonious we only keep enough exons for all the matching parts of reads to "land on"
#                 ## ie if they start in an exon they will end up in the same one
#                 ## rather than computing all the overlaps necessary for this
#                 ## we just approximate by removing exons for which a smaller one exists that accomodates all of the
#                 ## reads that start on its start or end on its end
#                 new.ex[ , maxwidth.start := (end-start) > ustart.maxwidth[list(new.ex$start), V1]]
#                 new.ex[ , maxwidth.end := (end-start) > uend.maxwidth[list(new.ex$end), V1]]
#
#                 new.ex[ , shorter.start.exists := !(1:length(seqnames) %in% which.min(end-start)), by = start]
#                 new.ex[ , shorter.end.exists := !(1:length(seqnames) %in% which.min(end-start)), by = end]
#
#                 new.ex = seg2gr(new.ex[!shorter.start.exists | !shorter.end.exists, ], seqlengths = seqlengths(sp.reads))[, c()]
#
#                 new.ex$type = 'inferred'
#                 if (verbose)
#                     cat(sprintf('inferred %s total unique exons with min.reads %s, min.exon.width %s, and max.exon.width %s\n', length(new.ex), min.reads, min.exon.width, max.exon.width))
#                 if (!is.null(ex))
#                     {
#                         old.ex = ex;
#                         ex = sort(unique(gUtils::grbind(ex, new.ex)))
#                         if (verbose)
#                             cat(sprintf('Added %s additional exons to yield %s total\n', length(ex)-length(old.ex), length(ex)))
#                     }
#                 else
#                    ex = sort(new.ex)
#
#                 ex.ov = gr.findoverlaps(ex, ex)
#             }
#
#         sp.reads = sort(sp.reads[sp.reads$type != 'N'])
#         sp.reads = sp.reads[order(sp.reads$rid), ]
#
#         tmp =  grdt(sp.reads)
#         tmp[, spid := 1:length(sp.reads)]
#         tmp[, riid := 1:length(seqnames), by = rid]
#         setkey(tmp, spid)
#         sp.reads$riid = tmp[list(1:length(sp.reads)), riid]
#
#         ov = gr.findoverlaps(ex, sp.reads, scol = c('rid', 'riid'), verbose = verbose, return.type = 'data.table')
#         if (length(ov)==0)
#             {
#                 ex$expr = 0
#                 ex$log.expr = log(ex$expr)
#                 return(gTrack(ex, y.field = 'log.expr', name = 'Log Read Density', height = 30))
#             }
#         ov[, width := end-start]
#         ov[, match.left.exon := as.numeric(start == start(ex)[query.id])]
#         ov[, match.right.exon := as.numeric(end == end(ex)[query.id])]
#         ov[, match.left.read := as.numeric(start == start(sp.reads)[subject.id])]
#         ov[, match.right.read := as.numeric(end == end(sp.reads)[subject.id])]
#         ov = ov[(match.left.read | match.left.exon) & (match.right.read | match.right.exon), ]
#         setkeyv(ov, c('rid', 'riid'))
#
#         sp.rid = unique(ov$rid[ov$riid>1])
#         ij = ov[list(sp.rid), list(from = rep(query.id, each = length(query.id)), to = rep(query.id, length(query.id)),
#             val = width[rep(1:length(query.id), each = length(query.id))] + width[rep(1:length(query.id), length(query.id))],
#             from.riid = rep(riid*match.right.exon*match.right.read, each = length(query.id)),
#             to.riid = rep(riid*match.left.exon*match.left.read, length(query.id))), by = rid]
#
#         ij = ij[(ij$to.riid - ij$from.riid) == 1 & ij$to.riid !=0 & ij$from.riid != 0, ]
#         setkeyv(ij, c('from', 'to'))
#         ij = ij[!list(ex.ov$query.id, ex.ov$subject.id), ]
#
# #        edges = ij[, list(val = sum(val)/(width(ex)[from] + width(ex)[to])), keyby = list(from, to)]
#         edges = ij[, list(val = sum(val)), keyby = list(from, to)]
#                                         #edges = edges[, val := 0]
#         edges = edges[, lwd := affine.map(log(val+1), c(0, 6), cap = TRUE)]
#         edges = edges[, v := 5]
#         edges = edges[, h := 2]
#         edges = edges[, col := alpha('gray10', affine.map(log(val+1), c(0.1,0.8), cap = TRUE))]
#         edges = edges[, cex.arrow := 0]
#         ex$expr = ov[(match.right.exon*match.right.read + match.left.exon*match.right.read) | (match.left.read==1 & match.right.read==1), sum(width), keyby = query.id][list(1:length(ex)), V1]/width(ex)
#
#         ex$log.expr = round(log10(ex$expr), 1)
#         ex$ywid = 0.8
#         ex$border = 'black'
#         ex$col = alpha('blue', 0.4)
#         td.ex = gTrack(ex, y.field = 'log.expr', edges = edges, name = 'Log Read Density', height = 30)
#         return(td.ex)
#     }

##########################
#' @name draw.ranges
#' @title draw.ranges
#' @description
#'
#' Draws (genomic, Iranges, Rle, or df with start and end field) intervals
#' as rectangles of thickness "lwd", col "col".
#'
#' If "col" or "y" are set to null then will examine "x" to see if the respective fields are contained within (e.g. as RangedData or data.frame)
#'
#' Adds optional backbone between beginning of first range
#' and end of last on existing plot
#' @keywords internal
#' @author Marcin Imielinski
draw.ranges = function(x, y = NULL, lwd = 0.5, col = "black", border = col, label = NULL,
                       draw.backbone = F, # if TRUE lines will be drawn linking all ranges at a single y-level +/- (a single y x group level if group is specified),
                       group = NULL, # a vector of group labels same size as x (char or numeric), non-NA groups with more than 1 member will be connected by backbone (if draw.backbone = T)
                       col.backbone = col[1],
                       cex.label = 1,
                       srt.label = 45,
                       adj.label = c(0.5, 0),
                       angle, # vector of angles (0 = rectangle -45 leftward barb, +45 rightward barb
                       circles = FALSE, # if TRUE will draw circles at range midpoints instead of polygons
                       bars = FALSE, # if TRUE will draw bars from y0.bar to the y locations
                       points = NA, # if integer then will draw points at the midpoint with pch value "points"
                       lines = FALSE, # if TRUE will connect endpoints of ranges with lines
                       y0.bar = NULL, # only relevant if bars = T, will default to pars('usr')[3] if not specified
                       strand.dir = T,  # if so, will interpret strand field as "direction sign" for interval, + -> right, - -> left, only makes difference if draw.pencils = T
                       xlim = NULL, # only relevant if new.plot = T, can be ranges object
                       ylim = NULL, # only relevant if new.plot = T
                       clip = NULL, # GRanges / IRanges or 1 row df with $start, $end denoting single region to clip to .. will not draw outside of this region

                       lwd.border = 1,
                       lwd.backbone = 1,
                       verbose=FALSE,
                       ...)
{
  PAD = 0 ##0.5;

  if (is.null(y))
    #if (is.data.frame(x))
    #  y = stack_segs(x)
    #else
    y = IRanges::disjointBins(IRanges::ranges(x))

  #if (inherits(x, 'Rle'))
  #  x = RangedData(x)

  if (inherits(x, 'GappedAlignments'))
    x = as(x, 'GRanges')

  if (inherits(x, 'Ranges') | inherits(x, 'GRanges') | inherits(x, 'RangedData'))
    x = standardize_segs(x)

  if (nrow(x)==0)
    return()

  if (!is.null(y))
    x$y = y;

  x$group = group;

  if (!is.null(col) & is.null(x$col))
    x$col = col;

  if (!is.null(border) & is.null(x$border))
    x$border = border;

  if (!is.null(lwd.border) & is.null(x$lwd.border))
    x$lwd.border = lwd.border;

  if (!is.null(label) & is.null(x$label))
    x$label = label;

  if (!is.null(circles) & is.null(x$circles))
    x$circles = circles

  if (!is.null(bars) & is.null(x$bars))
    x$bars = bars

  if (is.null(x$points))
    x$points = NA

  if (!is.na(points) & all(is.na(x$points)))
    x$points = points

  if (!is.null(clip))
  {
    clip = standardize_segs(clip);
    x = clip.seg(x, clip)
  }

  if (is.null(srt.label))
    srt.label = NA

  if (is.na(srt.label))
    srt.label = 45


  if (strand.dir & !is.null(strand))
    x$direction = x$strand

  x$lwd = lwd

  if (nrow(x)>0)
  {
    x$group = paste(x$group, x$y);

    if (draw.backbone)
    {
      gcount = table(x$group);
      x2 = x[x$group %in% names(gcount)[gcount>1], ]

      if (nrow(x2)>0)
      {
        b.st = aggregate(formula = pos1 ~ group, data = x2, FUN = min)
        b.en = aggregate(formula = pos2 ~ group, data = x2, FUN = max)
        b.y = aggregate(formula = y ~ group, data = x2, FUN = function(x) x[1])

        seg.coord = data.frame(x0 = b.st[,2], x1 = b.en[,2], y = b.y[,2]);

        if (!is.null(clip))
        {
          seg.coord$x0 = pmax(seg.coord$x0, clip$pos1)
          seg.coord$x1 = pmin(seg.coord$x1, clip$pos2)
        }

        if (verbose)
          print('Finished computing backbones')
        segments(seg.coord$x0, seg.coord$y, seg.coord$x1, seg.coord$y, col = col.backbone, lwd = lwd.backbone);
      }
    }

    if (is.null(x$direction))
      x$direction = "*"

    if (is.na(circles))
      circles = F

    if (is.na(lines))
      lines = F


    bars = ifelse(is.na(bars), FALSE, bars)
    circles = ifelse(is.na(bars), FALSE, circles)
    lines = ifelse(is.na(bars), FALSE, lines)

    
    if (!circles & !lines & !bars)
    {
      main.args = list(x0 = x$pos1 - PAD, y0 = x$y-x$lwd/2, x1 = x$pos2 + PAD, y1 = x$y + x$lwd/2, col = x$col, border = x$border , angle = angle*c('*'=0, '+'=1, '-'=-1)[x$direction], lwd = x$lwd.border)
      other.args = list(...);
      other.args = other.args[setdiff(names(other.args), names(main.args))]
      do.call(barbs, c(main.args, other.args))
    }
    else
    {
      if (!is.na(points))
      {
        points((x$pos1 + x$pos2)/2, x$y-x$lwd/2, pch = points, border = x$border, col = x$col, cex = x$lwd.border)
      }

      if (circles)
      {
        points((x$pos1 + x$pos2)/2, x$y-x$lwd/2, pch = points, col = x$col, cex = x$lwd.border)
        points((x$pos1 + x$pos2)/2, x$y-x$lwd/2, pch = 1, col = x$border, cex = x$lwd.border)
      }

      if (bars)
      {
        if (is.null(y0.bar))
          y0.bar = par('usr')[3]

        rect(x$pos1, rep(y0.bar, nrow(x)), x$pos2, x$y, col = x$col, border = x$col, lwd = x$lwd.border/2)
        #              rect(x$pos1, rep(y0.bar, nrow(x)), x$pos2, x$y, col = x$col, border = NA, lwd = NA)
      }

      if (lines)
        lines(as.vector(t(cbind(x$pos1, x$pos2))), as.vector(t(cbind(x$y, x$y))), col = as.vector(t(cbind(x$col, x$col))),
              lwd = as.vector(t(cbind(x$lwd.border, x$lwd.border))))

      if (any(!is.na(x$points)))
      {
        points((x$pos1 + x$pos2)/2, x$y-x$lwd/2, pch = x$points, col = x$col, cex = x$lwd.border)
      }

      if (circles)
      {
        points((x$pos1 + x$pos2)/2, x$y-x$lwd/2, pch = 19, col = x$col, cex = x$lwd.border)
        points((x$pos1 + x$pos2)/2, x$y-x$lwd/2, pch = 1, col = x$border, cex = x$lwd.border)
      }
    }

    if (!is.null(x$label))
    {
      has.label = !is.na(x$label)
      if (any(has.label))
      {
        if (nrow(x)>1)
          OFFSET = min(diff(sort(x$y)))/2
        else
          OFFSET = diff(par('usr')[3:4])/100

        if (adj.label[2] == 1)
            {
                text(rowMeans(x[has.label, c('pos1', 'pos2')]), x$y[has.label]+x$lwd[has.label]/2+OFFSET, as.character(x$label[has.label]),
                     adj = adj.label, cex = 0.5*cex.label, srt = srt.label)
            }
        else if (adj.label[2] == 0.5)
            {            
                text(rowMeans(x[has.label, c('pos1', 'pos2')]), x$y[has.label], as.character(x$label[has.label]),
                     adj = adj.label, cex = 0.5*cex.label, srt = srt.label)            
            }
        else
            {
                text(rowMeans(x[has.label, c('pos1', 'pos2')]), x$y[has.label]-x$lwd[has.label]/2-OFFSET, as.character(x$label[has.label]),
                     adj = adj.label, cex = 0.5*cex.label, srt = srt.label)
            }
    }
  }
} 
}

#' @name draw.grl
#' @title  draw.grl
#' @description
#'
#' wrapper function around draw.ranges to draw GRangesList across a set of genomic windows
#' By default, these windows correspond to the smallest set of "covered" regions (ie containing at least one
#' grl in their span, but can also be specified as an argument.
#'
#' basically same idea as draw.granges, except can more compactly draw / stack paired / group reads or
#' gene models / alignments and deal with "grouped" intervals
#'
#' additional values features of input grl that can impact drawing:
#' ywid = vertical thickness of y segments
#' border.col = color border (can override colormap)
#'
#' if optional arg draw.var == T
#' will either call varbase on grl (e.g. if grl is a GRanges or GappedAlignment object)
#' to get variant ranges or will take vars object (GRangesList that is of same length as grl)
#' and will display variant ranges according to color scheme in arg "var.col"
#' "var.col" is a named vector that maps variant types to colors, with the following syntax:
#'
#' XA, XT, XG, XC --> colors of base substitutions
#' S, D, I --> colors of sooutft clipped, deletions, insertions
#' e.g.
#' var.col = c(XA = 'blue', S = 'yellow', 'I' = green') would replace the default variant coloring of
#' "A" substitutions, soft clipped regions, and insertions with blue, yellow, and green colors, respectively.
#'
#' @keywords internal
#' @author Marcin IMielinski
#' @importFrom GenomicRanges GRanges values ranges width strand values<- strand<- seqnames coverage ranges<-
#' @importFrom data.table := setkeyv
#' @importFrom GenomeInfoDb Seqinfo seqinfo keepSeqlevels seqlevels seqlengths seqlevels<- seqlengths<- genome<- seqnames
draw.grl = function(grl,
                    y = NULL,  # can be either vector of length grl, or data.frame row / list with fields $start and $end
                    # specifying y coordinates to "pack" the granges into (or just length 2 list)
                    # note this is different from ylim, which determines the size of the canvas
                    ylim = NULL,  # if null and y provided, will find range of y and plot
                    ywid = NULL,
                    edges = NULL, ## data.frame specifying edges connecting items of grl, with fields $from $to indexing grl and optional fields $lwd, $col, $lty specifying formatting options for connectors, for gr the from arrow will leave the right side of a + range and left side of a - range, and enter the left side of a + range and right side of a - range.   For grl, from arrows will be drawn from the last range in the "from" list item to the first range in the "to" list item
                    draw.paths = F, # if draw.paths = T will treat grl's as paths / contigs,
                    # connecting intervals joined by arrowed curved arcs using alternate stacking algo (ignoring y information)

                    draw.var = F, # if true, then varbase will be called on grl or unlist(grl)
                    var = NULL, # optional GRangesList of same length as grl specifying variant bases with value cols $type, $varbase (ie output of varbase)
                    # corresponding to each grl item
                    var.col = NULL, # named vector with variant colors to override - valid names are XA, XT, XG, XC, S, D, I
                    var.soft = T, ## whether to draw soft clips for varbase
                    windows,  # windows specifies windows, can have optional meta data features $col and $border
                    win.gap = NULL,
                    stack.gap,
                    min.gapwidth = 1, # only applies if windows is not specified
                    col = NULL, border = NA, # all of these variables can be scalar or vectors of length(grl),
                    # can be over-ridden by values / columns in grl
                    col.backbone = alpha('gray', 0.8),
                    gr.colorfield, ## values field in the gr from which colors can be mapped
                    gr.colormap = NULL, ## named vector mapping fields in the gr.colorfield to colors, if unspecified brewer.master() will be applied
                    gr.labelfield = NULL, ## field of gr labels to draw.
                    grl.labelfield = NULL, ## field of grl to draw as label
                    leg.params,
                    labels = NULL, # vector of length(grl)
                    labels.suppress = F,
                    labels.suppress.grl = labels.suppress,
                    labels.suppress.gr = labels.suppress,
                    spc.label = 0.05, # number between 0 and 1 indicating spacing of label
                    adj.label = c(0, 1),
                    cex.label = 1,
                    gr.cex.label = 0.8 * cex.label,
                    gr.srt.label = 0,
                    gr.adj.label = c(0,0.5),
                    new.plot, new.axis, sep.draw, sep.lty, sep.lwd,
                    sep.bg.col,
                    y.pad,  # this is the fractional padding to put on top and bottom of ranges if y is specified as $start and $end pair (def was 0.05)
                    xaxis.prefix, xaxis.suffix, xaxis.unit, xaxis.round, xaxis.interval, xaxis.pos,
                    xaxis.pos.label, xaxis.cex.label, xaxis.newline, xaxis.chronly, xaxis.ticklen,
                    xaxis.width, xaxis.label.angle, xaxis.cex.tick,
                    ylim.grid = ylim, # can be also specified directly for plots with multiple axes and/or non-numeric tracks
                    y.grid = NA, # if non NA, then the number specifies the spacing between y grid (drawn from ylim[1] to ylim[2]), can also be named vector mapping grid.tick labels to plot locations
                    y.grid.col = alpha('gray', 0.5),
                    y.grid.pretty = 5,
                    y.grid.cex = 1,
                    y.grid.lty = 2,
                    y.grid.lwd = 1,
                    path.col = 'black',
                    path.col.arrow = 'black',
                    path.cex.arrow = 1,
                    path.stack.y.gap = 1,
                    path.stack.x.gap = 0,
                    path.cex.v = 1,
                    path.cex.h = 1,
                    draw.backbone = NULL,
                    xlim = c(0, 20000), # xlim of canvas
                    points = NA, ## if non NA then will draw a given point with pch style
                    circles = F, ## only one of these should be true, however if multiple are true then they will be interpreted in this order
                    bars = F,
                    y0.bar = NULL,
                    lines = F,
                    angle, # angle of barbs to indicate directionality of ranges
                    verbose=FALSE,
                    triangle=FALSE, # trigger a triangle matrix plot
                    ylim.parent=NULL, ## ylim of the full thing. This is importat for angle preseveration
                    legend.params = list(plot=TRUE),
                    bg.col = 'white', ## background of whole plot
                    ...)
{
  now = Sys.time();
  ylim.subplot = NULL
  empty.plot = FALSE

  if (length(grl)>0)
  {

    if (is.null(draw.backbone))
      draw.backbone = TRUE

    # if (inherits(grl, 'GappedAlignments'))
    #   grl = ga2gr(grl)

    if ((inherits(grl, 'GRanges')))
    {
      if (!is.null(windows)) ## hack to get over stupid GRanges speed issues when we have a giant GRanges input (eg coverage)
      {
        strand(windows) <- rep("*", length(windows)) ## don't need strand on windows, mess up intersections

        ix <- gUtils::gr.in(grl, windows)
        grl = grl[ix]

        if (!is.null(col))
          if (length(col)!=1)
            col = col[ix]

        if (!is.null(border))
          if (length(border)!=1)
            border = border[ix]

        if (!is.null(ywid))
          if (length(ywid)!=1)
            ywid = ywid[ix]

        if (!is.null(y))
          if (!is.list(y))
            if (length(y)!=1)
              y = y[ix]

        if (!is.null(labels))
          if (!is.list(labels))
            if (length(labels)!=1)
              labels = labels[ix]

        if (!is.null(edges))
          if (nrow(edges)>0 & all(c('from', 'to') %in% colnames(edges)))
          {
            if (data.table::is.data.table(edges))
              edges = as.data.frame(edges)

            ix.i = which(ix)
            edges = edges[edges$from %in% ix.i | edges$to %in% ix.i,];
            edges$from = match(edges$from, ix.i)
            edges$to = match(edges$to, ix.i)
          }
      }
      gr.names = names(grl);
      names(grl) = NULL;

      if (length(grl)>0)
        grl = GenomicRanges::split(grl, 1:length(grl))
      else
        grl = GenomicRanges::GRangesList(grl)

      names(grl) = gr.names;
    }

    if (is.null(labels))
      if (is.null(values(grl)$labels))
        labels = names(grl)
      else
        labels = values(grl)$labels
    
      # make sure names are unique
      names(grl) = 1:length(grl);

      if (is.na(labels.suppress.grl))
        labels.suppress.grl = FALSE

      if (is.na(labels.suppress.gr))
        labels.suppress.gr = FALSE

      if (is.na(draw.var))
        draw.var = F

      if (is.na(draw.paths))
        draw.paths = F

      if (!is.null(gr.colorfield))
        if (is.na(gr.colorfield))
          gr.colorfield = NULL


      if (!is.null(gr.colormap))
        if (all(is.na(gr.colormap)))
          gr.colormap = NULL
      else if (is.null(names(gr.colormap)))
        gr.colormap = NULL

      if (!is.null(col))
        if (all(is.na(col)))
          col = NULL

      if (is.null(col))
        if (!is.null(values(grl)$col))
          col = values(grl)$col
      else if (is.null(gr.colorfield) & is.null(gr.colormap))
        col = 'black'
      else
        col = NA

      if (is.na(col.backbone))
        col.backbone = NULL

      if (is.null(col.backbone))
        col.backbone = col

      if (is.na(border))
        border = NULL

      if (is.null(border))
        border = col

      grl.props = cbind(data.frame(group = names(grl), col = col, stringsAsFactors = F), as.data.frame(values(grl)))
      grl.props$border = border;
      grl.props$ywid = ywid;

      if (!is.null(y) & !is.list(y) & length(y)>0) ## if y coordinate is specified for the ranges
      {
        grl.props$y = y;
        bad = is.na(y)
        bad[which(is.infinite(y))] = TRUE
        if (any(bad))
        {
          grl.props = grl.props[!bad, ]
          grl = grl[!bad]

          if (!is.null(labels))
            labels = labels[!bad]

          y = grl.props$y
        }
      }
      if (!is.null(grl.labelfield))
      {
        if (!is.na(grl.labelfield))
            {
                if (grl.labelfield %in% names(grl.props))
                    grl.props$grl.labels = grl.props[, grl.labelfield]
            }        
        else if (!is.null(labels))
            grl.props$grl.labels = labels ## use $grl.labels to allow labeling of individual grs
    }
      else if (!is.null(labels))
        if (!is.na(labels[1])) ## will only be null if labels is NULL and names(grl) was NULL
          grl.props$grl.labels = labels ## use $grl.labels to allow labeling of individual grs

      gr = tryCatch(grl.unlist(grl), error = function(e)
      {
        ## ghetto solution if we get GRanges names snafu
        gr = unlist(grl);
        c = textConnection(names(gr));        
        cat('budget .. \n')
        if (length(gr)>0)
        {
          gr$grl.ix = read.delim(c, sep = '.', header = F)[,1];
          #            gr$grl.iix = levapply(rep(1, length(gr)), gr$grl.ix, FUN = function(x) 1:length(x))
          gr$grl.iix = data.table::data.table(ix = gr$grl.ix)[, iix := 1:length(ix), by = ix][, iix]
        }
        close(c)
        return(gr)
      })

      gr$group = grl.props$group[gr$grl.ix]
      gr$group.ord = gr$grl.iix
      gr$first = gr$grl.iix == 1

      last = iix = NULL ## NOTE fix
      if (length(gr)>0)
        gr$last = data.table::data.table(iix = as.numeric(gr$grl.iix), ix = gr$grl.ix)[, last := iix == max(iix), by = ix][, last]

      #          gr$last = levapply(gr$grl.iix, gr$grl.ix, FUN = function(x) x == max(x))

      grl.props$group = as.character(grl.props$group)

      S4Vectors::values(gr) = cbind(as.data.frame(values(gr)),
                                    grl.props[match(values(gr)$group, grl.props$group), setdiff(colnames(grl.props), c(colnames(values(gr)), 'group', 'labels')), drop = FALSE])

      #####
      # variant drawing
      ####

    if (draw.var & is.null(var) & length(gr)>0)
        {
                                        #        var = bamUtils::varbase(gr[gr.in(gr, windows)], soft = var.soft)
            gix = which(gr.in(gr, windows))
            var = varbase(gr[gix], soft = var.soft)
            if (any(iix <- var$type == 'I'))
                end(var[ix]) = end(var[ix])+1
            names(var) = gix            
        }
    else
        gix = NULL
        

    if (!is.null(var))
        if (class(var)=='GRangesList')
      {
        VAR.COL = get.varcol()

        if (!is.null(var.col))
            VAR.COL[names(var.col)] = var.col;        

        names(var) = NULL

        if (!is.null(gix))
            names(var) = gix
        else
            names(var) = 1:length(var)
        
                                        #        var.group = as.numeric(as.data.frame(var)$element)
        var.gr = grl.unlist(var)
        if (length(var.gr)>0)
            {
                var.group = names(var)[var.gr$grl.ix]
                                        #            var.gr = gr.stripstrand(unlist(var))


                                        # inherit properties from gr

                values(var.gr) = cbind(as.data.frame(values(var.gr)),
                          as.data.frame(values(gr)[as.numeric(var.group), setdiff(names(values(gr)), c('labels'))]))
                
                if (!is.null(gr.labelfield))
                    if (!is.na(gr.labelfield))
                        values(var.gr)[, gr.labelfield] = NA

                var.gr$col.sig = as.character(var.gr$type);
                xix = var.gr$type == 'X'
                var.gr$col.sig[xix] = paste(var.gr$col.sig[xix], var.gr$varbase[xix], sep = '')
                var.gr$col = VAR.COL[var.gr$col.sig]
                var.gr$border = var.gr$col
                var.gr$first = FALSE
                var.gr$last = FALSE

                if (draw.paths) ## if we are drawing paths, then need to setdiff variant vs non variant parts of edges and re-order
                    {
                        ## remove soft clips
                        var.gr = var.gr[!(var.gr$type %in% c('H',  'S'))]
                        gr2 = gr;
                        GenomicRanges::strand(var.gr) == GenomicRanges::strand(gr)[var.gr$grl.ix]

                        gr$grl.iix = as.numeric(gr$grl.iix)

                        ## now doing some ranges acrobatics to find all the pieces of gr that are not in var.gr
                        ## TODO: speed up, remove awkawardness
                        ir = IRanges::ranges(gr)
                        var.ir = IRanges::ranges(var.gr)
                        tmp.ix = split(1:length(var.gr), var.gr$grl.ix)
                        tmp.l = lapply(names(tmp.ix), function(i) IRanges::disjoin(c(ir[as.numeric(i)], var.ir[tmp.ix[[i]]])))
                        tmp.ogix = rep(as.numeric(names(tmp.ix)), sapply(tmp.l, length))
                        tmp.ir = do.call(c, tmp.l)
                        tmp.gr = GenomicRanges::GRanges(seqnames(gr)[tmp.ogix], tmp.ir, seqlengths = GenomeInfoDb::seqlengths(gr), og.ix = tmp.ogix)
                        tmp.ov = gr.findoverlaps(tmp.gr, var.gr)
                        tmp.ov = tmp.ov[tmp.gr$og.ix[tmp.ov$query.id] == var.gr$grl.ix[tmp.ov$subject.id] ]                   
                        new.gr = tmp.gr[!(1:length(tmp.gr) %in% tmp.ov$query.id), ] ## only keep the non variant pieces
                        GenomicRanges::strand(new.gr) = GenomicRanges::strand(gr)[new.gr$og.ix]
                        values(new.gr) = cbind(as.data.frame(values(gr)[new.gr$og.ix, ]) , og.ix = new.gr$og.ix)                            
                        var.gr$og.ix = var.gr$grl.ix
                        GenomicRanges::strand(var.gr) = GenomicRanges::strand(gr)[var.gr$og.ix]
                                        #          var.gr$group = as.numeric(as.character(gr$group[var.gr$og.ix]))

                        new.gr = grbind(new.gr, var.gr, gr[setdiff(1:length(gr), var.gr$grl.ix)])
                        new.gr$grl.iix = as.numeric(gr$grl.iix[new.gr$og.ix])
                        new.ord = mapply(function(x, y, z) if (y[1]) x[order(z)] else rev(x[order(z)]),
                            split(1:length(new.gr), new.gr$og.ix), split(as.logical(GenomicRanges::strand(new.gr)=='+'), new.gr$og.ix), split(start(new.gr), new.gr$og.ix))
                        new.ix = unlist(lapply(new.ord, function(x) ((1:length(x))-1)/length(x)))
                        new.gr$grl.iix[unlist(new.ord)] = new.gr$grl.iix[unlist(new.ord)] + new.ix

                        gr = new.gr[order(new.gr$group, new.gr$grl.iix), ]
                    }
            }

        else
        {
          gr = grbind(gr, var.gr)
        }
      }

      if (length(gr)>0)
        names(gr) = 1:length(gr)

      if (is.null(windows)) ## find covered windows in provided grl
      {
        seqlevels(gr) = seqlevels(gr)[seqlevels(gr) %in% as.character(seqnames(gr))]
        windows = as(coverage(gr), 'GRanges');
        windows = windows[values(windows)$score!=0]
        windows = reduce(windows, min.gapwidth = min.gapwidth);
      }

      else if (!is(windows, 'GRanges'))
        if (is(windows, 'GRangesList'))
          windows = unlist(windows)
      else  ## assume it's a seqinfo object or an object that has a seq
        windows = si2gr(windows)

      if (is.null(win.gap))
        win.gap = mean(width(windows))*0.2

      if (sum(as.numeric(width(windows)))==0)
        stop('Windows have width 0')

      if (verbose) {
        print('Before flatmap')
        print(Sys.time() - now)
    }
    
    ## add 1 bp to end for visualization .. ranges avoids weird width < 0 error
    if (length(gr)>0)
        {
            IRanges::ranges(gr) = IRanges::IRanges(start(gr), pmax(end(gr), pmin(end(gr)+1, GenomeInfoDb::seqlengths(gr)[as.character(seqnames(gr))], na.rm = T), na.rm = T)) ## jeremiah commented
                                        #        end(gr) = pmax(end(gr), pmin(end(gr)+1, seqlengths(gr)[as.character(seqnames(gr))], na.rm = T), na.rm = T)
        }

      suppressWarnings(end(windows) <- end(windows) + 1) ## shift one needed bc gr.flatmap has continuous convention, we have categorical (e.g. 2 bases is width 2, not 1)
      mapped = gr.flatmap(gr, windows, win.gap);
  
      grl.segs = mapped$grl.segs;
      window.segs = mapped$window.segs;

      if (verbose) {
        print('After flatmap')
        print(Sys.time() - now)
      }

      grl.segs$border = as.character(grl.segs$border)
      grl.segs$col = as.character(grl.segs$col)
      grl.segs$group = as.character(grl.segs$group)
      grl.segs$strand = as.character(grl.segs$strand)

      if (!is.null(gr.labelfield))
        if (!is.na(gr.labelfield))
          if (gr.labelfield %in% names(grl.segs))
            grl.segs$label = grl.segs[, gr.labelfield]

      if (nrow(grl.segs)==0)
      {
        #           warning('No ranges intersecting window');
        #          return()
      }

      if (is.list(y))
        if (all(c('start', 'end') %in% names(y)))
          ylim.subplot = c(y$start[1], y$end[1])
      else if (length(y)>1)
        ylim.subplot = c(y[[1]], y[[length(y)]])

      if (nrow(grl.segs)>0)
      {
        if (!draw.paths)
        {
          if (is.null(y) | !is.null(ylim.subplot))
          {
            pos1  = aggregate(formula = pos1 ~ group, data = grl.segs, FUN = min);
            pos2  = aggregate(formula = pos2 ~ group, data = grl.segs, FUN = max);
            pos1 = structure(pos1[,2], names = pos1[,1]) - round(stack.gap/2);
            pos2 = structure(pos2[,2]-1, names = pos2[,1]) + round(stack.gap/2);

            ix = order(as.numeric(names(pos1)))
            pos1 = pos1[ix]
            pos2 = pos2[ix]

            ## FIX HERE .. avoid using disjointBins
            ## these gymnastics allow us to use disjointBins (which works on IRanges)
            ## without getting integer overflow
            if (max(c(pos1, pos2))>2e9)
            {
              m = max(c(pos1,pos2));
              pos1 = ceiling(pos1/m*2e9)
              pos2 = floor(pos2/m*2e9);
            }

            # bin ranges
            y.bin = IRanges::disjointBins(IRanges::IRanges(pos1, pos2))

            m.y.bin = max(y.bin)
            if (is.null(ylim))
              ylim = c(1, m.y.bin) + c(-0.5*m.y.bin, 0.5*m.y.bin)

            ## squeeze y coordinates into provided (or inferred) ylim
            if (!is.null(ylim.subplot))
              tmp.ylim = ylim.subplot
            else
              tmp.ylim = ylim

            ## provide bottom and top padding of y.bin
            y.pad = max(c(0, min(0.49, y.pad)));
            #                    tmp.ylim = tmp.ylim + c(1, -1)*y.pad*diff(tmp.ylim);
            y.pad = pmin(1/(m.y.bin+1)/2, 0.125)
            tmp.ylim = tmp.ylim + c(1, -1)*y.pad*diff(tmp.ylim);

            y = structure(affine.map(y.bin, tmp.ylim), names = names(pos1));

            grl.segs$y = y[grl.segs$group]
          }
          else
          {
            if (is.null(ylim))
              if (any(!is.na(y[!is.infinite(y)])))
              {
                tmp.ylim = range(y[!is.infinite(y)], na.rm = T)
                ylim = tmp.ylim + c(-1, 1)*0.2*diff(tmp.ylim) + c(-1, 0.2)*y.pad*diff(tmp.ylim)
              }
            else
              ylim = c(0,10)
          }
        }
        else  ## draw.paths = T -->  will treat each grl as a sequence, which will be joined by connectors
        {
          ix.l = lapply(split(1:nrow(grl.segs), grl.segs$group), function(x) x[order(grl.segs$group.ord[x])])
          grl.segs$y.relbin = NA

          ## we want to layout paths so that we prevent collissions between different paths
          grl.segs$y.relbin[unlist(ix.l)] = unlist(lapply(ix.l, function(ix)
          {
            # find runs where start[i+1]>end[i] and strand[i] == strand[i+1] = '+'
            # and end[i+1]<start[i] and strand[i] == strand[i+1] = '-'
            if (length(ix)>1)
            {
              iix = 1:(length(ix)-1)
              concordant = ((grl.segs$pos1[ix[iix+1]] >= grl.segs$pos2[ix[iix]]
                             & grl.segs$strand[ix[iix+1]] != '-' & grl.segs$strand[ix[iix]] != '-') |
                              (grl.segs$pos1[ix[iix+1]] <= grl.segs$pos2[ix[iix]]
                               & grl.segs$strand[ix[iix+1]] == '-' & grl.segs$strand[ix[iix]] == '-'))
              return(c(0, cumsum(!concordant)))
            }
            else
              return(0)
          }))

          contig.lim = data.frame(
            group = names(vaggregate(formula = y.relbin ~ group, data = grl.segs, FUN = max)),
            pos1  = vaggregate(formula = pos1 ~ group, data = grl.segs, FUN = min) - round(stack.gap)/2,
            pos2  = vaggregate(formula = pos2~ group, data = grl.segs, FUN = max) + round(stack.gap)/2,
            height = vaggregate(formula = y.relbin ~ group, data = grl.segs, FUN = max)
          );
          contig.lim$width = contig.lim$pos2 - contig.lim$pos1
          contig.lim$y.bin = 0;

          contig.lim = contig.lim[order(-contig.lim$width), ]

          if (nrow(contig.lim)>1)
            for (i in 2:nrow(contig.lim))
            {

              # find lowest level at which there is no clash with this and previously stacked segments


              ir1 = IRanges::IRanges(contig.lim[1:(i-1), 'pos1'], contig.lim[1:(i-1), 'pos2'])
              ir2 = IRanges::IRanges(contig.lim[i, 'pos1'], contig.lim[i, 'pos2'])
#              clash = which(gUtils::gr.in(ir1, ir2 + path.stack.x.gap))
#              clash = which(gUtils::gr.in(ir1, ir2 + path.stack.x.gap))
              clash = which(ir1 %over% (ir2 + path.stack.x.gap))
              pick = clash[which.max(contig.lim$y.bin[clash] + contig.lim$height[clash])]
              contig.lim$y.bin[i] = c(contig.lim$y.bin[pick] + contig.lim$height[pick] + path.stack.y.gap, 0)[1]
            }

          grl.segs$y.bin = contig.lim$y.bin[match(grl.segs$group, contig.lim$group)] + grl.segs$y.relbin + 1

          m.y.bin = max(grl.segs$y.bin)
          if (is.null(ylim))
            ylim = c(1, m.y.bin) + c(-0.5*m.y.bin, 0.5*m.y.bin)

          ## squeeze y coordinates into provided (or inferred) ylim
          if (!is.null(ylim.subplot))
            tmp.ylim = ylim.subplot
          else
            tmp.ylim = ylim

          ## provide bottom and top padding of y.bin
          #
          y.pad = 1/(m.y.bin+1)/2
          y.pad = pmin(1/(m.y.bin+1)/2, 0.125)
          tmp.ylim = tmp.ylim + c(1, -1)*y.pad*diff(tmp.ylim);

          ## make final y coordinates by squeezing y.bin into tmp.ylim
          grl.segs$y = affine.map(grl.segs$y.bin, tmp.ylim)

        }
        if (verbose) {
          print('After agg')
          print(Sys.time() - now)
        }              #        xlim = c(min(window.segs$start), max(window.segs$end));
      }
      else
        empty.plot = TRUE

      #window.segs$end <- window.segs$end + 1 ##debug
      ## now collapse everything to 0, 1 based on windows.segs
      winlim = range(c(window.segs$start, window.segs$end))

      ## xlim here is arbitrary, just needs to be > 1, helps us scale overlayed plots to window.segs boundaries (if many plots with different
      ## windows are drawn on the same canvas
      grl.segs$pos1 = affine.map(grl.segs$pos1, winlim, ylim = xlim)
      grl.segs$pos2 = affine.map(grl.segs$pos2, winlim, ylim = xlim)
      window.segs$start = affine.map(window.segs$start, winlim, ylim = xlim)
      window.segs$end = affine.map(window.segs$end, winlim, ylim = xlim)
  }
  else
    empty.plot = TRUE

  if (empty.plot)
  {
    if (is.null(windows))
      stop('Either nonempty range data or windows must be provided')

    mapped = gr.flatmap(GRanges(), windows, win.gap)
    window.segs = mapped$window.segs
    winlim = range(c(window.segs$start, window.segs$end))
    window.segs$start = affine.map(window.segs$start, winlim, ylim = xlim)
    window.segs$end = affine.map(window.segs$end, winlim, ylim = xlim)
    #        xlim = c(min(window.segs$start), max(window.segs$end));

    if (is.null(ylim))
      ylim = c(0, 1)

    if (is.list(y) & is.null(ylim.subplot))
      if (all(c('start', 'end') %in% names(y)))
        ylim.subplot = c(y$start[1], y$end[1])
    else
      ylim.subplot = c(y[[1]], y[[2]])

    if (is.null(ylim.subplot))
      ylim.subplot = ylim
  }

  # if new plot will add (optional) axes and vertical lines separating windows
  if (new.plot)
  {
    if (verbose) {
      print('Before axis draw')
      print(Sys.time() - now)
    }
    plot.blank(xlim = xlim, ylim = ylim, bg.col = bg.col);
    new.axis = TRUE
  }

  if (new.plot || new.axis)
  {
    if (is.null(xaxis.pos)) {
      if (!is.null(ylim.subplot))
        xaxis.pos = ylim.subplot[1]-0.05*diff(ylim.subplot)
      else
        xaxis.pos = ylim[1]+0.12*diff(ylim)
    }

    if (is.null(window.segs$col))
      window.segs$col = sep.bg.col

    if (is.null(xaxis.pos.label)) {
      if (!is.null(ylim.subplot))
        xaxis.pos.label = xaxis.pos - 0.04*diff(ylim.subplot)
      else
        xaxis.pos.label = xaxis.pos - 0.04*diff(ylim)
    }

        if (sep.draw && length(windows)>1)
        {
          #rect(window.segs$end[1:(nrow(window.segs)-1)], rep(ylim[1], nrow(window.segs)-1),
          #window.segs$start[2:(nrow(window.segs))], rep(ylim[2], nrow(window.segs)-1), border = 'white', col = sep.col)
          if (any(width(windows)<=0))
            warning('Some windows are width 0')

          sep.loc = c(window.segs$start, window.segs$end);
          if (!is.null(ylim.subplot)) ## limit separator drawing to actual data limits
          {

            if (is.null(window.segs$border))
              window.segs$border = 'white'

            sep.x0 = window.segs$start[1:(nrow(window.segs))]
            sep.x1 = window.segs$end[1:(nrow(window.segs))]
            sep.y0 = rep(xaxis.pos, nrow(window.segs))
            bgcol.l <- as.character(window.segs$col) ## JEREMIAH
            bgborder.l <- as.character(window.segs$border) ## JEREMIAH
            #rep(min(xaxis.pos.label, xaxis.pos, ylim.subplot[1]), nrow(window.segs))
            sep.y1 = rep(ylim.subplot[2], nrow(window.segs))
            rect(sep.x0, sep.y0, sep.x1, sep.y1, border = bgborder.l, col = bgcol.l) ## JEREMIAH added bgcol.l
            segments(sep.x0, sep.y0, sep.x0, sep.y1, lty = sep.lty, lwd = sep.lwd)
            segments(sep.x1, sep.y0, sep.x1, sep.y1, lty = sep.lty, lwd = sep.lwd)
          }
          else
          {
            rect(window.segs$start[1:(nrow(window.segs))], rep(ylim[1], nrow(window.segs)),
                 #window.segs$end[1:(nrow(window.segs))], rep(ylim[2], nrow(window.segs)), border = 'white', col = sep.bg.col) ## MARCIN
                 window.segs$end[1:(nrow(window.segs))], rep(ylim[2], nrow(window.segs)), border = 'white', col = as.character(window.segs$col)) ## JEREMIAH
            abline(v = sep.loc, col = 'gray', lty = sep.lty, lwd = sep.lwd);
          }
        }

        if (new.axis)
        {
          nwin = length(windows);

          ## draw the actual x axis
          segments(window.segs$start, rep(xaxis.pos[1], nwin), window.segs$end, rep(xaxis.pos[1], nwin));

          # if (!is.null(xaxis.suffix))
          #   if (is.na(xaxis.suffix) | nchar(xaxis.suffix)==0)
          #     xaxis.suffix = NULL

          draw_x_ticks(xaxis.interval, windows, mapped, winlim, xlim, ylim, xaxis.pos, xaxis.suffix, xaxis.unit, xaxis.cex.tick, xaxis.ticklen, xaxis.round)

          # then (label) text
          newline <- ifelse(xaxis.newline, '\n', '')

          width.text = ''
          if (xaxis.width)
          {
            if (!is.null(xaxis.suffix))
              width.text = paste('(', paste(prettyNum(ifelse(rep(xaxis.unit == 1, length(windows)),
                                                             width(windows), round(width(windows)/xaxis.unit, 2)), big.mark = ','), xaxis.suffix),  ')', sep = '')
            else
              width.text = paste('(', prettyNum(ifelse(rep(xaxis.unit == 1, length(windows)),
                                                       width(windows), round(width(windows)/xaxis.unit, 2)), big.mark = ','),  ')', sep = '')
          }

          begin.text = prettyNum(pmax(floor(1/xaxis.unit),
                                      ifelse(rep(xaxis.unit == 1, length(windows)), start(windows), round(start(windows)/xaxis.unit, xaxis.round))),
                                 big.mark = ',')

          end.text = prettyNum(ifelse(rep(xaxis.unit == 1, length(windows)), end(windows),
                                      round(end(windows)/xaxis.unit, xaxis.round)), big.mark = ',')

          if (!xaxis.chronly) {
            text(rowMeans(window.segs[, c('start', 'end')]), rep(xaxis.pos.label, nwin),
                 paste(xaxis.prefix, ' ',  seqnames(windows), ':',newline,
                       begin.text,'-', newline,
                       end.text, ' ', xaxis.suffix, newline, width.text, sep = ''),
                 cex = xaxis.cex.label*0.8, srt = 0, adj = c(0.5, 0), srt=xaxis.label.angle)
          } else {
            text(rowMeans(window.segs[, c('start', 'end')]), rep(xaxis.pos.label, nwin),
                 paste(xaxis.prefix, ' ',  seqnames(windows),
                       sep = ''),
                 cex = xaxis.cex.label*0.8, srt = 0, adj = c(0.5, 0), srt=xaxis.label.angle)
          }
        }
  }

  if (empty.plot)
  {
    if (verbose) {
      print('Returning ..')
      print(Sys.time() - now)
    }
    return(window.segs)
  }

  line.loc = NULL
  if (!is.na(y.grid[1]))
  {
    if (is.logical(y.grid))
      line.loc = pretty(ylim.grid,y.grid.pretty)
    else if (length(y.grid)==1) ## only interval is specified
      line.loc = seq(floor(ylim.grid[1]/y.grid)*y.grid, ceiling(ylim.grid[2]/y.grid)*y.grid, y.grid)
    else ## specific grid lines are specified
      line.loc = y.grid

    if (is.null(names(line.loc)))
      names(line.loc) = line.loc;

    #        abline(h = line.loc, col = y.grid.col, lty = y.grid.lty, lwd = y.grid.lwd)
    segments(xlim[1], line.loc,  xlim[2], line.loc, col = y.grid.col, lty = y.grid.lty, lwd = y.grid.lwd)

    if (is.null(y.grid.cex))
      y.grid.cex = NA

    if (is.na(y.grid.cex))
      y.grid.cex = 1

    axis(2, at = line.loc, labels = names(line.loc), tick = TRUE, pos = line.loc[1], las = 2, cex.axis = y.grid.cex)
  }

  if (length(grl)>0)
  {
    if (is.null(grl.segs$ywid))
      grl.segs$ywid = 1

    if (any(nix <- is.na(grl.segs$ywid)))
      grl.segs$ywid[nix] = 1

    if (!is.null(ylim.subplot))
      tmp.ydiff = diff(ylim.subplot)*(1-2*y.pad)
    else if (!is.null(line.loc))
      tmp.ydiff = diff(range(line.loc))
    else
      tmp.ydiff = diff(range(grl.segs$y, na.rm = T))

    if (tmp.ydiff==0)
      tmp.ydiff = ylim

    fact = 1.5*(1+length(unique(grl.segs$y)))
    if (!is.null(line.loc))
      fact = pmax(10, pmin(1000, fact))

    #        if (!is.null(line.loc))
    #            grl.segs$ywid = pmin(tmp.ydiff / fact, min(diff(line.loc))/2) * grl.segs$ywid
    #            grl.segs$ywid = pmin(tmp.ydiff / fact, min(diff(line.loc))/2) * grl.segs$ywid
    #        else

    ywid = tmp.ydiff / fact

    if (!is.null(line.loc))
      ywid = pmin(ywid, min(c(1, diff(sort(unique(grl.segs$y))))*2)) ## want to be able to see a minimal difference between data points

    if (!is.null(line.loc)) ## don't want segments to be fatter than a grid unit
      ywid = pmin(min(diff(line.loc))/2, ywid)

    grl.segs$ywid  = ywid * grl.segs$ywid

    #            grl.segs$ywid[na.ix] = min((max(grl.segs$y[na.ix])-min(grl.segs$y[na.ix]))/length(unique(grl.segs$y))/1.25, diff(ylim.subplot)/5)
    #          }

    if (is.null(gr.colorfield))
        gr.colorfield = NA

    
    if (gr.colorfield %in% names(grl.segs))
        {            
            if (is.null(gr.colormap))
                {
                    uval = unique(as.character(grl.segs[, gr.colorfield]))
                    gr.colormap = structure(alpha(brewer.master(length(uval)), 0.5), names = uval)
                }
            
            cols = gr.colormap[as.character(grl.segs[, gr.colorfield])];
            grl.segs$col[!is.na(cols)] = cols[!is.na(cols)]
            
            if (is.null(grl.segs$border))
                grl.segs$border[!is.na(cols)] = cols[!is.na(cols)]
            else if (any(is.na(grl.segs$border))) ## only override border if unspecified
                {
                    ix = !is.na(cols) & is.na(grl.segs$border)
                    grl.segs$border[ix] = cols[ix]
                }
        }
            
    if (leg.params$plot && length(gr.colormap)>0)
        {
            leg.params$x = leg.params$xpos * diff(xlim) + xlim[1]
            leg.params$y = leg.params$ypos * diff(par('usr')[3:4]) + par('usr')[3]
            leg.params$legend = names(gr.colormap)
            if (circles) {
                leg.params$col = gr.colormap
                leg.params$pch = 16
            }
            else
                leg.params$fill = gr.colormap
            leg.params$border = gr.colormap
            leg.params$xpos = leg.params$ypos = NULL

                                        # if (length(gr.colormap)>legend.maxitems & legend.maxitems > 0)
                                        #   gr.colormap = gr.colormap[intersect(names(sort(table(grl.segs[, gr.colorfield]), decreasing = T)[1:legend.maxitems]),
                                        #     names(gr.colormap))]

            do.call(graphics::legend, leg.params)
                                        #if (circles)
                                        #    graphics::legend(legend.pos[1]*diff(xlim) + xlim[1], legend.pos[2]*diff(par('usr')[3:4]) + par('usr')[3], legend = names(gr.colormap), col = gr.colormap, border = gr.colormap, cex = legend.cex * 0.5, ncol = legend.ncol, xjust = legend.xjust, pch = 16, yjust = legend.yjust)
                                        #else
                                        #    graphics::legend(legend.pos[1]*diff(xlim) + xlim[1], legend.pos[2]*diff(par('usr')[3:4]) + par('usr')[3], legend = names(gr.colormap), fill = gr.colormap, border = gr.colormap, cex = legend.cex * 0.5, ncol = legend.ncol, xjust = legend.xjust, yjust = legend.yjust)
        }

    if (draw.paths)
      draw.backbone = FALSE

    if (labels.suppress.gr)
      grl.segs$label = NULL

    if (!is.na(lines))
      if (lines)
        grl.segs = grl.segs[order(grl.segs$pos1), ]

    draw.ranges(grl.segs, y = grl.segs$y, group = grl.segs$group, col = grl.segs$col, border = grl.segs$border, ylim = ylim, xlim = xlim, lwd = grl.segs$ywid, draw.backbone = draw.backbone, angle = angle, col.backbone = col.backbone, points = points, circles = circles, bars = bars, y0.bar = y0.bar, lines = lines, cex.label = gr.cex.label, srt.label = gr.srt.label, adj.label = gr.adj.label, ...)

    ## if draw.paths, will now draw connectors
    ##
    if (draw.paths)
    {
      grl.segs = grl.segs[order(grl.segs$grl.ix, grl.segs$grl.iix), ]
      ix.l = split(1:nrow(grl.segs), grl.segs$group)
      grl.segs$ctype = NA;  ## connector type

      # keep track to see if next index is missing --> will allow us to draw dotted connectors for these ..
      # "missing" indices occur when we are focusing on windows and potentially removing some of the
      # ranges in the contig
      grl.segs$next.missing = !((grl.segs$query.id+1) %in% grl.segs$query.id) & !grl.segs$last
      grl.segs$prev.missing = !((grl.segs$query.id-1) %in% grl.segs$query.id) & !grl.segs$first

      if (is.null(grl.segs$is.cycle))
        grl.segs$is.cycle = FALSE;

      connector.args = do.call('rbind', lapply(ix.l, function(ix)
      {
        # find runs where start[i+1]>end[i] and strand[i] == strand[i+1] = '+'
        # and end[i+1]<start[i] and strand[i] == strand[i+1] = '-'
        out = NULL;
        if (length(ix)>1)
        {
          out = data.frame(ix0 = ix[-length(ix)], ix1 = ix[-1], type = 'U', sign = '+', cyclic = F, stringsAsFactors = F);
          discordant = grl.segs$y.bin[ix[-length(ix)]] != grl.segs$y.bin[ix[-1]]
          out$type[discordant & grl.segs$strand[ix[-1]] == grl.segs$strand[ix[-length(ix)]]] = 'S'
          out$type[discordant & grl.segs$strand[ix[-1]] != grl.segs$strand[ix[-length(ix)]]] = 'S'
        }

        if (grl.segs$next.missing[ix[length(ix)]]) ## make bridge to nowhere
          out = rbind(out, data.frame(ix0 = ix[length(ix)], ix1 = NA, type = 'S', cyclic = F,
                                      sign = grl.segs$strand[ix[length(ix)]], stringsAsFactors = F))
        if (grl.segs$prev.missing[ix[1]]) ## make bridge from nowhere
          out = rbind(out, data.frame(ix0 = NA, ix1 = ix[1], type = 'S', cyclic = F,
                                      sign = grl.segs$strand[ix[1]], stringsAsFactors = F))

        if (grl.segs$is.cycle[ix[length(ix)]])
        {
          if (grl.segs$y.bin[ix[length(ix)]] == grl.segs$y.bin[ix[1]])
            out = rbind(out, data.frame(ix0 = ix[length(ix)], ix1 = ix[1], type = 'U', cyclic = T, sign = '-', stringsAsFactors = F))
          else if (grl.segs$strand[ix[length(ix)]] == '-' & grl.segs$strand[ix[1]] == '+')
            out = rbind(out, data.frame(ix0 = ix[length(ix)], ix1 = ix[1], type = 'S', cyclic = T, sign = '+', stringsAsFactors = F))
          else if (grl.segs$strand[ix[length(ix)]] == '+' & grl.segs$strand[ix[1]] == '-')
            out = rbind(out, data.frame(ix0 = ix[length(ix)], ix1 = ix[1], type = 'S', cyclic = T, sign = '-', stringsAsFactors = F))
          else if (grl.segs$strand[ix[length(ix)]] == '-' & grl.segs$strand[ix[1]] == '-')
            out = rbind(out, data.frame(ix0 = ix[length(ix)], ix1 = ix[1], type = 'S', cyclic = T, sign = '+', stringsAsFactors = F))
          else if (grl.segs$strand[ix[length(ix)]] == '+' & grl.segs$strand[ix[1]] == '+')
            out = rbind(out, data.frame(ix0 = ix[length(ix)], ix1 = ix[1], type = 'S', cyclic = T, sign = '-', stringsAsFactors = F))
        }
        return(out)
      }))

      if (!is.null(connector.args))
      {
        path.h = path.cex.h * rep(diff(par('usr')[1:2])/50, nrow(connector.args))
        if (any(connector.args$type == 'S'))
          path.h[connector.args$type == 'S'] = 2*path.h[connector.args$type == 'S']

        path.v = rep(path.cex.v, nrow(connector.args))
        path.v[is.na(connector.args$ix0) | is.na(connector.args$ix1)] = path.cex.v*2*grl.segs$ywid[connector.args$ix0[is.na(connector.args$ix0) | is.na(connector.args$ix1)]]

        #                path.v[!is.na(connector.args$ix0)] = path.cex.v*2*grl.segs$ywid[connector.args$ix0[!is.na(connector.args$ix0)]]
        #               path.v[!is.na(connector.args$ix1)] = pmax(path.v[!is.na(connector.args$ix1)],
        #                      path.cex.v*2*grl.segs$ywid[connector.args$ix1[!is.na(connector.args$ix1)]])

        connector.args$strand0[!is.na(connector.args$ix0)] = grl.segs$strand[connector.args$ix0[!is.na(connector.args$ix0)]]
        connector.args$strand1[!is.na(connector.args$ix1)] = grl.segs$strand[connector.args$ix1[!is.na(connector.args$ix1)]]
        connector.args$strand0[is.na(connector.args$ix0)] = grl.segs$strand[connector.args$ix1[is.na(connector.args$ix0)]]
        connector.args$strand1[is.na(connector.args$ix1)] = grl.segs$strand[connector.args$ix0[is.na(connector.args$ix1)]]

        connector.args$strand0 = ifelse(connector.args$strand0 == '*', '+', connector.args$strand0)
        connector.args$strand1 = ifelse(connector.args$strand1 == '*', '+', connector.args$strand1)

        connector.args$x0[!is.na(connector.args$ix0) & connector.args$strand0=='+'] =
          grl.segs$pos2[connector.args$ix0[!is.na(connector.args$ix0) & connector.args$strand0=='+']]
        connector.args$x0[!is.na(connector.args$ix0) & connector.args$strand0=='-'] =
          grl.segs$pos1[connector.args$ix0[!is.na(connector.args$ix0) & connector.args$strand0=='-']]
        connector.args$x1[!is.na(connector.args$ix1) & connector.args$strand1=='+'] =
          grl.segs$pos1[connector.args$ix1[!is.na(connector.args$ix1) & connector.args$strand1=='+']]
        connector.args$x1[!is.na(connector.args$ix1) & connector.args$strand1=='-'] =
          grl.segs$pos2[connector.args$ix1[!is.na(connector.args$ix1) & connector.args$strand1=='-']]

        # bridges from  nowhere
        connector.args$x0[is.na(connector.args$ix0) & connector.args$strand1=='+'] =
          grl.segs$pos1[connector.args$ix1[is.na(connector.args$ix0) & connector.args$strand1=='+']] - path.h[is.na(connector.args$ix0) & connector.args$strand1=='+']
        connector.args$x0[is.na(connector.args$ix0) & connector.args$strand1=='-'] =
          grl.segs$pos2[connector.args$ix1[is.na(connector.args$ix0) & connector.args$strand1=='-']] + path.h[is.na(connector.args$ix0) & connector.args$strand1=='-']

        # bridges to nowhere
        connector.args$x1[is.na(connector.args$ix1) & connector.args$strand0=='+'] =
          grl.segs$pos2[connector.args$ix0[is.na(connector.args$ix1) & connector.args$strand0=='+']] + path.h[is.na(connector.args$ix1) & connector.args$strand0=='+']
        connector.args$x1[is.na(connector.args$ix1) & connector.args$strand0=='-'] =
          grl.segs$pos1[connector.args$ix0[is.na(connector.args$ix1) & connector.args$strand0=='-']] - path.h[is.na(connector.args$ix1) & connector.args$strand0=='-']

        connector.args$y0[!is.na(connector.args$ix0)] = grl.segs$y[connector.args$ix0[!is.na(connector.args$ix0)]]
        connector.args$y1[!is.na(connector.args$ix1)] = grl.segs$y[connector.args$ix1[!is.na(connector.args$ix1)]]
        connector.args$y0[is.na(connector.args$ix0)] = grl.segs$y[connector.args$ix1[is.na(connector.args$ix0)]] - 0.25*path.v[is.na(connector.args$ix0)]
        connector.args$y1[is.na(connector.args$ix1)] = grl.segs$y[connector.args$ix0[is.na(connector.args$ix1)]] + 0.25*path.v[is.na(connector.args$ix1)]

        connector.args$lty = 1;

        connector.args$lty[grl.segs$next.missing[connector.args$ix0[!is.na(connector.args$ix0)]] &
                             !connector.args$cyclic[!is.na(connector.args$ix0)]] = 3
        connector.args$lty[grl.segs$prev.missing[connector.args$ix1[!is.na(connector.args$ix1)]] &
                             connector.args$cyclic[!is.na(connector.args$ix1)]] = 3
        connector.args$lty[is.na(connector.args$ix0) | is.na(connector.args$ix1)] = 3; ## label all bridge to / from nowhere with dotted line

        ##
        path.v[connector.args$y0 == connector.args$y1] = 0
        path.h[connector.args$y0 == connector.args$y1] = 0
        path.v[connector.args$cyclic] = path.cex.v*2*grl.segs$ywid[connector.args$ix1[connector.args$cyclic]]
        path.h[connector.args$cyclic] = path.cex.h * diff(par('usr')[1:2])/20

        ## workaround current lines() limitation in connectors
        if (any(lty3 <- connector.args$lty == 3))
          connectors(connector.args$x0[lty3], connector.args$y0[lty3], connector.args$strand0[lty3],
                     connector.args$x1[lty3], connector.args$y1[lty3], ifelse(connector.args$strand1[lty3] == '+', '-', '+'),
                     type = connector.args$type[lty3],
                     h = path.h[lty3], v = path.v[lty3],
                     lty = connector.args$lty[lty3], col = path.col, col.arrow = path.col.arrow,
                     cex.arrow = grl.segs$ywid[1]*path.cex.arrow, f.arrow = T)

        if (any(lty1 <- connector.args$lty == 1))
          connectors(connector.args$x0[lty1], connector.args$y0[lty1], connector.args$strand0[lty1],
                     connector.args$x1[lty1], connector.args$y1[lty1], ifelse(connector.args$strand1[lty1] == '+', '-', '+'),
                     type = connector.args$type[lty1],
                     h = path.h[lty1], v = path.v[lty1],
                     lty = connector.args$lty[lty1], col = path.col, col.arrow = path.col.arrow,
                     cex.arrow = grl.segs$ywid[1]*path.cex.arrow, f.arrow = T)


        #text(connector.args$x1, connector.args$y1, paste(connector.args$ix0, ' (', grl.segs$group[connector.args$ix0], '), ', connector.args$ix1, ' (', grl.segs$group[connector.args$ix1], '), ', connector.args$type, ' ', connector.args$sign, sep = ''))
        #            text(connector.args$x1, connector.args$y1-path.v, paste(grl.segs$group[connector.args$ix0], connector.args$type, ' ', connector.args$sign, sep = ''), )
      }
    }

    if (!is.null(edges))
      if (nrow(edges)>0 & all(c('from', 'to') %in% colnames(edges)))
      {
        if (data.table::is.data.table(edges))
          edges = as.data.frame(edges)

        if (is.null(edges$col))
          edges$col = NA

        if (is.null(edges$lty))
          edges$lty = NA

        if (is.null(edges$lwd))
          edges$lwd = NA

        if (is.null(edges$h))
          edges$h = 1

        if (is.null(edges$cex.arrow))
          edges$cex.arrow = 1

        edges.og = edges
        edges = edges[edges$to %in% grl.segs$group | edges$from %in% grl.segs$group, ]

        if (nrow(edges)>0)
        {
          if (any(ix <- is.na(edges$col)))
            edges$col[ix] = 'black'

          if (any(ix <- is.na(edges$lwd)))
            edges$lwd[ix] = 1

          if (any(ix <- is.na(edges$lty)))
            edges$lty[ix] = 1

          ## now translate $from $to from group to gr id, and choosing last range in group for "from" indices
          ## and first range in group for "to" indices putting NA's when gr not in grl.segs

          first.ix = which(grl.segs$first)
          last.ix = which(grl.segs$last)

          grl.segs$to.gr = grl.segs$from.gr = 1:nrow(grl.segs)
          edges$edge.id = 1:nrow(edges)
          edges = merge(merge(edges, grl.segs[last.ix, c('group', 'from.gr')], by.x = 'from', by.y = 'group', all.x = T),
                        grl.segs[first.ix, c('group', 'to.gr')], by.x = 'to', by.y = 'group', all.x = T)

          
          #                  edges$to.gr = first.ix[match(edges$to, grl.segs[first.ix, ]$group)]
          #                  edges$from.gr = last.ix[match(edges$from, grl.segs[last.ix, ]$group)]

          ## replicate edges that connect input gr's that have multiple instantiations

          ## figure out which grs are clipped
          grl.segs$clipped.start = !(start(gr)[grl.segs$query.id] == grl.segs$start)
          grl.segs$clipped.end = !(end(gr)[grl.segs$query.id] == grl.segs$end)

          clipped.to = ((grl.segs$clipped.start[edges$to.gr] & grl.segs$strand[edges$to.gr] == '+') |
                          (grl.segs$clipped.end[edges$to.gr] & grl.segs$strand[edges$to.gr] == '-'))
          clipped.from = ((grl.segs$clipped.start[edges$from.gr] & grl.segs$strand[edges$from.gr] == '-') |
                            (grl.segs$clipped.end[edges$from.gr] & grl.segs$strand[edges$from.gr] == '+'))

          clipped.to[is.na(clipped.to)] = F
          clipped.from[is.na(clipped.from)] = F
          ## now determine start and end points based on signs of connectors
          to.ix = !is.na(edges$to.gr) & !clipped.to
          from.ix = !is.na(edges$from.gr) & !clipped.from
          edges$x.pos.from = edges$x.pos.to = edges$y.pos.from = edges$y.pos.to = edges$dir.from = edges$dir.to =  NA

          #                  from.ix = !is.na(edges$fromx.gr)
          #                  to.ix = !is.na(edges$to.gr)

          if (any(to.ix)) ## connect to left end if + and right end if -
          {
            edges$x.pos.to[to.ix] = ifelse(grl.segs$strand[edges$to.gr[to.ix]] != '-',
                                           grl.segs$pos1[edges$to.gr[to.ix]], grl.segs$pos2[edges$to.gr[to.ix]])
            edges$y.pos.to[to.ix] = grl.segs$y[edges$to.gr[to.ix]]
            edges$dir.to[to.ix] = ifelse(grl.segs$strand[edges$to.gr[to.ix]] != '-',
                                         '-', '+')
          }

          if (any(from.ix)) ## connect from right end if + and left end if -
          {
            edges$x.pos.from[from.ix] = ifelse(grl.segs$strand[edges$from.gr[from.ix]] != '-',
                                               grl.segs$pos2[edges$from.gr[from.ix]], grl.segs$pos1[edges$from.gr[from.ix]])
            edges$y.pos.from[from.ix] = grl.segs$y[edges$from.gr[from.ix]]
            edges$dir.from[from.ix] = ifelse(grl.segs$strand[edges$from.gr[from.ix]] != '-',
                                             '+', '-')
          }


          ## now let NA endpoints "dangle"
          dangle.w = diff(par('usr')[1:2])/100
          edges$dangle = F

          if (is.null(edges$dangle.w))
            edges$dangle.w = 1

          edges$dangle.w  = edges$dangle.w * dangle.w

          if (any(!from.ix))
          {
              edges$x.pos.from[!from.ix] =
                  edges$x.pos.to[!from.ix] + ifelse(grl.segs$strand[edges$to.gr[!from.ix]] != '-', -1, 1)*edges$dangle.w[!from.ix]
            edges$y.pos.from[!from.ix] = edges$y.pos.to[!from.ix]
            edges$from.gr[!from.ix] =  edges$from.gr[!from.ix]
            edges$dir.from[!from.ix] =  ifelse(edges$dir.to[!from.ix] != '-', '-', '+')
            edges$dangle[!from.ix] = T
          }

          if (any(!to.ix))
          {
              edges$x.pos.to[!to.ix] =
                  edges$x.pos.from[!to.ix] + ifelse(grl.segs$strand[edges$from.gr[!to.ix]] != '-', 1, -1)*edges$dangle.w[!to.ix]
            edges$y.pos.to[!to.ix] = edges$y.pos.from[!to.ix]
            edges$to.gr[!to.ix] =  edges$to.gr[!to.ix]
            edges$dir.to[!to.ix] =  ifelse(edges$dir.from[!to.ix] != '-', '-', '+')
            edges$dangle[!to.ix] = T
          }

          edges$from.ix = from.ix
          edges$to.ix = to.ix
          edges = edges[order(edges$dangle), ]; ## hack to make sure we prefer non-dangling versions of edges when we have a choice (TODO: MAKE NON HACKY)
          dup = duplicated(cbind(edges$edge.id, edges$from.gr)) | duplicated(cbind(edges$edge.id, edges$to.gr))
          edges = edges[(edges$from.ix | edges$to.ix) & !dup, ]

          if (nrow(edges)>0)
          {

            if (!is.null(ylim.subplot))
              tmp.ydiff = diff(ylim.subplot)*(1-2*y.pad)
            else
              tmp.ydiff = diff(range(grl.segs$y, na.rm = T))

            #                      edges$type = ifelse(edges$y.pos.from == edges$y.pos.to, 'U', 'S')
            uthresh = max(grl.segs$ywid)
            edges$type = ifelse(abs(edges$y.pos.from - edges$y.pos.to) <= uthresh, 'U', 'S')

            edges$f.arrow = T
            edges$b.arrow = F
            if (is.null(edges$v))
              edges$v = 1


            edges$v = ifelse(abs(edges$y.pos.from - edges$y.pos.to)<=uthresh,
                             2*uthresh*affine.map(abs(edges$x.pos.to-edges$x.pos.from), xlim = c(0, diff(par('usr')[1:2])), ylim = c(0.5, 1.0)),
                             abs(edges$y.pos.from-edges$y.pos.to)/2)*edges$v
            #                  edges$v[is.na(edges$v)] = 0
            edges$h = dangle.w*edges$h * affine.map(abs(edges$x.pos.to-edges$x.pos.from), xlim = c(0, diff(par('usr')[1:2])), ylim = c(0.5, 1.0))
            make.flat = edges$type == 'U' & edges$dir.to != edges$dir.from

            if (!is.null(edges$not.flat)) ## allow edges specified as $not.flat to have height, unless dangling
              if (any(ix <- edges$not.flat & !edges$dangle))
                make.flat[ix] = F

            if (any( ix <- make.flat ))
            {
              edges$v[ix] = 0
              edges$h[ix] = 0
            }

            if (!is.null(edges$v.override))
              if (any(ix <- (!is.na(edges$v.override) & edges$to.ix & edges$from.ix)))
                edges$v[ix] = edges$v.override[ix]

            if (!is.null(edges$h.override))
              if (any(ix <- (!is.na(edges$h.override) & edges$to.ix & edges$from.ix)))
                edges$h[ix] = dangle.w*edges$h.override[ix]


            connectors(edges$x.pos.from, edges$y.pos.from, edges$dir.from,
                       edges$x.pos.to, edges$y.pos.to, edges$dir.to,
                       v = edges$v, h = edges$h, type = edges$type,
                       f.arrow = T, b.arrow = F,
                       cex.arrow = 0.2*edges$cex.arrow,
                       col.arrow = edges$col,
                       lwd = edges$lwd, lty = edges$lty, col = edges$col)

            if (!is.null(edges$label))
            {
              text((edges$x.pos.from + edges$x.pos.to)/2, edges$y.pos.from + 0.5*edges$v*sign(edges$y.pos.to - edges$y.pos.from + 0.01),
                   edges$label, adj = c(0.5, 0.5), col = 'black')
            }
          }
        }
      }

    if (verbose) {
      print('After draw')
      print(Sys.time() - now)
    }

    if (nrow(grl.segs)>0 & !is.null(grl.props$grl.labels) & !labels.suppress.grl)
    {
      if (!draw.paths)
      {
        pos1  = vaggregate(formula = pos1 ~ group, data = grl.segs, FUN = min);
        pos2  = vaggregate(formula = pos2 ~ group, data = grl.segs, FUN = max);
        ywid  = vaggregate(formula = ywid ~ group, data = grl.segs, FUN = max);
        y  = vaggregate(formula = y ~ group, data = grl.segs, FUN = mean);
        grl.segs.u = data.frame(group = names(pos1), pos1, pos2, y, ywid);
        grl.segs.u$grl.labels = grl.props$grl.labels[match(grl.segs.u$group, grl.props$group)]
        rownames(grl.segs.u) = grl.segs.u$group;

        if (adj.label[2] == 0)
          text(grl.segs.u$pos1, grl.segs.u$y+grl.segs.u$ywid + 0.005*diff(ylim),
               grl.segs.u$grl.labels, adj = adj.label, cex = cex.label)
        if  (adj.label[2] == 0.5)
          text(grl.segs.u$pos1-0.005*diff(xlim), grl.segs.u$y,
               grl.segs.u$grl.labels, adj = adj.label, cex = cex.label)
        else
          text(grl.segs.u$pos1, grl.segs.u$y-grl.segs.u$ywid - 0.005*diff(ylim),
               grl.segs.u$grl.labels, adj = adj.label, cex = cex.label)
      }
      else
      {
        pos1  = vaggregate(formula = pos1 ~ group, data = grl.segs, FUN = min);
        pos2  = vaggregate(formula = pos2 ~ group, data = grl.segs, FUN = max);
        y0  = vaggregate(formula = y ~ group, data = grl.segs, FUN = min);
        y1  = vaggregate(formula = y ~ group, data = grl.segs, FUN = max);
        grl.segs.u = data.frame(group = names(pos1), pos1, pos2, y0, y1);
        grl.segs.u$grl.labels = grl.props$grl.labels[match(grl.segs.u$group, grl.props$group)]

        text(adj.label[1]*grl.segs.u$pos1 + (1-adj.label[1])*grl.segs.u$pos2, adj.label[2]*grl.segs.u$y0 + (1-adj.label[2])*grl.segs.u$y1,
             grl.segs.u$grl.labels, adj = adj.label, cex = cex.label)

      }

    }
  }

  return(window.segs)
}


#####################################
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

#' @name connectors
#' @title connectors
#' @description
#'
#' connectors
#'
#' draws bezier connectors between pairs of signed points points a = (x0, y0, s0) and b(x1, y1, s1)
#'
#' "S" connectors will use bezier control points that are located in an intermediate y location
#' between y0 and y1, and "U" connectors will use bezier control points that are located in an
#' intermediate x location between x0 and x1.
#'
#' type can be "U" or "S", if "S", then sign of corresponding v value is ignored (ie control points will always
#' be vertically in between a and b, in addition, neither intermeidate control will never be more than 1/3 of the y distance
#' between a and b
#'
#' The sign of the endpoints refers to the direction from which the connector approaches the point, s0 = 1 (or '+') means from the right
#' and s0 = 0 means from the left, similarly for s1
#'
#' h = horizonal "bulge" of signed connector, ie depending on the sign the extent to which the connector will extend to the left
#' or the right of points an and b
#' v = vertical bulge (only applies to U connectors), how far above max(y0, y1) the U connector will bulge if v>0 or how far below min(y0, y1)
#' the connector will bulge.
#'
#' nsteps determines how many segments used to approximate the curve
#'
#' all args except nsteps can be vecorized
#'
#' @name Marcin Imielinski
#' @keywords internal
connectors = function(x0, y0, s0 = 1, x1, y1, s1 = 1, v = 0.1, h = 0.1,
                      type = "S", f.arrow = T, b.arrow = F, nsteps = 20,
                      cex.arrow = 1, col.arrow = 'black', lwd = 1, lty = 1, col = 'black')

{
  B = rbind(bernsteinp(6, nsteps), NA); # bernstein basis for bezier with 6 control points

  dat = data.frame(x0, y0, s0, x1, y1, s1, type, v, h, f.arrow, b.arrow, cex.arrow, lwd, lty, stringsAsFactors = F, col = col, col.arrow = col.arrow)

  cpoints.names =
    matrix(paste('c', rep(1:ncol(B), each = 2), c('.x', '.y'), sep = ''), nrow = 2, ncol = ncol(B)-2, dimnames = list(c('x', 'y'), 1:(ncol(B)-2)))

  dat[, as.vector(cpoints.names)] = NA

  if (length(dat) == 0)
    return(NULL)

  scon = dat$type == "S"
  ucon = dat$type == "U"

  if (is.numeric(dat$s0))
    dat$s0 = sign(dat$s0)
  else
    dat$s0 = sign((dat$s0 == '+')-0.5)

  if (is.numeric(dat$s1))
    dat$s1 = sign(dat$s1)
  else
    dat$s1 = sign((dat$s1 == '+')-0.5)

  if (any(scon))
  {
    dat.s = dat[scon, , drop = FALSE]

    dat.s$v = pmin(abs(dat.s$v), abs(dat.s$y1-dat.s$y0)/2)
    dat.s$h = abs(dat.s$h)

    dat.s[, cpoints.names['x', 1]] = dat.s$x0 + dat.s$s0*dat.s$h
    dat.s[, cpoints.names['y', 1]] = dat.s$y0

    ## middle control points will form an "S" or "C" depending on s0 and s1
    dat.s[, cpoints.names['x', 2]] = dat.s[, cpoints.names['x', 1]]
    dat.s[, cpoints.names['y', 2]] = dat.s[, cpoints.names['y', 1]] + ifelse(dat.s$y0>dat.s$y1, -1, 1)*dat.s$v

    dat.s[, cpoints.names['x', 4]] = dat.s$x1 + dat.s$s1*dat.s$h
    dat.s[, cpoints.names['y', 4]] = dat.s$y1

    dat.s[, cpoints.names['x', 3]] = dat.s[, cpoints.names['x', 4]]
    dat.s[, cpoints.names['y', 3]] = dat.s[, cpoints.names['y', 4]] + ifelse(dat.s$y0>dat.s$y1, 1, -1)*dat.s$v

    dat[scon, ] = dat.s
  }

  if (any(ucon))
  {
    dat.u = dat[ucon, , drop = FALSE]

    #        dat.u$h = pmin(abs(dat.u$h), abs(dat.u$x1-dat.u$x0)/3)
    dat.u$h = abs(dat.u$h)

    dat.u[, cpoints.names['x', 1]]  = ifelse(dat.u$x0 < dat.u$x1, pmin((dat.u$x0 + dat.u$x1)/2, dat.u$x0 + dat.u$s0*dat.u$h), pmax((dat.u$x0 + dat.u$x1)/2, dat.u$x0 + dat.u$s0*dat.u$h))
    #        dat.u[, cpoints.names['x', 1]] = dat.u$x0 + dat.u$s0*dat.u$h
    dat.u[, cpoints.names['y', 1]] = dat.u$y0

    ## middle control points will form verticle bubble portion of U (if v<0) or upside down U (if v>0)
    dat.u[, cpoints.names['x', 2]] = dat.u$x0
    dat.u[, cpoints.names['y', 2]] = ifelse(dat.u$v>=0, pmax(dat.u$y0, dat.u$y1) + dat.u$v, pmin(dat.u$y0, dat.u$y1) + dat.u$v)

    dat.u[, cpoints.names['x', 3]] = dat.u$x1
    dat.u[, cpoints.names['y', 3]] = dat.u[, cpoints.names['y', 2]]

    dat.u[, cpoints.names['x', 4]] = ifelse(dat.u$x0 < dat.u$x1, pmax((dat.u$x0 + dat.u$x1)/2, dat.u$x1 + dat.u$s1*dat.u$h), pmin((dat.u$x0 + dat.u$x1)/2, dat.u$x1 + dat.u$s1*dat.u$h))
    #        dat.u[, cpoints.names['x', 4]] = dat.u$x1 + dat.u$s1*dat.u$h
    dat.u[, cpoints.names['y', 4]] = dat.u$y1

    dat[ucon, ] = dat.u
  }

  ## draw bezier around control points .. draw separately for each unique combo of lwd, lty, col
  formats = factor(paste(dat$lwd, dat$lty, dat$col))

  lapply(levels(formats), function(lev)
  {
    lev.ix = which(formats == lev)
    x = as.vector(t(as.matrix(dat[lev.ix, c('x0', cpoints.names['x', ], 'x1')])))
    y = as.vector(t(as.matrix(dat[lev.ix, c('y0', cpoints.names['y', ], 'y1')])))
    ix = seq(1, length(x), ncol(B))
    XY = do.call('rbind', lapply(ix, function(i) B %*% cbind(x[i:(i+ncol(B)-1)], y[i:(i+ncol(B)-1)])))
    lines(XY, lwd = dat$lwd[lev.ix], lty = dat$lty[lev.ix], col = as.character(dat$col)[lev.ix])
  })

  x.wid = diff(par('usr')[1:2])
  y.wid = diff(par('usr')[3:4])
  arrow.height = cex.arrow*0.05*y.wid
  x.offset = tan(120)*arrow.height/y.wid*par('pin')[2]*x.wid/par('pin')[1]

  ## add forward arrow heads if any
  if (any(dat$f.arrow))
  {
    tmp = dat$x1[dat$f.arrow] + dat$s1[dat$f.arrow]*x.offset
    xcoord = rbind(dat$x1[dat$f.arrow], tmp, tmp, NA);
    ycoord = rbind(dat$y1[dat$f.arrow], dat$y1[dat$f.arrow] + arrow.height/2, dat$y1[dat$f.arrow] - arrow.height/2, NA);
    polygon(xcoord, ycoord, col = as.character(dat$col.arrow[dat$f.arrow]), border = as.character(dat$col.arrow[dat$f.arrow]))
  }

  if (any(dat$b.arrow))
  {
    tmp = dat$x0[dat$b.arrow] + dat$s0[dat$b.arrow]*x.offset
    xcoord = rbind(dat$x0[dat$b.arrow], tmp, tmp, NA);
    ycoord = rbind(dat$y0[dat$b.arrow], dat$y0[dat$b.arrow] + arrow.height/2, dat$y0[dat$b.arrow] - arrow.height/2, NA);
    polygon(xcoord, ycoord, col = as.character(dat$col.arrow[dat$b.arrow]), border = as.character(dat$col.arrow[dat$b.arrow]))
  }
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
#' @export
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

#' @name draw.triangle
#' @title draw.triangle
#' @description
#' internal draw.triangle function
#'
#' @author Jeremiah Wala
#' @keywords internal
draw.triangle <- function(grl,mdata,y,
                          ylim.parent=NULL,
                          windows = NULL,
                          win.gap = NULL,
                          m.bg.col,
                          sigma = NA, ## if not NA then will blur with a Gaussian filter using a sigma value of this many base pairs
                          col.min='white',
                          col.max='red',
                          gr.colormap = NA,
                          cmap.min, cmap.max, m.sep.lwd,
                          leg.params,
                          min.gapwidth = 1,
                          islog = FALSE) {

  # if (!is.null(gr.colormap)) {
  #   if (!is.null(names(gr.colormap)))
  #     palette.colors = names(gr.colormap)
  # }

  ylim.subplot = NULL

  xlim = c(0, 20000)
  if (is.list(y)) {
    if (all(c('start', 'end') %in% names(y)))
      ylim.subplot = c(y$start[1], y$end[1])
  } else {
    ylim.subplot = c(y[[1]], y[[2]])
  }

  if (inherits(grl, "GRangesList")) {
    gr <- grl.unlist(grl)
  } else {
    gr <- grl
  }

  ##assume grl is always non zero
  if (is.null(mdata)) {
    mdata = matrix(nrow=length(gr), ncol=length(gr), 0)
  } else if (nrow(mdata) != length(gr)) {
    warning('draw.triangle: square matrix should be same dim as gr')
  }

  ## find covered windows in provided grl
  if (is.null(windows)) {
    seqlevels(gr) = seqlevels(gr)[seqlevels(gr) %in% as.character(seqnames(gr))]
    windows = as(coverage(gr), 'GRanges');
    windows = windows[values(windows)$score!=0]
    windows = reduce(windows, min.gapwidth = min.gapwidth);
  }

  if (is.null(win.gap))
    win.gap = mean(width(windows))*0.2

  ## calculate the background
  mapped = gr.flatmap(gr, windows, win.gap);
  grl.segs = mapped$grl.segs;
  window.segs = mapped$window.segs;
  winlim = range(c(window.segs$start, window.segs$end + 1)) ## 3/20/16 added +1
  mdata = as.matrix(mdata)[as.numeric(grl.segs$query.id), as.numeric(grl.segs$query.id)]

  if (!is.na(sigma)) ## MARCIN: if blur use spatstat library to blur matrix for n base pairs but making sure we don't bleed across windows
  {
    uwin = unique(grl.segs$window)
    win.map = split(1:length(grl.segs$window), grl.segs$window)
    uwin.pairs = cbind(rep(uwin, length(uwin)), rep(uwin, each = length(uwin)))
    for (p in 1:nrow(uwin.pairs))
    {
      print(paste('...blurring', p, 'of', nrow(uwin.pairs), 'windows'))
      tmp.i = win.map[[uwin.pairs[p, 1]]]
      tmp.j = win.map[[uwin.pairs[p, 2]]]
      mdata[tmp.i, tmp.j] = as.matrix(spatstat::blur(spatstat::im(mdata[tmp.i, tmp.j], xcol = grl.segs$start[tmp.j], yrow = grl.segs$end[tmp.i]), sigma = sigma))
    }
  }
  if (nrow(mdata) != nrow(grl.segs))
    warning('problem after flatmap. Should have trimmed matrix to len(gr)')

  if (nrow(grl.segs) > 0)
    if (nrow(grl.segs) != nrow(mdata)) {
      warning('problemo')
    }

  ## transform into plot coordinates
  grl.segs$pos1 = affine.map(grl.segs$pos1, winlim, ylim = xlim) #d 4 down
  grl.segs$pos2 = affine.map(grl.segs$pos2 + 1, winlim, ylim = xlim) ## +1 added

  window.segs$start = affine.map(window.segs$start, winlim, ylim = xlim)
  window.segs$end = affine.map(window.segs$end + 1, winlim, ylim = xlim) ## 3/20/16 added +1 to make [5,6][7,8] adjacent, etc

  #####################
  ## send to plotting device
  #####################

  ## should have already called draw.grl!
  ## make the plot
  wid = 20000
  hgt.subp = abs(diff(ylim.subplot))
  asp.ratio <- (par('pin')[2])/(par('pin')[1])
  hgt.plot <- abs(diff(ylim.parent))
  dlim <- c(0, wid * asp.ratio * (hgt.subp / hgt.plot))
  y0 <- dlim[1]
  y1 <- dlim[2]

  ## draw blank background
  rect(xlim[1]-diff(xlim)*0.1, ylim.subplot[1], xlim[2], ylim.subplot[2], border = NA, col = 'white')

  ## draw the background BOXES
  if (nrow(window.segs) > 1) {
    bgx = .all.xpairs(window.segs$start, window.segs$end)
    out <- diamond(bgx[,1], bgx[,2], bgx[,3], bgx[,4], y0, y1)

    ## affine map to local coorinates
    out$y[!is.na(out$y)] <- affine.map(out$y[!is.na(out$y)], xlim=dlim, ylim=ylim.subplot)
    if (m.sep.lwd > 0)
      polygon(out$x, out$y, col=m.bg.col, lwd=m.sep.lwd)
    else
      polygon(out$x, out$y, col=m.bg.col, border=NA)
  }

  ## draw the background triangles
  out = triangle(x1=window.segs$start, x2=window.segs$end, y=y0, y0=y0, y1=y1)

  out$y[!is.na(out$y)] <- affine.map(out$y[!is.na(out$y)], xlim=dlim, ylim=ylim.subplot)
  #out$y[seq(from=2, to=length(out$y), by=4)] <- affine.map(out$y[seq(from=2, to=length(out$y), by=4)], xlim=dlim, ylim=ylim.subplot)
  if (m.sep.lwd > 0)
    polygon(out$x, out$y, col=m.bg.col, lwd=m.sep.lwd)
  else
    polygon(out$x, out$y, col=m.bg.col, border=NA)

  if (nrow(grl.segs) == 0)
    return(window.segs)

  ################
  ## plot the data
  ################

  ## set the color scale
  bgx = .all.xpairs(grl.segs$pos1, grl.segs$pos2)
  col = mdata[matrix(nrow=nrow(bgx), ncol=2, c(bgx[,5], bgx[,6]))]
  out <- diamond(bgx[,1], bgx[,2], bgx[,3], bgx[,4], y0, y1, col)

  ## set the color scale
  if (is.na(cmap.min))
    #      cmap.min = min(mdata)
    cmap.min = quantile(mdata, 0.01)

  if(is.na(cmap.max))
    cmap.max = quantile(mdata, 0.99)


  #cs <- col.scale(seq(cmap.min, cmap.max), val.range=c(cmap.min, cmap.max), col.min=col.min, col.max=col.max)
  if (is.na(gr.colormap))
    gr.colormap = c("light green", "yellow", "orange", "red")
  else if (is.list(gr.colormap))
    gr.colormap = unlist(gr.colormap)

  cs <- colorRampPalette(gr.colormap)(length(seq(cmap.min, cmap.max, by=(cmap.max-cmap.min)/100)))

  ## affine map to local coorinates
  out$y[!is.na(out$y)] <- affine.map(out$y[!is.na(out$y)], xlim=dlim, ylim=ylim.subplot)
  ix.min <- out$col < cmap.min
  ix.max <- out$col > cmap.max
  out$col[ix.min | ix.max] <- NA
  cr <- rep('black', length(out$col))
  cr[!is.na(out$col)] <- cs[floor(affine.map(out$col[!is.na(out$col)], xlim = c(cmap.min, cmap.max), ylim=c(1,100)))]
  cr[ix.min] <- 'white' ## deal with cmap.min and cmap.max
  cr[ix.max] <- 'black' ## deal with cmap.min and cmap.max

  ## plot the triangles part
  polygon(out$x, out$y, col=cr, border=NA)

  ## plot the triangles part
  col = diag(mdata)
  out.t = triangle(grl.segs$pos1, grl.segs$pos2, y0, y0, y1, col=col)
  out.t$y[!is.na(out.t$y)] <-
    affine.map(out.t$y[!is.na(out.t$y)], xlim=dlim, ylim=ylim.subplot)
  ix.min <- out.t$col < cmap.min
  ix.max <- out.t$col > cmap.max
  out.t$col[ix.min | ix.max] <- NA
  cr <- rep('black', length(out.t$col))
  cr[!is.na(out.t$col)] <- cs[floor(affine.map(out.t$col[!is.na(out.t$col)], xlim = c(cmap.min, cmap.max), ylim=c(1,100)))]
  cr[ix.min] <- 'white' ## deal with cmap.min and cmap.max
  cr[ix.max] <- 'black' ## deal with cmap.min and cmap.max
  #cr <- cs[ceiling(out.t$col) - cmap.min + 1]
  polygon(out.t$x, out.t$y, col=cr, border=NA)

  legend.cex = leg.params$cex
  if (is.null(legend.cex))
    legend.cex = 1

  ## plot the legend
  if (!islog)
    txt = format(c(cmap.min, 0.5*(cmap.max-cmap.min) + cmap.min, cmap.max), digits=1)
  else
    txt = format(c(10^cmap.min, 10^(0.3333*((cmap.max-cmap.min) + cmap.min)), 10^(0.66667*((cmap.max-cmap.min) + cmap.min)), 10^(cmap.max)), digits=1)
  color.bar(lut = cs, xpos=0, ypos=ylim.subplot[1] + diff(ylim.subplot)*0.3, width=500, height=diff(ylim.subplot)*0.3,
            text=txt, cex=legend.cex)

  return(window.segs)
}


#' @name clip_polys
#' @title clip_polys
#' @description
#' Clip polygons
#' @author Jeremiah Wala
#' @keywords internal
clip_polys <- function(dt, y0, y1) {

  ## fix global def NOTE
  y2 <- y3 <- y4 <- y5 <- y6 <- NULL
  left <- x3.tmp <- x3 <- y3.tmp <- lean.sign <- NULL

  ## leave early if not necessary
  miny = min(dt[, c(y1, y2, y3, y4, y5, y6)], na.rm = TRUE)
  maxy = max(dt[, c(y1, y2, y3, y4, y5, y6)], na.rm = TRUE)
  ##print(paste('MinY:', miny, 'MaxY:', maxy, "Y0:", y0, 'Y1:', y1))
  if (y0 <= miny && y1 >= maxy)
    return(dt)

  dt[, left := y3 > y6]
  dt[, x3.tmp := x3]
  dt[, y3.tmp := y3]
  dt$x3[dt$left] <- dt$x6[dt$left]
  dt$y3[dt$left] <- dt$y6[dt$left]
  dt$x6[dt$left] <- dt$x3.tmp[dt$left]
  dt$y6[dt$left] <- dt$y3.tmp[dt$left]
  dt[, lean.sign := ifelse(left, -1, 1)]

  ## remove stuff all the way out
  yc1 = y1
  yc0 = y0
  dt <- subset(dt, y1 < yc1 & y4 > yc0)

  #############
  # deal with the top clip
  #############
  ## clip to only bottom corner
  ix <- dt$y3 >= y1
  dt$y3[ix] <- dt$y6[ix] <- y1
  dt$x3[ix] <- dt$x1[ix] - dt$lean.sign[ix] * (y1 - dt$y1[ix])
  dt$x4[ix] <- dt$x5[ix] <- dt$x3[ix]
  dt$y4[ix] <- dt$y5[ix] <- dt$y3[ix]
  dt$x6[ix] <- dt$x1[ix] + dt$lean.sign[ix] * (y1 - dt$y1[ix])
  ## clip off 4,5,6
  ix <- dt$y6 >= y1 & !ix
  dt$y4[ix] <- dt$y5[ix] <- dt$y6[ix] <- y1
  dt$x4[ix] <- dt$x5[ix] <- dt$x3[ix] + dt$lean.sign[ix] * (y1 - dt$y3[ix])
  dt$x6[ix] <- dt$x1[ix] + dt$lean.sign[ix] * (y1 - dt$y1[ix])
  ## clip off 4,5
  ix <- dt$y4 > y1 & !ix
  dt$y4[ix] <- dt$y5[ix] <- y1
  dt$x4[ix] <- dt$x3[ix] + dt$lean.sign[ix] * (y1 - dt$y3[ix])
  dt$x5[ix] <- dt$x6[ix] - dt$lean.sign[ix] * (y1 - dt$y6[ix])

  ##################
  # deal with the bottom clip
  ##################
  ## clip to only top corner
  ix <- dt$y6 <= y0
  dt$x3[ix] <- dt$x3[ix] + dt$lean.sign[ix] * (y0 - dt$y3[ix])
  dt$x6[ix] <- dt$x6[ix] - dt$lean.sign[ix] * (y0 - dt$y6[ix])
  dt$y3[ix] <- dt$y6[ix] <- y0
  dt$y1[ix] <- dt$y2[ix] <- dt$y3[ix]
  dt$x1[ix] <- dt$x2[ix] <- dt$x3[ix]
  ## clip off 1,2,3
  ix <- dt$y3 < y0 & !ix
  dt$x3[ix] <- dt$x3[ix] + dt$lean.sign[ix] * (y0 - dt$y3[ix])
  dt$x1[ix] <- dt$x1[ix] + dt$lean.sign[ix] * (y0 - dt$y1[ix])
  dt$y1[ix] <- dt$y3[ix] <- y0
  dt$y2[ix] <- dt$y1[ix]
  dt$x2[ix] <- dt$x1[ix]
  ## clip off 1,2
  ix <- dt$y1 < y0 & !ix
  dt$x1[ix] <- dt$x1[ix] - dt$lean.sign[ix] * (y0 - dt$y1[ix])
  dt$x2[ix] <- dt$x2[ix] + dt$lean.sign[ix] * (y0 - dt$y2[ix])
  dt$y1[ix] <- dt$y2[ix] <- y0

  return(dt)

}

#' @name diamond
#' @title diamond
#' @description
#' internal draw.triangle function
#' @author Jeremiah Wala
#' @keywords internal
diamond <- function (x11, x12, x21, x22, y0, y1, col=NULL) {


  i1 = i2 = .geti(x12, x21)
  i3 = .geti(x11, x21)
  i4 = i5 = .geti(x11, x22)
  i6 = .geti(x12, x22)

  ## dummy
  if (is.null(col))
    col <- i1$x

  ## setup the data table
  dt <- data.table::data.table(x1=i1$x, x2=i2$x, x3=i3$x, x4=i4$x, x5=i5$x, x6=i6$x,
                   y1=i1$y, y2=i2$y, y3=i3$y, y4=i4$y, y5=i5$y, y6=i6$y, col=col)
  dt <- clip_polys(dt, y0, y1)
  iN <- rep(NA, nrow(dt))

  out.x <- as.numeric(t(matrix(c(dt$x1, dt$x2, dt$x3, dt$x4, dt$x5, dt$x6, iN), nrow=nrow(dt))))
  out.y <- as.numeric(t(matrix(c(dt$y1, dt$y2, dt$y3, dt$y4, dt$y5, dt$y6, iN), nrow=nrow(dt))))

  return(list(x=out.x, y=out.y, col=dt$col))
}

#' @name triangle
#' @title triangle
#' @description
#' internal draw.triangle function
#' @author Marcin Imielinski
#' @keywords internal
triangle <- function(x1, x2, y, y0, y1, col=NULL) {

  if (length(y) == 1)
    y = rep(y, length(x1))
  else if (length(y) != length(x1))
    warning('triangle: expecnting length(y) == length(x1) or length(y) == 1')
  if (length(x1) != length(x2))
    warning('triangle: expecting length(x1) == length(x2)')

  i1 = i2 = i3 = list(x=x1, y=y)
  i4 = i5 = .geti(x1, x2)
  i6 = list(x=x2, y=y)

  ## dummy
  if (is.null(col))
    col <- i1$x

  dt <- data.table::data.table(x1=i1$x, x2=i2$x, x3=i3$x, x4=i4$x, x5=i5$x, x6=i6$x,
                   y1=i1$y, y2=i2$y, y3=i3$y, y4=i4$y, y5=i5$y, y6=i6$y, col=col)
  dt <- clip_polys(dt, y0, y1)
  iN <- rep(NA, nrow(dt))

  out.x <- as.numeric(t(matrix(c(dt$x1, dt$x2, dt$x3, dt$x4, dt$x5, dt$x6, iN), nrow=nrow(dt))))
  out.y <- as.numeric(t(matrix(c(dt$y1, dt$y2, dt$y3, dt$y4, dt$y5, dt$y6, iN), nrow=nrow(dt))))

  return(list(x=out.x, y=out.y, col=col))
}


# @description
# internal draw.triangle function
# @keywords internal
.geti <- function(x1, x2) {
  if (length(x1) != length(x2))
    stop('x1 and x2 need to be same length')
  dt = data.frame(x1, x2)
  x = rowMeans(dt)
  dt$x1 = (-1) * dt$x1
  y = abs(rowMeans(dt)) # / slope
  list(x=x, y=y)
}

#' @keywords internal
.getpoints <- function (b1, b2) {

  i1 = .geti(b1[1], b2[1])
  i2 = .geti(b1[1], b2[2])
  i3 = .geti(b1[2], b2[2])
  i4 = .geti(b1[2], b2[1])

  return(list(x=c(i1[1], i2[1], i3[1], i4[1]), y=c(i1[2], i2[2], i3[2], i4[2])))
}

# internal draw.triangle function
# @keywords internal
.all.xpairs <- function(start, end) {
  if (length(start) != length(end)) {
    warning('.all.diamonds: lengths of input need to be same. Reducing')
    start <- start[min(length(start), length(end))]
    end <- end[min(length(start), length(end))]
  }

  if (length(start) == 1) {
    return(matrix())
  }

  num.boxes = sum(seq(1, length(start)-1))

  sr = seq(1, length(start)-1)
  rsr = rev(sr)
  ir = rep(sr, rsr)
  jr = unlist(lapply(seq(2,length(start)), function(x) seq(x,length(start))))

  bgx = matrix(c(start[ir], end[ir], start[jr], end[jr], ir, jr), nrow=num.boxes, ncol=6)

  return(bgx)

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
readData = function(..., environment = NULL)
{
  my.env  = new.env()
  data(..., envir = my.env);
  return(as.list(my.env))
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

gr.stripstrand = function(gr)
{
  GenomicRanges::strand(gr) = "*"
  return(gr)
}

format_windows <- function(windows, .Object) {
    if (is(windows, 'character'))
        windows = unlist(parse.grl(windows, seqlengths(seqinfo(.Object))))

    if (is(windows, 'Seqinfo'))
        windows = si2gr(windows)

    if (is(windows, 'GRangesList'))
        windows = unlist(windows)

    windows = windows[width(windows)>0]  ## otherwise will get non-matching below
    
    if (is.null(windows$col))
        windows$col = 'gray98'

#    if (is.list(windows))
#        windows = do.call('GRangesList', windows)

    ## collapse and match metadata back to original
    tmp = reduce(gr.stripstrand(windows))
    ix = gr.match(tmp, windows)
    values(tmp) = values(windows)[ix, ]    
    windows = tmp
##    if (!inherits(windows, 'GRangesList')) ## GRangesList windows deprecated
##        windows = GenomicRanges::GRangesList(windows)

    if (sum(as.numeric(width(grl.unlist(windows))))==0)
        {

            if (length(seqinfo(.Object))) {
                windows = si2gr(seqinfo(.Object))
            } else {
                warning("no windows provided and no seqinfo. Drawing blank plot")
                return(GRanges())
            }
        }

  return (grl.unlist(windows))
}

prep_defaults_for_plotting <- function(.Object) {
  # if (is.null(.Object@formatting$triangle))
  #   .Object@formatting$triangle = NA
  #
  # if (any(ix <- is.na(.Object@formatting$triangle)))
  #   .Object@formatting$triangle[ix] = FALSE
  #
  if (is.null(formatting(.Object)$source.file.chrsub))
    formatting(.Object)$source.file.chrsub = TRUE

  # layout legends if colorfield or colormap is NA and legends have no xpos set
  # leg.ix = which(.Object@formatting$legend & (!is.na(.Object@formatting$gr.colorfield) | !sapply(.Object@colormap, is.null)))
  # if (is.null(.Object@formatting$legend.xpos))
  #   .Object@formatting$legend.xpos = NA
  #
  # if (is.null(.Object@formatting$legend.xjust))
  #   .Object@formatting$legend.xjust = NA
  #
  # if (any(is.na(.Object@formatting$legend.xpos[leg.ix])))
  # {
  #   .Object@formatting$legend.xpos[leg.ix] = seq(0.1, 1, length.out = length(leg.ix))
  #   .Object@formatting$legend.xjust[leg.ix] = ifelse(.Object@formatting$legend.xpos[leg.ix]<0.5, 0, 1)
  # }

  # layout y coordinates
  sumh = sum(formatting(.Object)$height + formatting(.Object)$ygap)
  formatting(.Object)$height  = formatting(.Object)$height/sumh
  formatting(.Object)$ygap  = formatting(.Object)$ygap/sumh
  #            formatting(.Object)$ywid  = formatting(.Object)$ywid/sumh

  if (is.null(formatting(.Object)$max.ranges))
    formatting(.Object)$max.ranges = NA

  return(.Object)
}

extract_data_from_tmp_dat <- function(.Object, j, this.windows) {

  tmp.dat = dat(.Object)[[j]]

  # if (inherits(tmp.dat, 'RleList'))
  # {
  #   tmp.score = rle.query(tmp.dat, this.windows)
  #   tmp.dat = gUtils::gr.dice(this.windows)
  #   tmp.dat$score = as.numeric(tmp.score)
  #
  #   if (is.na(formatting(.Object)$y0[j]))
  #     formatting(.Object)$y0[j] = min(tmp.dat$score, na.rm = T)
  #
  #   if (is.na(formatting(.Object)$y1[j]))
  #     formatting(.Object)$y1[j] = max(tmp.dat$score, na.rm = T)
  #
  #   formatting(.Object)$y.field[j] = 'score'
  #   pre.filtered = TRUE
  # }
  ## else if
  if (is.character(tmp.dat))
  {
    if (file.exists(tmp.dat))
    {
      bed.style = FALSE
      if (grepl('(\\.bw)|(\\.bigwig)', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::BigWigFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.wig', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::WIGFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.bed', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::BEDFile(normalizePath(tmp.dat))
        #                            si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
        bed.style = TRUE
      }
      else if (grepl('\\.gff', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::GFFFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.2bit', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::TwoBitFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.bedgraph', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::BEDGraphFile(normalizePath(tmp.dat))
        #                           si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
        bed.style = TRUE
      }

      if (!bed.style) ## bed.style file objects do not have seqinfo option
      {
        si = tryCatch(seqinfo(f), error = function(x) '')
        pre.filtered = T

        if (!is.character(si))
        {
          if (formatting(.Object)[j, 'source.file.chrsub'])
            sel = gr.fix(gr.chr(this.windows), si, drop = TRUE)

          tmp.dat = rtracklayer::import(f, selection = sel, asRangedData = FALSE)

          if (formatting(.Object)[j, 'source.file.chrsub'])
            tmp.dat = gUtils::gr.sub(tmp.dat, 'chr', '')

          if (is.na(formatting(.Object)$y.field[j]))
            formatting(.Object)$y.field[j] = 'score'

        }
      }
      else
      {
        tmp.dat = rtracklayer::import(f, asRangedData = FALSE)

        if (formatting(.Object)[j, 'source.file.chrsub'])
          tmp.dat = gUtils::gr.sub(tmp.dat, 'chr', '')

        if (is.na(formatting(.Object)$y.field[j]))
          formatting(.Object)$y.field[j] = 'score'
      }
    }
  }
  else if (is(tmp.dat, 'ffTrack'))
  {
    tmp.dat = tmp.dat[this.windows, gr = TRUE]
    formatting(.Object)$y.field[j] = 'score'
    pre.filtered = T
  }
  else
  {
    if (formatting(.Object)[j, 'source.file.chrsub'])
      tmp.dat = gUtils::gr.sub(tmp.dat, 'chr', '')
  }


  if (is.character(tmp.dat)) ## file was not found
  {
    warnings('Track bigwig file not found')
    tmp.dat = GRanges()
  }

  return(list(o=.Object, t=tmp.dat))
}

enforce_max_ranges <- function(.Object, pre.filtered, j, tmp.dat, this.windows) {

  ## enforce max.ranges (everything should be a GRanges at this point)
  if (!pre.filtered & nrow(edgs(.Object)[[j]])==0 & length(.Object@vars[[j]])==0) ## assume any track associated with edges is pre-filtered
  {
    if (length(tmp.dat)>formatting(.Object)$max.ranges[j])
      if (inherits(tmp.dat, 'GRangesList'))
      {
        vals = values(tmp.dat)
        nm = names(tmp.dat)
        tmp2 = grl.unlist(tmp.dat)
        tmp2 = tmp2[gUtils::gr.in(tmp2, this.windows)]
        ##tmp2 = tmp2[tmp2 %over% this.windows]
        tmp.dat = GenomicRanges::split(tmp2, tmp2$grl.ix)
        values(tmp.dat) = vals[as.numeric(names(tmp.dat)), ]
        names(tmp.dat) = nm[as.numeric(names(tmp.dat))]
      }
    else if (length(tmp.dat)>formatting(.Object)$max.ranges[j])
    {
      tmp.dat = tmp.dat[gUtils::gr.in(tmp.dat, this.windows)]
      ##tmp.dat = tmp.dat[tmp.dat %over% this.windows]
      pre.filtered = TRUE
    }
  }

  if (length(tmp.dat)>formatting(.Object)$max.ranges[j] &
      nrow(edgs(.Object)[[j]])==0 & !.Object@formatting$triangle[j]) ## don't sample if there are edges or triangle
  {
    tmp.dat = sample(tmp.dat, ceiling(formatting(.Object)$max.ranges[j]))
  }

  return(list(p=pre.filtered, t=tmp.dat))

}

smooth_yfield <- function(.Object, j, tmp.dat) {

  tmp = S4Vectors::runmean(GenomicRanges::coverage(tmp.dat, weight = values(tmp.dat)[, formatting(.Object)$y.field[j]]),
                           k = floor(formatting(.Object)$smooth[j]/2)*2+1, endrule = 'constant', na.rm = TRUE)

  if (!is.na(formatting(.Object)$round[j]))
    tmp = round(tmp, formatting(.Object)$round[j])

  tmp = as(tmp, 'GRanges')
  tmp = tmp[gUtils::gr.in(tmp, tmp.dat)]
  ##tmp = tmp[tmp %over% tmp.dat]
  tmp.val = tmp$score
  values(tmp) = values(tmp.dat)[gr.match(tmp, tmp.dat), , drop = F]
  values(tmp)[, formatting(.Object)$y.field[j]] = tmp.val
  tmp.dat = tmp

}

format_yfield_limits <- function(.Object, j, tmp.dat, pre.filtered, this.windows) {

  strand(this.windows) <- rep("*", length(this.windows))

  if (!(formatting(.Object)$y.field[j] %in% names(values(tmp.dat))))
    stop('y.field missing from input granges')

  y0.global = min(values(tmp.dat)[, formatting(.Object)$y.field[j]], na.rm = TRUE)
  y1.global = max(values(tmp.dat)[, formatting(.Object)$y.field[j]], na.rm = TRUE)

  if (!pre.filtered)
    tmp.dat.r <- tmp.dat[gUtils::gr.in(tmp.dat, this.windows)]
    ##tmp.dat.r <- tmp.dat[tmp.dat %over% this.windows]
  else
    tmp.dat.r = tmp.dat

  val = values(tmp.dat.r)[, formatting(.Object)$y.field[j]]
  #                          r = range(val[!is.infinite(val)], na.rm = TRUE)

  p.quantile = 0.01
  if (!is.null(formatting(.Object)$y.quantile[j]) && !is.na(formatting(.Object)$y.quantile[j]))
    p.quantile = pmin(pmax(pmin(formatting(.Object)$y.quantile[j], 1-formatting(.Object)$y.quantile[j]), 0), 1)

  r = quantile(val[!is.infinite(val)], probs = c(p.quantile, 1-p.quantile), na.rm = TRUE)

  if (is.na(diff(r)))
    r = c(0, 0)

  if (is.na(formatting(.Object)$y0[j])) ## adjust y limits here if not specified
  {
    formatting(.Object)$y0[j] =  pmax(y0.global, r[1] - diff(r)*0.05)

    if (diff(r) == 0 & !is.na(formatting(.Object)$y1)[j])
      formatting(.Object)$y0[j] =
        r[1] - diff(c(r[1], formatting(.Object)$y1[j]))*0.05

  }

  if (is.na(formatting(.Object)$y1[j])) ## adjust y limits here if not specified
  {
    formatting(.Object)$y1[j] = pmin(y1.global, r[2] + diff(r)*0.05)
    if (diff(r) == 0 & !is.na(formatting(.Object)$y0[j]))
      formatting(.Object)$y1[j] =
        r[2] + diff(c(formatting(.Object)$y0[j], r[1]))*0.05
  }

  if (!is.null(formatting(.Object)$y0.bar[j]) && !is.na(formatting(.Object)$y0.bar[j]))
    formatting(.Object)$y0[j] = min(y0.global, formatting(.Object)$y0.bar[j])

  return(.Object)
}

draw_x_ticks <- function(xaxis.interval, windows, mapped, winlim, xlim, ylim, xaxis.pos, xaxis.suffix, xaxis.unit, xaxis.cex.tick, xaxis.ticklen, xaxis.round) {
  wid <- sum(as.numeric(width(windows)))
  xaxis.nticks = 20 ## default number of ticks
  if (!is.na(xaxis.interval) && xaxis.interval == 'auto') {
    if (wid > 100)
      xaxis.interval = 10^(ceiling(log10(wid/xaxis.nticks)))
    else
      xaxis.interval <- max(floor(wid/xaxis.nticks), 1)
  }

  # don't let xaxis.interval to be too small .. so tick drawing doesn't get too out of control
  xaxis.nticks = 1000 ## max number of ticks
  if (!is.na(xaxis.interval))
    xaxis.interval = max(xaxis.interval, 10^(ceiling(log10(sum(as.numeric(width(windows)))/xaxis.nticks))))

  seq.at.og = lapply(1:length(windows), function(x)
  {
    out = c(start(windows)[x], seq(ceiling(start(windows)[x]/xaxis.interval),
                                   floor(end(windows)[x]/xaxis.interval))*xaxis.interval, end(windows)[x])
    out[out>=start(windows)[x] & out<=end(windows)[x]]
  });

  winlim[2] <- winlim[2] + 1
  seq.at = lapply(1:length(seq.at.og), function(x) affine.map(seq.at.og[[x]]-start(windows)[x]+mapped$window.segs$start[x] + 0.5, ## added 0.5
                                                              xlim = winlim, ylim = xlim))
  seq.at.og = unlist(seq.at.og)
  seq.at = unlist(seq.at);
  dup.ix = duplicated(seq.at);
  seq.at.og = seq.at.og[!dup.ix]
  #seq.at.og <- seq.at.og[-length(seq.at.og)]## jeremiah
  seq.at = seq.at[!dup.ix]

  # then ticks
  seq.at <- seq.at[!is.na(seq.at)]
  tick.len = 0.01*diff(ylim)*xaxis.ticklen
  y0.tick = rep(xaxis.pos, length(seq.at))
  y1.tick = rep(xaxis.pos, length(seq.at)) - tick.len
  segments(seq.at[!is.na(seq.at)], y0.tick, seq.at, y1.tick)

  # then (tick) text
  if (xaxis.unit == 1)
    tick.text = prettyNum(paste(seq.at.og, xaxis.suffix), big.mark = ',')
  else
    tick.text = paste(round(seq.at.og/xaxis.unit, xaxis.round), xaxis.suffix)

  ##if (xaxis.nticks > 0)
  if (!is.na(xaxis.interval) && xaxis.interval > 0)
    text(seq.at, y1.tick-tick.len, tick.text,
         cex = xaxis.cex.tick, srt = 90, adj = c(1, 0.5))
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
                          slen = GenomeInfoDb::seqlengths(gr.fix(x))
                      else if (is(x, 'GRangesList'))
                          slen = GenomeInfoDb::seqlengths(gr.fix(unlist(x)))
                  }
          }
      else if (is(x, 'ffTrack'))
      {
        if (is.null(slen))
          slen = GenomeInfoDb::seqlengths(x)
        else
          slen[seqlevels(x)] = pmax(slen[seqlevels(x)], GenomeInfoDb::seqlengths(x), na.rm = TRUE)
        formatting(.Object)[i, 'yaxis'] = TRUE
      }
      else if (inherits(x, 'RleList'))
      {
        .Object@data[[i]] = as(x, 'RleList')
        formatting(.Object)[i, 'yaxis'] = TRUE
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
          formatting(.Object)[i, 'yaxis'] = TRUE
          if (is.null(slen))
            stop('External file must be a valid and existing .bigwig file')
        }
        else if (grepl('\\.wig', x, ignore.case = TRUE))
        {
          f = rtracklayer::WIGFile(normalizePath(x))
          slen = tryCatch(GenomeInfoDb::seqlengths(f), error = function(x) NULL)
          formatting(.Object)[i, 'yaxis'] = TRUE
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
          {
            if (.Object@formatting$source.file.chrsub)
              .Object@data[[i]] = gUtils::gr.sub(tmp.out, 'chr', '')
            else
              .Object@data[[i]] = tmp.out
          }

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
          {
            if (.Object@formatting$source.file.chrsub)
              .Object@data[[i]] = gUtilts::gr.sub(tmp.out, 'chr', '')
            else
              .Object@data[[i]] = tmp.out
          }

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
          {
            if (.Object@formatting$source.file.chrsub)
              .Object@data[[i]] = gr.sub(tmp.out, 'chr', '')
            else
              .Object@data[[i]] = tmp.out
          }

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
        slen = data.table(len = slen, nm = names(slen))[, max(len), by = nm][, structure(len, names = nm)]
            
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
      .Object@data = lapply(.Object@data, gr.fix, seqinfo)
    }

  }
  return(.Object)
}

## convert input data into a list of length 'len' of type FUN
listify <- function(x, FUN, len = 1) {
  if (is.null(x))
    return(rep(list(FUN()), len))
  if (!is.list(x))
    return(list(x))
  return(x)
}


