#' @name gaxis-class
#' @rdname gaxis-class
#' @title S4 class for x-axis params of \code{gTrack}
#'
#' Class \code{gaxis} is a simple container storing x-axis plotting information for gTrack
#'
#' @exportClass gaxis
#' @importFrom methods setClass setGeneric setMethod setRefClass
setClass('gaxis', representation(formatting = 'list'))

setMethod('initialize', 'gaxis', function(.Object, ...) ## only place NON formatting fields here. The rest passed with ...
{
  ## set the defaults
  .Object@formatting = list(...)
  .Object@formatting = append(.Object@formatting,
                              list(x.pos=0, y.pos=0, x.label.pos = 0,
                                   y.label.pos = 0))
  return(.Object)
})

#' @param suffix vector or scalar numeric specifying the suffix that will be used in describin x axis coordinates (TODO: move to display) (formatting)
#' @param unit vector or scalar numeric specifying the unit that will be used in describing x axis coordinates (TODO: move to display)  (formatting)
#' @param round vector or scalar non-neg integer specifying number of decimals to round axis coordinate labels (formatting)
#' @param nticks vector or scalar positive integer specifying how many axis ticks to optimally draw (formatting)
#' @param label.angle vector or scalar numeric between 0 and 360 specifying angle with which to draw axis coordinate labels (formatting)
#' @param newline vector or scalar logical specifying whether to draw a newline in the axis coordinate labels (formatting)
#' @param width Logical scalar specifying whether to add window width to axis window labels [TRUE]
#' @param cex.tick Scalar numeric specifying expansion factor for axis tick labels  (formatting)
#' @param ticklen Scalar numeric specifying lengths for axis ticks  (formatting)
#' @param quantile if \code{min} or \code{max} is not specified then will draw between quantile and 1-quantile of data
#' @param cap whether to cap values at min, min (only relevant if y.field specified)
#' Construct a new \code{gaxis}
#' @name gaxis
#' @rdname gaxis-class
#' @export
gaxis = function( tick.adjust.x = 1.0001, tick.adjust.y = 0.50001,
                  label.adjust.x = 0.5, label.adjust.y = 0,
                  prefix = "", unit = 1,
                  suffix = "", round = 3,
                  cex.label = 1, newline = TRUE,
                  chronly = FALSE, width= TRUE,
                  interval = 'auto', label.angle = 0,
                  ticklen = 1, cex.tick = 1,
                  min=NA, max=NA, min.bar = NA, pretty=TRUE,
                  plot = TRUE, quantile = 0.01, cap = TRUE,
                  tick.angle=90, log=FALSE,
                  breaks = NA, labels = NA) {

  if (!is.na(breaks[1])) {
    if (!is.list(breaks)) ## if only one break vec give, convert to list
      breaks = list(breaks)
  }

  ## set a different default for tick angle of 0
  if (tick.adjust.x == 1.0001 && tick.adjust.y == 0.50001 && tick.angle == 0) {
    tick.adjust.x = 0.5
    tick.adjust.y = 1
  }

  new('gaxis',
      breaks = breaks, labels = labels,
      prefix = prefix, unit = unit,
      suffix = suffix, round = round,
      cex.label = cex.label, newline = newline,
      chronly = chronly, width= width,
      interval = interval, label.angle = label.angle,
      ticklen = ticklen, cex.tick = cex.tick, min=min, max=max,
      min.bar = min.bar, pretty=pretty, plot = plot,
      quantile = quantile, cap = cap, tick.angle=tick.angle,
      log = log,
      tick.adjust.x = tick.adjust.x,
      tick.adjust.y = tick.adjust.y,
      label.adjust.x = label.adjust.x,
      label.adjust.y = label.adjust.y)
}

#' @name show
#' @title show
#' @description Display a \code{gaxis} object
#' @docType methods
#' @param object \code{gaxis} to display
#' @rdname gaxis-show-methods
#' @aliases show,gaxis-method
#' @export
#' @author Marcin Imielinski
setMethod('show', 'gaxis', function(object)
{
  cat(sprintf('gaxis object with params:\n'))
  print(object@formatting)
})

#' @name $
#' @title $
#' @description
#'
#' Accessing columns of gaxis formatting data.frame
#'
#' @param x \code{gaxis} object
#' @param name Name of the \code{formatting} field to view
#' @docType methods
#' @rdname gaxis-cash-methods
#' @aliases $,gaxis-method
#' @export
#' @author Marcin Imielinski
setMethod('$', 'gaxis', function(x, name)
{
  return(x@formatting[[name]])
})

#' @name $<-
#' @title $<-
#' @description
#'
#' Setting formats of gaxis object - ie modifying the formatting(gt) data frame after
#' an object has already been instantiated
#'
#' @param x \code{gaxis} object to alter \code{formatting} field of
#' @param name \code{formatting} field to alter
#' @param value New value
#' @docType methods
#' @rdname gaxis-cash-set-methods
#' @aliases $<-,gaxis-method
#' @export
setMethod('$<-', 'gaxis', function(x, name, value)
{
  x@formatting[[name]] = value
  return(x)
})


#' @keywords internal
get.end.text <- function(x, windows) {
  end.text = prettyNum(ifelse(rep(x$unit == 1, length(windows)), end(windows),
                              round(end(windows)/x$unit, x$round)), big.mark = ',')
  return (end.text)
}

#' @keywords internal
get.begin.text <- function(x, windows) {

  begin.text = prettyNum(pmax(floor(1/x$unit),
                              ifelse(rep(x$unit == 1, length(windows)),
                                     start(windows), round(start(windows)/x$unit,
                                     x$round))), big.mark = ',')
  return (begin.text)
}

#' @keywords internal
get.width.text <- function(x, windows) {
  if (x$width)
    return('')

  width.text = ''
    if (!is.null(x$suffix))
      width.text = paste('(', paste(prettyNum(ifelse(rep(x$unit == 1, length(windows)),
                                                     width(windows), round(width(windows)/x$unit, 2)),
                                              big.mark = ','), x$suffix),  ')', sep = '')
    else
      width.text = paste('(', prettyNum(ifelse(rep(x$unit == 1, length(windows)),
                                               width(windows), round(width(windows)/x$unit, 2)),
                                        big.mark = ','),  ')', sep = '')
    return (width.text)
}

#' @keywords internal
draw_x_ticks <- function(x, windows, mapped, winlim, xlim, ylim) {

  wid <- sum(as.numeric(width(windows)))
  xaxis.nticks = 20 ## default number of ticks
  if (!is.na(x$interval) && x$interval == 'auto') {
    if (wid > 100)
      x$interval = 10^(ceiling(log10(wid/xaxis.nticks)))
    else
      x$interval <- max(floor(wid/xaxis.nticks), 1)
  }

  # don't let xaxis.interval to be too small .. so tick drawing doesn't get too out of control
  xaxis.nticks = 1000 ## max number of ticks
  if (!is.na(x$interval))
    x$interval = max(x$interval, 10^(ceiling(log10(sum(as.numeric(width(windows)))/xaxis.nticks))))

  OFFSET = ifelse(wid > 100, 0, 0.5) ## for small data, think of as categorical
  if (!is.na(x$breaks[1])) { ## break are explicitly given
    if (length(x$breaks) != length(windows))
      stop("axis breaks must be list of length as num windows")
    seq.at.og <- lapply(x$breaks, function(z) return(z + OFFSET))
  } else {
    seq.at.og = lapply(1:length(windows), function(z) {
      out = c(start(windows)[z],
              seq(ceiling(start(windows)[z]/x$interval),
                  floor(end(windows)[z]/x$interval))*x$interval, end(windows)[z])
      out[out>=start(windows)[z] & out<=end(windows)[z]]
      return (out + OFFSET)
    })
  }

  seq.at = lapply(1:length(seq.at.og), function(z)
    affine.map(seq.at.og[[z]]-start(windows)[z]+mapped$window.segs$start[z], ## add 0.5 for offset
               xlim = winlim, ylim = xlim))
  seq.at.og = unlist(seq.at.og)
  seq.at = unlist(seq.at)
  dup.ix = duplicated(seq.at)
  seq.at.og = seq.at.og[!dup.ix]
  seq.at = seq.at[!dup.ix]

  # then ticks
  seq.at <- seq.at[!is.na(seq.at)]
  tick.len = 0.01*diff(ylim)*x$ticklen
  y0.tick = rep(x$y.pos, length(seq.at))
  y1.tick = rep(x$y.pos, length(seq.at)) - tick.len
  segments(seq.at[!is.na(seq.at)], y0.tick, seq.at, y1.tick)

  # then (tick) text
  seq.at.og = seq.at.og - OFFSET ## set it back
  if (!is.na(x$breaks[1])) {
    if (!is.null(x$breaks[[1]])) ## are the names provided in breaks?
      tick.text = names(seq.at.og)
    else
      names = as.character(seq.at.og)
  } else if (x$unit == 1)
    tick.text = prettyNum(paste(seq.at.og, x$suffix), big.mark = ',')
  else
    tick.text = paste(round(seq.at.og/x$unit, x$round), x$suffix)

  ##if (xaxis.nticks > 0)
  if (!is.na(x$interval) && x$interval > 0)
    text(seq.at, y1.tick-tick.len, tick.text,
         cex = x$cex.tick, srt = x$tick.angle,
         adj = c(x$tick.adjust.x, x$tick.adjust.y))
}

#' @keywords internal
draw_x_label <- function(xaxis, window.segs, windows) {
  newline <- ifelse(xaxis$newline, '\n', '')

  if (!is.na(xaxis$labels[1])) { ## label names explicitly provided
    if (length(xaxis$labels) != nrow(window.segs))
      stop("axis names must be same length as windows")
    labs <- xaxis$labels
  } else {  ## label names automatically generated
    width.text <- get.width.text(xaxis, windows)
    begin.text <- get.begin.text(xaxis, windows)
    end.text   <- get.end.text(xaxis, windows)
    if (!xaxis$chronly)
      labs <- paste(xaxis$prefix, ' ',  seqnames(windows), ':',newline,
                    begin.text,'-', newline,
                    end.text, ' ', xaxis$suffix, newline, width.text, sep = '')
    else
      paste(xaxis$prefix, ' ',  seqnames(windows), sep = '')
  }

  text(rowMeans(window.segs[, c('start', 'end')]), rep(xaxis$y.pos.label, length(windows)),
       labs,
       cex = xaxis$cex.label, adj = c(xaxis$label.adjust.x, xaxis$label.adjust.y),
       srt=xaxis$label.angle)

}

