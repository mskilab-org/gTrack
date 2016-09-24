#' @name gxaxis-class
#' @rdname gxaxis-class
#' @title S4 class for x-axis params of \code{gTrack}
#'
#' Class \code{gxaxis} is a simple container storing x-axis plotting information for gTrack
#'
#' @exportClass gxaxis
#' @importFrom methods setClass setGeneric setMethod setRefClass
setClass('gxaxis', representation(formatting = 'data.frame'))

setMethod('initialize', 'gxaxis', function(.Object, ...) ## only place NON formatting fields here. The rest passed with ...
{
  ## set the defaults
  .Object@formatting = as.data.frame(list(...), stringsAsFactors = FALSE)
  return(.Object)
})

#' @param xaxis.suffix vector or scalar numeric specifying the suffix that will be used in describin x axis coordinates (TODO: move to display) (formatting)
#' @param xaxis.unit vector or scalar numeric specifying the unit that will be used in describing x axis coordinates (TODO: move to display)  (formatting)
#' @param xaxis.round vector or scalar non-neg integer specifying number of decimals to round xaxis coordinate labels (formatting)
#' @param xaxis.nticks vector or scalar positive integer specifying how many xaxis ticks to optimally draw (formatting)
#' @param xaxis.label.angle vector or scalar numeric between 0 and 360 specifying angle with which to draw xaxis coordinate labels (formatting)
#' @param xaxis.newline vector or scalar logical specifying whether to draw a newline in the xaxis coordinate labels (formatting)
#' @param xaxis.width Logical scalar specifying whether to add window width to xaxis window labels [TRUE]
#' @param xaxis.cex.tick Scalar numeric specifying expansion factor for axis tick labels  (formatting)
#' @param xaxis.ticklen Scalar numeric specifying lengths for axis ticks  (formatting)
#' Construct a new \code{gxaxis}
#' @name gxaxis
#' @rdname gxaxis-class
#' @export
gxaxis = function(xaxis.prefix = "", xaxis.unit = 1,
                  xaxis.suffix = "", xaxis.round = 3,
                  xaxis.cex.label = 1, xaxis.newline = TRUE,
                  xaxis.chronly = FALSE, xaxis.width= TRUE,
                  xaxis.interval = 'auto', xaxis.label.angle = 0,
                  xaxis.ticklen = 1, xaxis.cex.tick = 1) {
  new('gxaxis', xaxis.prefix = xaxis.prefix, xaxis.unit = xaxis.unit,
      xaxis.suffix = xaxis.suffix, xaxis.round = xaxis.round,
      xaxis.cex.label = xaxis.cex.label, xaxis.newline = xaxis.newline,
      xaxis.chronly = xaxis.chronly, xaxis.width= xaxis.width,
      xaxis.interval = xaxis.interval, xaxis.label.angle = xaxis.label.angle,
      xaxis.ticklen = xaxis.ticklen, xaxis.cex.tick = xaxis.cex.tick)
}

#' @name show
#' @title show
#' @description Display a \code{gxaxis} object
#' @docType methods
#' @param object \code{gxaxis} to display
#' @rdname gxaxis-show-methods
#' @aliases show,gxaxis-method
#' @export
#' @author Marcin Imielinski
setMethod('show', 'gxaxis', function(object)
{
  cat(sprintf('gxaxis object with params:\n'))
  print(object@formatting)
})

#' @name $
#' @title $
#' @description
#'
#' Accessing columns of gxaxis formatting data.frame
#'
#' @param x \code{gxaxis} object
#' @param name Name of the \code{formatting} field to view
#' @docType methods
#' @rdname gxaxis-cash-methods
#' @aliases $,gxaxis-method
#' @export
#' @author Marcin Imielinski
setMethod('$', 'gxaxis', function(x, name)
{
  return(x@formatting[, name])
})

#' @keywords internal
get.end.text <- function(x, windows) {
  end.text = prettyNum(ifelse(rep(x$xaxis.unit == 1, length(windows)), end(windows),
                              round(end(windows)/x$xaxis.unit, x$xaxis.round)), big.mark = ',')
  return (end.text)
}

#' @keywords internal
get.begin.text <- function(x, windows) {

  begin.text = prettyNum(pmax(floor(1/x$xaxis.unit),
                              ifelse(rep(x$xaxis.unit == 1, length(windows)),
                                     start(windows), round(start(windows)/x$xaxis.unit,
                                     x$xaxis.round))), big.mark = ',')
  return (begin.text)
}

#' @keywords internal
get.width.text <- function(x, windows) {
  if (x$xaxis.width)
    return('')

  width.text = ''
    if (!is.null(x$xaxis.suffix))
      width.text = paste('(', paste(prettyNum(ifelse(rep(x$xaxis.unit == 1, length(windows)),
                                                     width(windows), round(width(windows)/x$xaxis.unit, 2)),
                                              big.mark = ','), x$xaxis.suffix),  ')', sep = '')
    else
      width.text = paste('(', prettyNum(ifelse(rep(x$xaxis.unit == 1, length(windows)),
                                               width(windows), round(width(windows)/x$xaxis.unit, 2)),
                                        big.mark = ','),  ')', sep = '')
    return (width.text)
}
