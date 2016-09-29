#' @name glegend-class
#' @rdname glegend-class
#' @title S4 class for legend params of \code{gTrack}
#'
#' Class \code{glegend} is a simple container storing legend plotting information for gTrack
#'
#' @exportClass glegend
#' @importFrom methods setClass setGeneric setMethod setRefClass
setClass('glegend', representation(formatting = 'list'))

setMethod('initialize', 'glegend', function(.Object, ...)
{
  ## set the defaults
  .Object@formatting = list(...)
  return(.Object)
})


#' Construct a new \code{glegend}
#'
#' @param maxitems Scalar positive integer specifying
#' what is the maximum number of items to include in legend [Inf]
#' @name glegend
#' @rdname glegend-class
#' @export
glegend = function(xpos = 0, ypos = 1,
                   plot = TRUE, title = "", maxitems = Inf,
                   x= xpos, y=ypos, legend = NA,
                   ...) {
  new('glegend', xpos=xpos, ypos=ypos, plot=plot, title=title,
      x=x, y=y, legend = legend, ...)
}

## provide a legend from a color map if one is not explitily
## provided.
#' @keywords internal
calc.legs <- function(leg.params, xlim, gr.colormap, plot.type) {

  leg.params$x = leg.params$xpos * diff(xlim) + xlim[1]
  leg.params$y = leg.params$ypos * diff(par('usr')[3:4]) + par('usr')[3]
  if (is.na(leg.params$legend)) {
    if (is.na(gr.colormap[1]))
      leg.params$legend = ""
    else
      leg.params$legend = names(gr.colormap)
  }

  if (plot.type == 'scatter') {
    leg.params$col = gr.colormap
    leg.params$pch = 16
  }
  else
    leg.params$fill = gr.colormap
  leg.params$border = gr.colormap
  leg.params$xpos = leg.params$ypos = NULL
  return (leg.params)
}

#' @name show
#' @title show
#' @description Display a \code{glegend} object
#' @docType methods
#' @param object \code{glegend} to display
#' @rdname glegend-show-methods
#' @aliases show,glegend-method
#' @export
#' @author Marcin Imielinski
setMethod('show', 'glegend', function(object)
{
  cat(sprintf('glegend object with params:\n'))
  print(object@formatting)
})

#' @name $<-
#' @title $<-
#' @description
#'
#' Setting formats of glegend object - ie modifying the formatting(gt) data frame after
#' an object has already been instantiated
#'
#' @param x \code{glegend} object to alter \code{formatting} field of
#' @param name \code{formatting} field to alter
#' @param value New value
#' @docType methods
#' @rdname glegend-cash-set-methods
#' @aliases $<-,glegend-method
#' @export
setMethod('$<-', 'glegend', function(x, name, value)
{
  x@formatting[[name]] = value
  return(x)
})

#' @name $
#' @title $
#' @description
#'
#' Accessing columns of glegend formatting data.frame
#'
#' @param x \code{glegend} object
#' @param name Name of the \code{formatting} field to view
#' @docType methods
#' @rdname glegend-cash-methods
#' @aliases $,glegend-method
#' @export
#' @author Marcin Imielinski
setMethod('$', 'glegend', function(x, name)
{
  return(x@formatting[[name]])
})

get.legend.params <- function(x) {
  return (x@formatting[!names(x@formatting) %in% "maxitems"])
}
