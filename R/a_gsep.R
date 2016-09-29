#' @name gsep-class
#' @rdname gsep-class
#' @title S4 class for legend params of \code{gTrack}
#'
#' Class \code{gsep} is a simple container storing legend plotting information for gTrack
#'
#' @exportClass gsep
#' @importFrom methods setClass setGeneric setMethod setRefClass
setClass('gsep', representation(formatting = 'list'))

setMethod('initialize', 'gsep', function(.Object, ...) {
  ## set the defaults
  .Object@formatting = list(...)
  return(.Object)
})


#' Construct a new \code{gsep}
#'
#' @param lty Scalar integer specifying line style for window separators [2] (dashed line)
#' @param lwd Scalar numeric specifying line thickness for window separators [1]
#' @param bg.col Color (supplied by name or hex code) for the background of the windows [gray95]
#' @param plot Logical allowing separators to be turned off [FALSE]
#' what is the maximum number of items to include in legend [Inf]
#' @name gsep
#' @rdname gsep-class
#' @export
gsep = function(plot = TRUE, lty = 1,
                lwd = 1, bg.col = "gray95",
                ...) {
  new('gsep', plot = plot, lty = lty, lwd = lwd, bg.col = bg.col, ...)
}

#' @name show
#' @title show
#' @description Display a \code{gsep} object
#' @docType methods
#' @param object \code{gsep} to display
#' @rdname gsep-show-methods
#' @aliases show,gsep-method
#' @export
setMethod('show', 'gsep', function(object) {
  cat(sprintf('gsep object with params:\n'))
  print(object@formatting)
})

#' @name $<-
#' @title $<-
#' @description
#'
#' Setting formats of gsep object - ie modifying the formatting(gt) data frame after
#' an object has already been instantiated
#'
#' @param x \code{gsep} object to alter \code{formatting} field of
#' @param name \code{formatting} field to alter
#' @param value New value
#' @docType methods
#' @rdname gsep-cash-set-methods
#' @aliases $<-,gsep-method
#' @export
setMethod('$<-', 'gsep', function(x, name, value) {
  x@formatting[[name]] = value
  return(x)
})

#' @name $
#' @title $
#' @description
#'
#' Accessing columns of gsep formatting data.frame
#'
#' @param x \code{gsep} object
#' @param name Name of the \code{formatting} field to view
#' @docType methods
#' @rdname gsep-cash-methods
#' @aliases $,gsep-method
#' @export
setMethod('$', 'gsep', function(x, name) {
  return(x@formatting[[name]])
})
