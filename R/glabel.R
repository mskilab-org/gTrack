#' @name glabel-class
#' @rdname glabel-class
#' @title S4 class for legend params of \code{gTrack}
#'
#' Class \code{glabel} is a simple container storing legend plotting information for gTrack
#'
#' @exportClass glabel
#' @importFrom methods setClass setGeneric setMethod setRefClass
setClass('glabel', representation(formatting = 'list'))

setMethod('initialize', 'glabel', function(.Object, ...) {
  ## set the defaults
  .Object@formatting = list(...)
  return(.Object)
})

#' Construct a new \code{glabel}
#'
#' @name glabel
#' @rdname glabel-class
#' @export
glabel <- function(gr.cex = 1, grl.cex = 1, grl.labelfield = NA,
                  gr.labelfield = NA, gr.angle = 0, grl.angle = 0,
                  hadj = 1, vadj = 0.5,
                ...) {
  new('glabel', gr.cex = gr.cex, grl.cex = grl.cex, grl.labelfield = grl.labelfield,
      gr.labelfield = gr.labelfield, gr.angle = gr.angle, grl.angle = grl.angle,
      hadj = hadj, vadj = vadj, ...)
}

#' @name show
#' @title show
#' @description Display a \code{glabel} object
#' @docType methods
#' @param object \code{glabel} to display
#' @rdname glabel-show-methods
#' @aliases show,glabel-method
#' @export
setMethod('show', 'glabel', function(object) {
  cat(sprintf('glabel object with params:\n'))
  print(object@formatting)
})

#' @name $<-
#' @title $<-
#' @description
#'
#' Setting formats of glabel object - ie modifying the formatting(gt) data frame after
#' an object has already been instantiated
#'
#' @param x \code{glabel} object to alter \code{formatting} field of
#' @param name \code{formatting} field to alter
#' @param value New value
#' @docType methods
#' @rdname glabel-cash-set-methods
#' @aliases $<-,glabel-method
#' @export
setMethod('$<-', 'glabel', function(x, name, value) {
  x@formatting[[name]] = value
  return(x)
})

#' @name $
#' @title $
#' @description
#'
#' Accessing columns of glabel formatting data.frame
#'
#' @param x \code{glabel} object
#' @param name Name of the \code{formatting} field to view
#' @docType methods
#' @rdname glabel-cash-methods
#' @aliases $,glabel-method
#' @export
setMethod('$', 'glabel', function(x, name) {
  return(x@formatting[[name]])
})
