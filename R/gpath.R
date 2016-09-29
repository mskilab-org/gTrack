#' @name gpath-class
#' @rdname gpath-class
#' @title S4 class for legend params of \code{gTrack}
#'
#' Class \code{gpath} is a simple container storing legend plotting information for gTrack
#'
#' @exportClass gpath
#' @importFrom methods setClass setGeneric setMethod setRefClass
setClass('gpath', representation(formatting = 'list'))

setMethod('initialize', 'gpath', function(.Object, ...)
{
  ## set the defaults
  .Object@formatting = list(...)
  return(.Object)
})

#' Construct a new \code{gpath}
#'
#' @param col scalar character specifying color of path (only applicable for tracks with draw.paths = TRUE)
#' @param col.arrow scalar character specifying color of arrow of path
#' @param arrow scalar numeric > 0 specifying expansion factor of arrow of path
#' @param stack.y.gap scalar numeric > 0 specifying y stack gap of paths
#' @param stack.x.gap scalar numeric > 0 specifying x stack gap for paths
#' @param cex.v scalar numeric > 0 specifying vertical bulge of sline in paths
#' @param cex.h scalar numeric > 0 specifying horizontal bulge of spline in paths
#' @name gpath
#' @rdname gpath-class
#' @export
gpath = function(cex.v = 1, cex.h = 1, cex.arrow = 1, col = 'black',
                 col.arrow = 'black', stack.y.gap = 1, stack.x.gap = 0,
                 ...) {
  new('gpath', cex.v = cex.v, cex.h = cex.h, cex.arrow = cex.arrow,
      col = col, col.arrow = col.arrow,
      stack.y.gap = stack.y.gap, stack.x.gap = stack.x.gap, ...)
}



#' @name show
#' @title show
#' @description Display a \code{gpath} object
#' @docType methods
#' @param object \code{gpath} to display
#' @rdname gpath-show-methods
#' @aliases show,gpath-method
#' @export
#' @author Marcin Imielinski
setMethod('show', 'gpath', function(object)
{
  cat(sprintf('gpath object with params:\n'))
  print(object@formatting)
})

#' @name $<-
#' @title $<-
#' @description
#'
#' Setting formats of gpath object - ie modifying the formatting(gt) data frame after
#' an object has already been instantiated
#'
#' @param x \code{gpath} object to alter \code{formatting} field of
#' @param name \code{formatting} field to alter
#' @param value New value
#' @docType methods
#' @rdname gpath-cash-set-methods
#' @aliases $<-,gpath-method
#' @export
setMethod('$<-', 'gpath', function(x, name, value)
{
  x@formatting[[name]] = value
  return(x)
})

#' @name $
#' @title $
#' @description
#'
#' Accessing columns of gpath formatting data.frame
#'
#' @param x \code{gpath} object
#' @param name Name of the \code{formatting} field to view
#' @docType methods
#' @rdname gpath-cash-methods
#' @aliases $,gpath-method
#' @export
#' @author Marcin Imielinski
setMethod('$', 'gpath', function(x, name)
{
  return(x@formatting[[name]])
})

draw.connector.paths <- function(grl.segs, path) {

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

if (!is.null(connector.args)) {
  path.h = path$cex.h * rep(diff(par('usr')[1:2])/50, nrow(connector.args))
  if (any(connector.args$type == 'S'))
    path.h[connector.args$type == 'S'] = 2*path.h[connector.args$type == 'S']

  path.v = rep(path$cex.v, nrow(connector.args))
  path.v[is.na(connector.args$ix0) | is.na(connector.args$ix1)] = path$cex.v*2*grl.segs$ywid[connector.args$ix0[is.na(connector.args$ix0) | is.na(connector.args$ix1)]]

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
  path.v[connector.args$cyclic] = path$cex.v*2*grl.segs$ywid[connector.args$ix1[connector.args$cyclic]]
  path.h[connector.args$cyclic] = path$cex.h * diff(par('usr')[1:2])/20

  ## workaround current lines() limitation in connectors
  if (any(lty3 <- connector.args$lty == 3))
    connectors(connector.args$x0[lty3], connector.args$y0[lty3], connector.args$strand0[lty3],
               connector.args$x1[lty3], connector.args$y1[lty3], ifelse(connector.args$strand1[lty3] == '+', '-', '+'),
               type = connector.args$type[lty3],
               h = path.h[lty3], v = path.v[lty3],
               lty = connector.args$lty[lty3], col = path$col, col.arrow = path$col.arrow,
               cex.arrow = grl.segs$ywid[1]*path$cex.arrow, f.arrow = T)

  if (any(lty1 <- connector.args$lty == 1))
    connectors(connector.args$x0[lty1], connector.args$y0[lty1], connector.args$strand0[lty1],
               connector.args$x1[lty1], connector.args$y1[lty1], ifelse(connector.args$strand1[lty1] == '+', '-', '+'),
               type = connector.args$type[lty1],
               h = path.h[lty1], v = path.v[lty1],
               lty = connector.args$lty[lty1], col = path$col, col.arrow = path$col.arrow,
               cex.arrow = grl.segs$ywid[1]*path$cex.arrow, f.arrow = T)


  #text(connector.args$x1, connector.args$y1, paste(connector.args$ix0, ' (', grl.segs$group[connector.args$ix0], '), ', connector.args$ix1, ' (', grl.segs$group[connector.args$ix1], '), ', connector.args$type, ' ', connector.args$sign, sep = ''))
  #            text(connector.args$x1, connector.args$y1-path.v, paste(grl.segs$group[connector.args$ix0], connector.args$type, ' ', connector.args$sign, sep = ''), )
}

}
