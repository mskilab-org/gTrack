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
