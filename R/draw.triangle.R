

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
                          separator,
                          sigma = NA, ## if not NA then will blur with a Gaussian filter using a sigma value of this many base pairs
                          col.min='white',
                          col.max='red',
                          gr.colormap = NA,
                          cmap.min, cmap.max,
                          leg.params,
                          min.gapwidth = 1,
                          islog = FALSE) {

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
  winlim = range(c(window.segs$start, window.segs$end))
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
  grl.segs$pos1 = affine.map(grl.segs$pos1, winlim, ylim = xlim)
  grl.segs$pos2 = affine.map(pmin(grl.segs$pos2+1, winlim[2]), winlim, ylim = xlim) ## 160925 added +1
  window.segs$start = affine.map(window.segs$start, winlim, ylim = xlim)
  window.segs$end = affine.map(window.segs$end, winlim, ylim = xlim)

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

  ## empty, so just draw bounding box and leave
  if (nrow(grl.segs) == 0) {
    draw_bounding_triangle(window.segs, y0, y1, dlim, ylim.subplot, separator)
    return(window.segs)
  }

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

  ## now plot the box around the data
  draw_bounding_triangle(window.segs, y0, y1, dlim, ylim.subplot, separator)

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

draw_bounding_triangle <- function(window.segs, y0, y1, dlim, ylim.subplot, separator) {

  ## draw the background BOXES
  if (nrow(window.segs) > 1) {
    bgx = .all.xpairs(window.segs$start, window.segs$end)
    out <- diamond(bgx[,1], bgx[,2], bgx[,3], bgx[,4], y0, y1)

    ## affine map to local coorinates
    out$y[!is.na(out$y)] <- affine.map(out$y[!is.na(out$y)], xlim=dlim, ylim=ylim.subplot)
    if (separator$lwd > 0)
      polygon(out$x, out$y, col=NA, lwd=separator$lwd) ## should be NA col because then it is see-through
    else
      polygon(out$x, out$y, col=NA, border=NA)
  }

  ## draw the background triangles
  out = triangle(x1=window.segs$start, x2=window.segs$end, y=y0, y0=y0, y1=y1)
  out$y[!is.na(out$y)] <- affine.map(out$y[!is.na(out$y)], xlim=dlim, ylim=ylim.subplot)
  if (separator$lwd > 0)
    polygon(out$x, out$y, col=NA, lwd=separator$lwd)
  else
    polygon(out$x, out$y, col=NA, border=NA)
}
