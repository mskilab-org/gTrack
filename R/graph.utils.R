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
