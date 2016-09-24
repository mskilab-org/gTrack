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

        rect(x$pos1, ifelse(y0.bar<x$y, rep(y0.bar, nrow(x)), x$y), x$pos2, ifelse(y0.bar<x$y, x$y, rep(y0.bar, nrow(x))),
             col = x$col, border = x$col, lwd = x$lwd.border/2)
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
