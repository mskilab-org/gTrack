#' @keywords internal
draw.links <- function(links, win.u, window.segs) {

  if (length(links)==0)
    return();

  # first map rearrangements to various windows>
  win.u$grl.ix = 1  ##holdover from grangeslist windows
  ##win.u = gr.stripstrand(grl.unlist(windows))
  window.segs.u = do.call(rbind, window.segs)
  window.segs.u$width = window.segs.u$end - window.segs.u$start + 1
  window.segs.xlim = do.call('rbind', lapply(window.segs, function(x) data.frame(start = min(x$start), end = max(x$end))))

  links.u = grl.unlist(links)

  if (any(table(links.u$grl.ix)!=2))
    stop('Links should be GRangesList of range pairs.')

  links.p = grl.pivot(links)

  ## flip strands to conform to connectors convention (- connection to left side and + is connection to right side)
  ## with ra specification convention (- refers to segment to left of breakpoint and + refers to segment to right)
  GenomicRanges::strand(links.p[[1]]) = c("-" = "+", "+" = "-")[as.character(GenomicRanges::strand(links.p[[1]]))]
  GenomicRanges::strand(links.p[[2]]) = c("-" = "+", "+" = "-")[as.character(GenomicRanges::strand(links.p[[2]]))]

  ## find overlaps with windows and calculate their window specific coordinates
  l1 = gr.findoverlaps(links.p[[1]], win.u)
  values(l1) = cbind(as.data.frame(values(l1)), as.data.frame(values(links)[l1$query.id, , drop = FALSE]))
  GenomicRanges::strand(l1) = GenomicRanges::strand(links.p[[1]])[l1$query.id]
  l1$stack.id = win.u$grl.ix[l1$subject.id]
  l1$y.pos = ylim.stacks$xaxis.pos[l1$stack.id]
  l1$x.pos = mapply(function(x,y,z,a) (y-z)*a + x, x = window.segs.u[l1$subject.id,]$start, y = start(l1),
                    z = start(win.u[l1$subject.id]), a = window.segs.u$width[l1$subject.id] / width(win.u)[l1$subject.id])

  l2 = gr.findoverlaps(links.p[[2]], win.u)
  values(l2) = cbind(as.data.frame(values(l2)), as.data.frame(values(links)[l2$query.id, , drop = FALSE]))
  GenomicRanges::strand(l2) = GenomicRanges::strand(links.p[[2]])[l2$query.id]
  l2$stack.id = win.u$grl.ix[l2$subject.id]
  l2$y.pos = ylim.stacks$xaxis.pos[l2$stack.id]
  l2$x.pos = mapply(function(x,y,z,a) (y-z)*a + x, x = window.segs.u[l2$subject.id,]$start, y = start(l2),
                    z = start(win.u[l2$subject.id]), a = window.segs.u$width[l2$subject.id] / width(win.u)[l2$subject.id])


  .fix.l = function(ll)
  {
    if (!is.null(links.feat))
      for (col in names(links.feat))
        values(ll)[, col] = links.feat[ll$query.id, col]

      # set up connectors
      if (is.null(ll$v))
        ll$v = y.gaps[ll$stack.id]/4
      else
        ll$v = y.gaps[ll$stack.id]*ll$v/2

      if (is.null(ll$h))
        ll$h = (window.segs.xlim$end[ll$stack.id] - window.segs.xlim$start[ll$stack.id])/20
      else
        ll$h = (window.segs.xlim$end[ll$stack.id] - window.segs.xlim$start[ll$stack.id])*ll$h

      if (is.null(ll$arrow))
        ll$arrow = TRUE

      if (is.null(ll$cex.arrow))
        ll$cex.arrow = 1

      if (is.null(ll$lwd))
        ll$lwd = 1


      if (is.null(ll$lty))
        ll$lty = 3


      if (is.null(ll$col))
        ll$col = 'red'

      if (is.null(ll$col.arrow))
        ll$col.arrow = ll$col

      return(ll)
  }

  if (length(l1)>0)
    l1 = .fix.l(l1)

  if (length(l2)>0)
    l2 = .fix.l(l2)


  ## now pair up / merge l1 and l2 using query.id as primary key
  if (length(l1)>0 & length(l2)>0)
  {
    pairs = merge(data.frame(l1.id = 1:length(l1), query.id = l1$query.id), data.frame(l2.id = 1:length(l2), query.id = l2$query.id))[, c('l1.id', 'l2.id')]
    l1.paired = GenomicRanges::as.data.frame(l1)[pairs[,1], ]
    l2.paired = GenomicRanges::as.data.frame(l2)[pairs[,2], ]
  }
  else
  {
    l2.paired = l1.paired = data.frame()
    pairs = data.frame()
  }


  ## some l1 and l2 will be unpaired
  l.unpaired = GRanges(seqlengths = GenomeInfoDb::seqlengths(links));
  p1 = p2 = c();
  if (nrow(pairs)>0)
  {
    p1 = pairs[,1]
    p2 = pairs[,2]
  }

  if (length(l1)>0)
    l.unpaired = c(l.unpaired, l1[setdiff(1:length(l1), p1)])

  if (length(l2)>0)
    l.unpaired = c(l.unpaired, l2[setdiff(1:length(l2), p2)])

  if (length(l.unpaired)>0)
  {
    l.unpaired$v = l.unpaired$v/4
    l.unpaired$h = l.unpaired$h/2
    l.unpaired$col = alpha(l.unpaired$col, 0.5)
    ## unpaired will get a "bridge to nowhere" in the proper direction (eg down)
    l.unpaired$y.pos = ylim.stacks$end[l.unpaired$stack.id]
    #                    l.unpaired$y.pos2 = l.unpaired$y.pos + top.gaps[l.unpaired$stack.id]
    l.unpaired$y.pos2 = l.unpaired$y.pos + l.unpaired$v

    connectors(l.unpaired$x.pos, l.unpaired$y.pos, as.character(GenomicRanges::strand(l.unpaired)),
               l.unpaired$x.pos, l.unpaired$y.pos2,
               as.character(GenomicRanges::strand(l.unpaired)),
               v = abs(l.unpaired$v), h = l.unpaired$h, type = 'S',
               f.arrow = l.unpaired$arrow, b.arrow = l.unpaired$arrow,
               cex.arrow = 0.2*l.unpaired$cex.arrow,
               col.arrow = l.unpaired$col.arrow,
               lwd = l.unpaired$lwd, lty = l.unpaired$lty, col = l.unpaired$col)

    if (!is.null(l.unpaired$label))
    {
      cex = 0.5;
      if (!is.null(l.unpaired$cex.label))
        cex = l.unpaired$cex.label

      #                        text(l.unpaired$x.pos, l.unpaired$y.pos2, l.unpaired$label, adj = c(0.5, 0.5), cex = cex)
      l.unpaired$text.y.pos =l.unpaired$y.pos - (ylim.stacks$end[l.unpaired$stack.id]-ylim.stacks$start[l.unpaired$stack.id])/100
      text(l.unpaired$x.pos, l.unpaired$text.y.pos, l.unpaired$label, adj = c(0.5, 1), cex = cex)
    }
  }

  ## now draw connectors for paired links
  ## pairs on the same level will get a "U" link,
  ## pairs on different levels will get "S" links

  ## fix y distances so that connections go from topmost part of bottom most connection to the bottom most part of
  ## top connection (ie the xaxis position)

  if (nrow(l1.paired)>0)
  {
    l1.paired$bottom = l1.paired$y.pos < l2.paired$y.pos

    if (any(l1.paired$bottom))
      l1.paired$y.pos[l1.paired$bottom] = ylim.stacks$end[l1.paired$stack.id[l1.paired$bottom]]

    if (any(!l1.paired$bottom))
      l2.paired$y.pos[!l1.paired$bottom] = ylim.stacks$end[l2.paired$stack.id[!l1.paired$bottom]]

    # also make all c type connections top connections with positive v
    ctype = ifelse(l1.paired$stack.id == l2.paired$stack.id, 'U', 'S')
    l1.paired$y.pos[ctype == 'U'] = ylim.stacks$end[l1.paired$stack.id[ctype == 'U']]
    l2.paired$y.pos[ctype == 'U'] = ylim.stacks$end[l2.paired$stack.id[ctype == 'U']]
    l1.paired$v[ctype == 'U'] = abs(l1.paired$v[ctype == 'U'])

    win.width = diff(par('usr')[1:2])
    l1.paired$v = l1.paired$v * affine.map(abs(l2.paired$x.pos - l1.paired$x.pos)/diff(par('usr')[1:2]), xlim = c(0, 1), ylim = c(0.5, 1)) ## make longer links taller
    connectors(l1.paired$x.pos, l1.paired$y.pos, l1.paired$strand, l2.paired$x.pos, l2.paired$y.pos, l2.paired$strand,
               v = l1.paired$v, h = l1.paired$h, type = ctype,
               f.arrow = l1.paired$arrow, b.arrow = l1.paired$arrow, cex.arrow = 0.2*l1.paired$cex.arrow, col.arrow = l1.paired$col.arrow,
               lwd = l1.paired$lwd, lty = l1.paired$lty, col = l1.paired$col)

    if (!is.null(l1.paired$label))
    {
      cex = 0.5;
      if (!is.null(l1.paired$cex.label))
        cex = l1.paired$cex.label

      ##                         text(l1.paired$x.pos, l1.paired$y.pos+l1.paired$v/2, l1.paired$label, adj = c(0.5, 0.5), cex = cex)
      ##                         text(l2.paired$x.pos, l2.paired$y.pos+l1.paired$v/2, l2.paired$label, adj = c(0.5, 0.5), cex = cex)

      l1.paired$text.y.pos = l1.paired$y.pos - (ylim.stacks$end[l1.paired$stack.id]-ylim.stacks$start[l1.paired$stack.id])/100
      l2.paired$text.y.pos = l2.paired$y.pos - (ylim.stacks$end[l2.paired$stack.id]-ylim.stacks$start[l2.paired$stack.id])/100

      text(l1.paired$x.pos, l1.paired$text.y.pos, l1.paired$label, adj = c(0.5, 1), cex = cex)
      text(l2.paired$x.pos, l2.paired$text.y.pos, l2.paired$label, adj = c(0.5, 1), cex = cex)
    }
  }

}
