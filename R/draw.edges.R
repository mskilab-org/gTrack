draw.edges <- function(gr, edges, grl.segs, ylim.subplot, y.pad) {

if (data.table::is.data.table(edges))
  edges = as.data.frame(edges)

if (is.null(edges$col))
  edges$col = NA

if (is.null(edges$lty))
  edges$lty = NA

if (is.null(edges$lwd))
  edges$lwd = NA

if (is.null(edges$h))
  edges$h = 1

if (is.null(edges$cex.arrow))
  edges$cex.arrow = 1

edges.og = edges
edges = edges[edges$to %in% grl.segs$group | edges$from %in% grl.segs$group, ]

if (nrow(edges)>0)
{
  if (any(ix <- is.na(edges$col)))
    edges$col[ix] = 'black'

  if (any(ix <- is.na(edges$lwd)))
    edges$lwd[ix] = 1

  if (any(ix <- is.na(edges$lty)))
    edges$lty[ix] = 1

  ## now translate $from $to from group to gr id, and choosing last range in group for "from" indices
  ## and first range in group for "to" indices putting NA's when gr not in grl.segs

  first.ix = which(grl.segs$first)
  last.ix = which(grl.segs$last)

  grl.segs$to.gr = grl.segs$from.gr = 1:nrow(grl.segs)
  edges$edge.id = 1:nrow(edges)
  edges = merge(merge(edges, grl.segs[last.ix, c('group', 'from.gr')], by.x = 'from', by.y = 'group', all.x = T),
                grl.segs[first.ix, c('group', 'to.gr')], by.x = 'to', by.y = 'group', all.x = T)


  #                  edges$to.gr = first.ix[match(edges$to, grl.segs[first.ix, ]$group)]
  #                  edges$from.gr = last.ix[match(edges$from, grl.segs[last.ix, ]$group)]

  ## replicate edges that connect input gr's that have multiple instantiations

  ## figure out which grs are clipped
  grl.segs$clipped.start = !(start(gr)[grl.segs$query.id] == grl.segs$start)
  grl.segs$clipped.end = !(end(gr)[grl.segs$query.id] == grl.segs$end)

  clipped.to = ((grl.segs$clipped.start[edges$to.gr] & grl.segs$strand[edges$to.gr] == '+') |
                  (grl.segs$clipped.end[edges$to.gr] & grl.segs$strand[edges$to.gr] == '-'))
  clipped.from = ((grl.segs$clipped.start[edges$from.gr] & grl.segs$strand[edges$from.gr] == '-') |
                    (grl.segs$clipped.end[edges$from.gr] & grl.segs$strand[edges$from.gr] == '+'))

  clipped.to[is.na(clipped.to)] = F
  clipped.from[is.na(clipped.from)] = F
  ## now determine start and end points based on signs of connectors
  to.ix = !is.na(edges$to.gr) & !clipped.to
  from.ix = !is.na(edges$from.gr) & !clipped.from
  edges$x.pos.from = edges$x.pos.to = edges$y.pos.from = edges$y.pos.to = edges$dir.from = edges$dir.to =  NA

  #                  from.ix = !is.na(edges$fromx.gr)
  #                  to.ix = !is.na(edges$to.gr)

  if (any(to.ix)) ## connect to left end if + and right end if -
  {
    edges$x.pos.to[to.ix] = ifelse(grl.segs$strand[edges$to.gr[to.ix]] != '-',
                                   grl.segs$pos1[edges$to.gr[to.ix]], grl.segs$pos2[edges$to.gr[to.ix]])
    edges$y.pos.to[to.ix] = grl.segs$y[edges$to.gr[to.ix]]
    edges$dir.to[to.ix] = ifelse(grl.segs$strand[edges$to.gr[to.ix]] != '-',
                                 '-', '+')
  }

  if (any(from.ix)) ## connect from right end if + and left end if -
  {
    edges$x.pos.from[from.ix] = ifelse(grl.segs$strand[edges$from.gr[from.ix]] != '-',
                                       grl.segs$pos2[edges$from.gr[from.ix]], grl.segs$pos1[edges$from.gr[from.ix]])
    edges$y.pos.from[from.ix] = grl.segs$y[edges$from.gr[from.ix]]
    edges$dir.from[from.ix] = ifelse(grl.segs$strand[edges$from.gr[from.ix]] != '-',
                                     '+', '-')
  }


  ## now let NA endpoints "dangle"
  dangle.w = diff(par('usr')[1:2])/100
  edges$dangle = F

  if (is.null(edges$dangle.w))
    edges$dangle.w = 1

  edges$dangle.w  = edges$dangle.w * dangle.w

  if (any(!from.ix))
  {
    edges$x.pos.from[!from.ix] =
      edges$x.pos.to[!from.ix] + ifelse(grl.segs$strand[edges$to.gr[!from.ix]] != '-', -1, 1)*edges$dangle.w[!from.ix]
    edges$y.pos.from[!from.ix] = edges$y.pos.to[!from.ix]
    edges$from.gr[!from.ix] =  edges$from.gr[!from.ix]
    edges$dir.from[!from.ix] =  ifelse(edges$dir.to[!from.ix] != '-', '-', '+')
    edges$dangle[!from.ix] = T
  }

  if (any(!to.ix))
  {
    edges$x.pos.to[!to.ix] =
      edges$x.pos.from[!to.ix] + ifelse(grl.segs$strand[edges$from.gr[!to.ix]] != '-', 1, -1)*edges$dangle.w[!to.ix]
    edges$y.pos.to[!to.ix] = edges$y.pos.from[!to.ix]
    edges$to.gr[!to.ix] =  edges$to.gr[!to.ix]
    edges$dir.to[!to.ix] =  ifelse(edges$dir.from[!to.ix] != '-', '-', '+')
    edges$dangle[!to.ix] = T
  }

  edges$from.ix = from.ix
  edges$to.ix = to.ix
  edges = edges[order(edges$dangle), ]; ## hack to make sure we prefer non-dangling versions of edges when we have a choice (TODO: MAKE NON HACKY)
  dup = duplicated(cbind(edges$edge.id, edges$from.gr)) | duplicated(cbind(edges$edge.id, edges$to.gr))
  edges = edges[(edges$from.ix | edges$to.ix) & !dup, ]

  if (nrow(edges)>0)
  {

    if (!is.null(ylim.subplot))
      tmp.ydiff = diff(ylim.subplot)*(1-2*y.pad)
    else
      tmp.ydiff = diff(range(grl.segs$y, na.rm = T))

    #                      edges$type = ifelse(edges$y.pos.from == edges$y.pos.to, 'U', 'S')
    uthresh = max(grl.segs$ywid)
    edges$type = ifelse(abs(edges$y.pos.from - edges$y.pos.to) <= uthresh, 'U', 'S')

    edges$f.arrow = T
    edges$b.arrow = F
    if (is.null(edges$v))
      edges$v = 1


    edges$v = ifelse(abs(edges$y.pos.from - edges$y.pos.to)<=uthresh,
                     2*uthresh*affine.map(abs(edges$x.pos.to-edges$x.pos.from), xlim = c(0, diff(par('usr')[1:2])), ylim = c(0.5, 1.0)),
                     abs(edges$y.pos.from-edges$y.pos.to)/2)*edges$v
    #                  edges$v[is.na(edges$v)] = 0
    edges$h = dangle.w*edges$h * affine.map(abs(edges$x.pos.to-edges$x.pos.from), xlim = c(0, diff(par('usr')[1:2])), ylim = c(0.5, 1.0))
    make.flat = edges$type == 'U' & edges$dir.to != edges$dir.from

    if (!is.null(edges$not.flat)) ## allow edges specified as $not.flat to have height, unless dangling
      if (any(ix <- edges$not.flat & !edges$dangle))
        make.flat[ix] = F

    if (any( ix <- make.flat ))
    {
      edges$v[ix] = 0
      edges$h[ix] = 0
    }

    if (!is.null(edges$v.override))
      if (any(ix <- (!is.na(edges$v.override) & edges$to.ix & edges$from.ix)))
        edges$v[ix] = edges$v.override[ix]

    if (!is.null(edges$h.override))
      if (any(ix <- (!is.na(edges$h.override) & edges$to.ix & edges$from.ix)))
        edges$h[ix] = dangle.w*edges$h.override[ix]


    connectors(edges$x.pos.from, edges$y.pos.from, edges$dir.from,
               edges$x.pos.to, edges$y.pos.to, edges$dir.to,
               v = edges$v, h = edges$h, type = edges$type,
               f.arrow = T, b.arrow = F,
               cex.arrow = 0.2*edges$cex.arrow,
               col.arrow = edges$col,
               lwd = edges$lwd, lty = edges$lty, col = edges$col)

    if (!is.null(edges$label))
    {
      text((edges$x.pos.from + edges$x.pos.to)/2, edges$y.pos.from + 0.5*edges$v*sign(edges$y.pos.to - edges$y.pos.from + 0.01),
           edges$label, adj = c(0.5, 0.5), col = 'black')
    }
  }
}

}
