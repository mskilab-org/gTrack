library(gTrack)
library(gUtils)

library(testthat)

context('gTrack tests')

HI=300
WI=600
##gr <- GRanges(c(1,1), IRanges(c(5,6), width=1), strand=c("+","-"), A=c(5,8), B=c(3,2), seqinfo= gUtils::si)

grl <- GRangesList(list(GRanges(c(1,2), IRanges(c(6,10), width=1), strand=c("+","-"),seqinfo=gUtils::si)))
set.seed(137)
mdat <- matrix(sample(10, 4, replace=TRUE), ncol=2, nrow=2)
mdat[upper.tri(mdat)] <- mdat[lower.tri(mdat)]

#test_that("gTrack receives input vals into format", {
  #g <- gTrack(gr, height=7, xaxis.chronly = TRUE, y.field="A", bars=TRUE, xaxis.suffix = "bp")
  #png("rtdocs/figures/basic_background_col.png", height=HI, width=WI);
  #plot(g,windows=reduce(g))
  #dev.off()

  #g <- gTrack(gr, xaxis.chronly = TRUE, mdata=mdat)
  #g2 <- gTrack(gr, height=7, xaxis.chronly = TRUE)

  #png("rtdocs/figures/basic_background_col.png", height=HI, width=WI);
  ##plot(g2, GRanges(1, IRanges(4,6)))
  #dev.off()
#})

##create a GRanges object storing 10 sequences. These sequences will serve as nodes for the graph.
gr = GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3"), c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19), end = c(2,4,6,8,10,12,14,16,18,20),  names = head(letters,10)), GC=seq(1,10,length=10), name=seq(5,10,length=10))

## Specify links between nodes using a matrix. Numeric 1s refer to a connection while conversely with 0s.
##create an N*N matrix filled with 0s.
graph = matrix(0 , nrow = 10 , ncol = 10)
##set certain indices to 1.
graph[1,3]=1
graph[1,10]=1
graph[2,5]=1
graph[2,8]=1
graph[3,5]=1
graph[4,1]=1
graph[4,2]=1
graph[4,6]=1
graph[4,9]=1
graph[5,1]=1
graph[5,2]=1
graph[5,4]=1
graph[8,1]=1
graph[8,2]=1
graph[9,1]=1
graph[10,1]=1
    




#### gTrack testing


##    new('gTrack', data = data, y.field = y.field, mdata = mdata, name = name, format = formatting,
##      edges = edges, vars = vars, draw.paths = draw.paths, colormaps = colormaps, height = height, ygap = ygap,
##      stack.gap = stack.gap, col = col, border = border, angle = angle, draw.var = draw.var,
##      gr.colorfield = gr.colorfield, y.quantile = y.quantile,
##      cex.label = cex.label, gr.cex.label.gr = gr.cex.label, gr.srt.label = gr.srt.label,
##      y.cap = y.cap, lwd.border = lwd.border, hadj.label = hadj.label, vadj.label = vadj.label, smooth = smooth,
##      round = round, ywid = ywid, ypad = ypad, seqinfo = seqinfo, circles = circles, lines = lines,
##      bars = bars, triangle = triangle, ylab = ylab, max.ranges = max.ranges, source.file.chrsub = source.file.chrsub,
##      y0.bar = y0.bar, yaxis = yaxis, yaxis.pretty = yaxis.pretty, yaxis.cex = yaxis.cex,
##      chr.sub = chr.sub, edgevars = edgevars, gr.labelfield = gr.labelfield,
##      grl.labelfield = grl.labelfield, xaxis.prefix = xaxis.prefix, xaxis.unit = xaxis.unit,
##      xaxis.suffix = xaxis.suffix, xaxis.round = xaxis.round, xaxis.interval = xaxis.interval,
##      xaxis.cex.label = xaxis.cex.label, xaxis.newline = xaxis.newline,
##      xaxis.chronly = xaxis.chronly, xaxis.width = xaxis.width,
##      labels.suppress = labels.suppress, labels.suppress.gr = labels.suppress.gr, labels.suppress.grl = labels.suppress.grl,
##      xaxis.label.angle = xaxis.label.angle, xaxis.ticklen = xaxis.ticklen,
##      xaxis.cex.tick = xaxis.cex.tick, sep.lty = sep.lty, sep.lwd = sep.lwd, sep.bg.col = sep.bg.col,
##      sep.draw = sep.draw, y0 = y0, y1 = y1, m.sep.lwd = m.sep.lwd, m.bg.col = m.bg.col,
##      cmap.min = cmap.min, cmap.max = cmap.max, bg.col = bg.col)



## test_that('gTrack(), edges ',  {   
##     plot(gTrack(gr , edges = graph , stack.gap = 5))
## })




test_that('gTrack(), stack.gap', {

    foo = gTrack(gr , stack.gap = 2)
    expect_equal(foo$stack.gap, 2)   ## test here
    expect_equal(as.logical(foo$y.field), NA)
    expect_equal(foo$angle, 15)
    expect_equal(foo$draw.paths, FALSE)
    expect_equal(foo$height, 10)
    expect_equal(foo$ygap, 2)

})




test_that('gTrack(), y.field', {

    foo = gTrack(gr , y.field = 'GC') 
    expect_equal(foo$y.field, 'GC')  ## test here
    expect_equal(foo$draw.paths, FALSE)
    expect_equal(foo$height, 10)
    expect_equal(foo$ygap, 2)
    expect_equal(foo$stack.gap, 0)

})



test_that('gTrack(), bars', {
    
    foo = gTrack(gr , y.field = 'GC' , bars = TRUE , col = 'light blue')
    expect_equal(foo$y.field, 'GC')  ## test here
    expect_equal(foo$draw.paths, FALSE)
    expect_equal(foo$bars, TRUE)
    expect_equal(foo$col, 'light blue')

})




test_that('gTrack(), lines', {
    
    foo = gTrack(gr , y.field = 'GC' , lines = TRUE , col = 'purple')
    expect_equal(foo$y.field, 'GC')  ## test here
    expect_equal(foo$draw.paths, FALSE)
    expect_equal(foo$bars, FALSE)
    expect_equal(foo$lines, TRUE)
    expect_equal(foo$col, 'purple')

})




test_that('gTrack(), circles', {
    
    foo = gTrack(gr , y.field = 'GC' , circles = TRUE , col = 'magenta' , border = '60')  ## why is 'border' a character string? 
    expect_equal(foo$y.field, 'GC')  ### test
    expect_equal(as.logical(foo$draw.paths), FALSE) 
    expect_equal(foo$height, 10)  
    expect_equal(foo$ygap, 2)
    expect_equal(foo$stack.gap, 0) 
    expect_equal(foo$col, 'magenta')   ## test
    expect_equal(foo$border, '60')  ## test
    expect_equal(foo$angle, 15)
    expect_equal(as.logical(foo$gr.labelfield), NA)

})





test_that('gTrack(), colorfield', {
    
    foo = gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , colormaps = list(GC = c("1"="red" , "2" = "blue" , "3"="magenta", "4"="light blue" ,"5"="black" , "6"="green", "7"="brown" , "8"="pink", "9"="yellow", "10" = "orange")) )
    expect_equal(foo$y.field, 'GC')  ### test
    expect_equal(as.logical(foo$draw.paths), FALSE) 
    expect_equal(foo$height, 10)  
    expect_equal(foo$ygap, 2)
    expect_equal(foo$stack.gap, 0) 
    expect_equal(as.logical(foo$col), NA)   ## test
    expect_equal(as.logical(foo$gr.colorfield), NA)  
    expect_equal(foo$angle, 15)
    expect_equal(as.logical(foo$gr.labelfield), NA)
    expect_equal(foo$bars, TRUE)  ## test
    ## how to test 'colormaps'?


})





test_that('gTrack(), gr.colorfield', {
    
    foo = gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC')
    expect_equal(foo$y.field, 'GC')  ### test
    expect_equal(as.logical(foo$draw.paths), FALSE) 
    expect_equal(foo$height, 10)  
    expect_equal(foo$ygap, 2)
    expect_equal(foo$stack.gap, 0) 
    expect_equal(as.logical(foo$col), NA)   ## test
    expect_equal(foo$gr.colorfield, 'GC')  ## test
    expect_equal(foo$angle, 15)
    expect_equal(as.logical(foo$gr.labelfield), NA)
    expect_equal(foo$bars, TRUE)

})





test_that('gTrack(), gr.labelfield', {
    
    foo = gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC' , gr.labelfield = 'name name _ name')
    expect_equal(foo$y.field, 'GC')
    expect_equal(as.logical(foo$draw.paths), FALSE) 
    expect_equal(foo$height, 10)  
    expect_equal(foo$ygap, 2)
    expect_equal(foo$stack.gap, 0) 
    expect_equal(as.logical(foo$col), NA)   
    expect_equal(foo$ gr.colorfield, 'GC')  
    expect_equal(foo$angle, 15)
    expect_equal(foo$gr.labelfield, 'name name _ name')

})





test_that('gTrack(), col', {
    
    graph = data.frame(from = 1:9, to = c(6,9,7,2,4,10,8,5,3) , col = c('red', 'blue', 'green'))
    foo = gTrack(gr, edges = graph , stack.gap = 5)  
    expect_equal(as.logical(foo$y.field), NA)
    expect_equal(as.logical(foo$draw.paths), FALSE) 
    expect_equal(foo$height, 10)  
    expect_equal(foo$ygap, 2)
    expect_equal(foo$angle, 15)
    ## test edges = graph
    expect_equal(foo$stack.gap, 5)  ## test here

})


#test_that("lwd", {
#    graph$lwd = 1.844941
#    gTrack_plot = gTrack(gr, edges = graph, stack.gap = 5)
#})

#test_that("lty", {
#    graph$lty = c(1,2,3)
#    plot(gTrack(gr , edges = graph , stack.gap = 5))
#})

#test_that("h", {
#    graph$h = 10 
#    plot(gTrack(gr , edges = graph , stack.gap = 5))
#})

test_that('gTrack(), mdata' ,{
    heatMap = matrix(runif(length(gr)^2), nrow = 10, ncol = 10)
    gTrack_heatMap = gTrack(gr, mdata = heatMap, stack.gap = 5)
    expect_equal(gTrack_heatMap$stack.gap, 5)
    expect_equal(gTrack_heatMap$triangle, TRUE)
})












#test_that("multipleTracks", {
#    ## It is also possible to add multiple plots to the same window. Use the concatenation operator.
#    plot(c(gTrack(gr, edges = graph, stack.gap = 5), gTrack(gr, mdata = heatMap, stack.gap = 5)))
#})

if (FALSE) {
test_that("draw.paths", {
    gene1 = sort(sample(gUtils::gr.tile(gUtils::parse.gr('1:1-5e3+'), 50), 5))
    gene2 = rev(sort(sample(gUtils::gr.tile(gUtils::parse.gr('2:1-5e3-'), 50), 12)))
    gene3 = sort(sample(gUtils::gr.tile(gUtils::parse.gr('3:1-5e3+'), 50), 8))

    ##Create a column that keeps a counter of the exon number.

    gene1$exon = 1:length(gene1)
    gene2$exon = 1:length(gene2)
    gene3$exon = 1:length(gene3)

    ## Combine into GRangesList
    grl = GRangesList(gene1 = gene1, gene2 = gene2, gene3 = gene3)

    gt.genes = gTrack(grl)

    ## Plot two graphs, one with and one without the draw.paths parameter.
    fusion = GRangesList(c(grl$gene1[1:3], grl$gene2[5:9], grl$gene3[7:8]))
    gt.fusion = gTrack(fusion, draw.paths = FALSE)
    gt.fusion.o = gTrack(fusion, draw.paths = TRUE)

    ## separating the windows for the graph.
    win = gUtils::parse.gr(c('1:1-1e4', '2:1-1e4', '3:1-1e4'))
    plot(c(gt.genes, gt.fusion, gt.fusion.o), win +1e3)
})
}

## need to fix d variable. 

if (FALSE) {
test_that("name", {
    ## create sequences from chromosomes 1-3.
    fake.genome = c('1'=1e4, '2'=1e3, '3'=5e3)
    tiles = gr.tile(fake.genome, 1)

    ## Choose 5 random indices. These indices will store the variants.
    hotspots = sample(length(tiles), 5)

    ## for each sequence, calculate the shortest distance to one of the hotspots.
    d = pmin(Inf, values(distanceToNearest(tiles, tiles[hotspots]))$distance, na.rm = TRUE)
    ## for sequences near the hotspots, the "prob" will be a higher positive number. It becomes smaller as it moves farther from the hotspot.
    prob = .05 + exp(-d^2/10000)

    ## sample 2000 of the sequences. the one nearer to the hotspots will "probably" be selected.
    mut = sample(tiles, 2000, prob = prob, replace = TRUE)

    ## graph with different degrees of stack.gap. The higher numeric supplied to stack.gap helps separate the data, visually.
    gt.mut0 = gTrack(mut, circle = TRUE, stack.gap = 0, name = "Track 0")
    gt.mut2 = gTrack(mut, circle = TRUE, stack.gap = 2, name = "Track 2")
    gt.mut10 = gTrack(mut, circle = TRUE, stack.gap = 10, name = "Track 10")
    gt.mut50 = gTrack(mut, circle = TRUE, stack.gap = 50, name = "Track 50")

    win = si2gr(fake.genome)
    plot(c(gt.mut0, gt.mut2, gt.mut10, gt.mut50), win)
})
}



### karyogram 
### karyogram = function(hg19 = TRUE, bands = TRUE, arms = TRUE, tel.width = 2e6, ... )





### track.gencode = function(gencode = NULL,
###     gene.collapse = TRUE,
###     genes = NULL,
###     grep = NULL,
###     grepe = NULL, ## which to exclude
###     bg.col = alpha('blue', 0.1), cds.col = alpha('blue', 0.6), utr.col = alpha('purple', 0.4),
###     st.col = 'green',
###     en.col = 'red',
###     grl.labelfield, ## Don't touch these
###     gr.labelfield,
###     col,
###     cached = T, ## if true will use cached version
###     cached.dir = Sys.getenv('GENCODE_DIR'),
###     cached.path = paste(cached.dir, "gencode.composite.rds", sep = '/'),  ## location of cached copy
###     cached.path.collapsed = paste(cached.dir, "gencode.composite.collapsed.rds", sep = '/'),  ## location of cached copy
###     gr.srt.label = 0,
###     gr.cex.label = 0.3,
###     cex.label = 0.5,
###     labels.suppress.gr = FALSE,
###     drop.rp11 = TRUE,
###     stack.gap = 1e6,
###     ...)
### 



test_that("track.gencode", {
    expect_error(track.gencode(), NA) ## check it runs correctly
})



## Draws (genomic, Iranges, Rle, or df with start and end field) intervals
## as rectangles of thickness "lwd", col "col".

## draw.ranges
## draw.ranges = function(x, y = NULL, lwd = 0.5, col = "black", border = col, label = NULL,
##                        draw.backbone = F, # if TRUE lines will be drawn linking all ranges at a single y-level +/- (a single y x group level if group is specified),
##                        group = NULL, # a vector of group labels same size as x (char or numeric), non-NA groups with more than 1 member will be connected by backbone (if draw.backbone = T)
##                        col.backbone = col[1],
##                        cex.label = 1,
##                        srt.label = 45,
##                        adj.label = c(0.5, 0),
##                        angle, # vector of angles (0 = rectangle -45 leftward barb, +45 rightward barb
##                        circles = FALSE, # if TRUE will draw circles at range midpoints instead of polygons
##                        bars = FALSE, # if TRUE will draw bars from y0.bar to the y locations
##                        points = NA, # if integer then will draw points at the midpoint with pch value "points"
##                        lines = FALSE, # if TRUE will connect endpoints of ranges with lines
##                        y0.bar = NULL, # only relevant if bars = T, will default to pars('usr')[3] if not specified
##                        strand.dir = T,  # if so, will interpret strand field as "direction sign" for interval, + -> right, - -> left, only makes difference if draw.pencils = T
##                        xlim = NULL, # only relevant if new.plot = T, can be ranges object
##                        ylim = NULL, # only relevant if new.plot = T
##                        clip = NULL, # GRanges / IRanges or 1 row df with $start, $end denoting single region to clip to .. will not draw outside of this region
## 
##                        lwd.border = 1,
##                        lwd.backbone = 1,
##                        verbose=FALSE,
##                        ...)




## draw.grl
## draw.grl = function(grl,
##                     y = NULL,  # can be either vector of length grl, or data.frame row / list with fields $start and $end
##                     # specifying y coordinates to "pack" the granges into (or just length 2 list)
##                     # note this is different from ylim, which determines the size of the canvas
##                     ylim = NULL,  # if null and y provided, will find range of y and plot
##                     ywid = NULL,
##                     edges = NULL, ## data.frame specifying edges connecting items of grl, with fields $from $to indexing grl and optional fields $lwd, $col, $lty specifying formatting options for connectors, for gr the from arrow will leave the right side of a + range and left side of a - range, and enter the left side of a + range and right side of a - range.   For grl, from arrows will be drawn from the last range in the "from" list item to the first range in the "to" list item
##                     draw.paths = F, # if draw.paths = T will treat grl's as paths / contigs,
##                     # connecting intervals joined by arrowed curved arcs using alternate stacking algo (ignoring y information)
## 
##                     draw.var = F, # if true, then varbase will be called on grl or unlist(grl)
##                     var = NULL, # optional GRangesList of same length as grl specifying variant bases with value cols $type, $varbase (ie output of varbase)
##                     # corresponding to each grl item
##                     var.col = NULL, # named vector with variant colors to override - valid names are XA, XT, XG, XC, S, D, I
##                     var.soft = T, ## whether to draw soft clips for varbase
##                     windows,  # windows specifies windows, can have optional meta data features $col and $border
##                     win.gap = NULL,
##                     stack.gap,
##                     min.gapwidth = 1, # only applies if windows is not specified
##                     col = NULL, border = NA, # all of these variables can be scalar or vectors of length(grl),
##                     # can be over-ridden by values / columns in grl
##                     col.backbone = alpha('gray', 0.8),
##                     gr.colorfield = NULL, ## values field in the gr from which colors can be mapped
##                     gr.colormap = NULL, ## named vector mapping fields in the gr.colorfield to colors, if unspecified brewer.master() will be applied
##                     gr.labelfield = NULL, ## field of gr labels to draw.
##                     grl.labelfield = NULL, ## field of grl to draw as label
##                     leg.params,
##                     labels = NULL, # vector of length(grl)
##                     labels.suppress = F,
##                     labels.suppress.grl = labels.suppress,
##                     labels.suppress.gr = labels.suppress,
##                     spc.label = 0.05, # number between 0 and 1 indicating spacing of label
##                     adj.label = c(0, 1),
##                     cex.label = 1,
##                     gr.cex.label = 0.8 * cex.label,
##                     gr.srt.label = 0,
##                     gr.adj.label = c(0,0.5),
##                     new.plot, new.axis, 
##                     sep.lty = 2,
##                     sep.lwd = 1,
##                     sep.bg.col = 'gray95',
##                     sep.draw = TRUE,
##                     y.pad,  # this is the fractional padding to put on top and bottom of ranges if y is specified as $start and $end pair (def was 0.05)
##                     xaxis.prefix = '', xaxis.suffix = 'MB', xaxis.unit = 1, xaxis.round = 3,
##                     xaxis.interval = 'auto', xaxis.pos = 1,
##                     xaxis.pos.label, xaxis.cex.label,
##                     xaxis.newline = FALSE,
##                     xaxis.chronly = FALSE,
##                     xaxis.ticklen = 1,
##                     xaxis.width = TRUE,
##                     xaxis.label.angle = 0,
##                     xaxis.cex.tick = 1,
##                     ylim.grid = ylim, # can be also specified directly for plots with multiple axes and/or non-numeric tracks
##                     y.grid = NA, # if non NA, then the number specifies the spacing between y grid (drawn from ylim[1] to ylim[2]), can also be named vector mapping grid.tick labels to plot locations
##                     ylab = NULL,
##                     y.grid.col = alpha('gray', 0.5),
##                     y.grid.pretty = 5,
##                     y.grid.cex = 1,
##                     y.grid.lty = 2,
##                     y.grid.lwd = 1,
##                     path.col = 'black',
##                     path.col.arrow = path.col,
##                     path.cex.arrow = 1,
##                     path.stack.y.gap = 1,
##                     path.stack.x.gap = 0,
##                     path.cex.v = 1,
##                     path.cex.h = 1,
##                     draw.backbone = NULL,
##                     xlim = c(0, 20000), # xlim of canvas
##                     points = NA, ## if non NA then will draw a given point with pch style
##                     circles = FALSE, ## only one of these should be true, however if multiple are true then they will be interpreted in this order
##                     bars = FALSE,
##                     y0.bar = NULL,
##                     lines = F,
##                     angle, # angle of barbs to indicate directionality of ranges
##                     verbose=FALSE,
##                     triangle=FALSE, # trigger a triangle matrix plot
##                     ylim.parent=NULL, ## ylim of the full thing. This is importat for angle preseveration
##                     legend.params = list(plot=TRUE),
##                     bg.col = NA, ## background of whole plot
##                     ...)
## 






## bernsteinp
## 



## connectors
## connectors = function(x0, y0, s0 = 1, x1, y1, s1 = 1, v = 0.1, h = 0.1,
##                      type = "S", f.arrow = T, b.arrow = F, nsteps = 100,
##                      cex.arrow = 1, col.arrow = 'black', lwd = 1, lty = 1, col = 'black')




## affine.map
## affine.map = function(x, ylim = c(0,1), xlim = c(min(x), max(x)), cap = F, cap.min = cap, cap.max = cap, clip = T, clip.min = clip, clip.max = clip)


## alpha
## alpha = function(col, alpha)



## 
## Blends colors
##



## col.scale
## col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white', invert = F # if T flips rgb.min and rgb.max
## 


## brewer.master


## lighten

## plot.blank


## draw.triangle
## function(grl,mdata,y, ylim.parent=NULL, windows = NULL, win.gap = NULL, m.bg.col = NA,
##         sigma = NA, ## if not NA then will blur with a Gaussian filter using a sigma value of this many base pairs
##         col.min='white', col.max='red', gr.colormap = NA, cmap.min = NA, cmap.max = NA, m.sep.lwd = NA, legend = TRUE, leg.params = list(),
##         min.gapwidth = 1,islog = FALSE) 
## 


## clip_polys <- function(dt, y0, y1)
## 


## diamond <- function (x11, x12, x21, x22, y0, y1, col=NULL) {


## triangle <- function(x1, x2, y, y0, y1, col=NULL) 




