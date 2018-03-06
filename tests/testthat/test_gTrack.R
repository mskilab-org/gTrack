library(gTrack)
library(gUtils)

library(testthat)

context('gTrack tests')


gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))




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


test_that("initialize setMethod", {

    expect_error(gTrack(1), "check input: gTrack objects can only be defined around GRanges, GRangesLists, RleLists, ffTrack, file paths to .rds files of the latter object types, or file paths to  UCSC format files")

    expect_error(gTrack(GRanges(1, IRanges(1, 100)), mdata = 1), "optional arg mdata be either empty matrix or a square matrix of same dimensions as data GRanges")

    sapply(list(matrix(c(2,3,4,2), nrow = 2, ncol = 2, byrow = TRUE), matrix(c(2,3,4,2), nrow = 2, ncol = 2, byrow = TRUE)), function(x) is.array(x) | inherits(x, "Matrix"))


    gtrack_y_field = gTrack(gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3") , c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) , end = c(2,4,6,8,10,12,14,16,18,20), names = head(letters,10)), GC=seq(1,10,length=10), name=seq(5,10,length=10)), y.field = 'GC')

    expect_equal(gtrack_y_field$y.field, 'GC')
    
})


test_that("[ method", {

    gt = gTrack(GRanges(1, IRanges(1,100)))

    expect_equal(gt[1]$height, 10)

    #expect_equal(gt[2]$height, NA)
    
})



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



test_that("c method", {

    gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3") ,c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) ,end = c(2,4,6,8,10,12,14,16,18,20), names = head(letters,10)),GC=seq(1,10,length=10), name=seq(5,10,length=10))

    heatMap = matrix(runif(length(gr)^2), nrow = 10, ncol = 10)

    ##c(gTrack(gr, edges = graph, stack.gap = 5), gTrack(gr, mdata = heatMap, stack.gap = 5))
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


test_that('mdata() function', {
    
    heatMap = matrix(runif(length(gr)^2), nrow = 10, ncol = 10)
    gTrack_heatMap = gTrack(gr, mdata = heatMap)

    gTrack_heatMap_matrices = mdata(gTrack_heatMap)

    expect_is(gTrack_heatMap_matrices[[1]], "matrix")

    expect_equal(mdata(gTrack_heatMap, GRanges()), NULL)

    expect_is(mdata(gTrack_heatMap, GRanges(1, IRanges(1,100))), "matrix")
    
    # si2gr(gTrack_heatMap)
    
})

test_that('reduce() function', {

    gt = gTrack(GRanges(1, IRanges(1,100)))    
    expect_identical(reduce(gt), GRanges(1, IRanges(1, 100)))
    
})


test_that('show setMethod', {
    
    gt = gTrack(GRanges(1, IRanges(1,100)))

    yfield_showed_gt = show(gt)[[1]]

    expect_equal(yfield_showed_gt, NA)

})


test_that("karyogram method", {
    
    expect_equal(karyogram(TRUE)$y.field, NA)

    expect_equal(karyogram(hg19 = TRUE, bands = TRUE)$y.field, NA)

    karyogram(FALSE)

    #karyogram(arms = TRUE)
    
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



## plot
## function(x,   y, windows = si2gr(seqinfo(x)), ## windows to plot can be Granges or GRangesList
##                    links = NULL, ## GRangesList of pairs of signed locations,
##                    gap = NULL,  ## spacing betwen windows (in bp)
##                    y.heights = NULL, # should be scalar or length(windows) if windows is a GRangesList
##                    y.gaps = NULL, # relative heights of gaps below and above stacks (xaxes will be drawn here)
##                    cex.xlabel = 1,
##                    cex.ylabel = 1,
##                    max.ranges = NA, # parameter for max ranges to draw on canvas in each track (overrides formatting)
##                    links.feat = NULL, # links features override for links (must be nrow 1 or length(links) data frame
##                    verbose=FALSE,
##                    legend.params = list(),
##                    ... ## additional args to draw.grl OR last minute formatting changes to gTrack object
##                   )



### karyogram 
### karyogram = function(hg19 = TRUE, bands = TRUE, arms = TRUE, tel.width = 2e6, ... )

test_that('karyogram() works', {
    
    ## default
    expect_equal(karyogram()$height, 10)
    expect_equal(karyogram()$angle, 10)
    expect_equal(karyogram()$y.quantile, 0.01)
    expect_equal(karyogram()$cex.label, 1)
    expect_equal(karyogram()$gr.cex.label.gr, 0.8)
    ## hg19
    expect_equal(karyogram(hg19 = FALSE)$height, 10)
    expect_equal(karyogram(hg19 = FALSE)$angle, 10)
    expect_equal(karyogram(hg19 = FALSE)$y.quantile, 0.01)
    expect_equal(karyogram(hg19 = FALSE)$cex.label, 1)
    expect_equal(karyogram(hg19 = FALSE)$gr.cex.label.gr, 0.8)
    ## bands
    ## nobands = karyogram(bands = FALSE)
    ## Error in validObject(.Object) : 
    ##     invalid class “GRanges” object: 'seqnames(x)' contains missing values
    expect_equal(karyogram(arms = FALSE)$height, 10)
    expect_equal(karyogram(arms = FALSE)$angle, 10)
    expect_equal(karyogram(arms = FALSE)$y.quantile, 0.01)
    expect_equal(karyogram(arms = FALSE)$cex.label, 1)
    expect_equal(karyogram(arms = FALSE)$gr.cex.label.gr, 0.8)   
    ### tel.width = 2e6
    expect_equal(karyogram(tel.width = 1e3)$height, 10)
    expect_equal(karyogram(tel.width = 1e3)$angle, 10)
    expect_equal(karyogram(tel.width = 1e3)$y.quantile, 0.01)
    expect_equal(karyogram(tel.width = 1e3)$cex.label, 1)
    expect_equal(karyogram(tel.width = 1e3)$gr.cex.label.gr, 0.8)   

})




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




## brewer.master


## track.straw


### seg.on.seg()

test_that('seg.on.seg() works', {

    expect_true(seg.on.seg(dt, dt)[1])
    expect_true(seg.on.seg(dt, dt)[2])
    expect_true(seg.on.seg(dt, dt)[3])

})



## dedup()

test_that('dedup() works', {

    expect_equal(dedup(c(rep(2, 10.5), rep(3, 20)))[30], "3.20")

})



## listify()

test_that("testing listify() works", {

    expect_true(is(listify(gr2, length()), 'list'))
    
})






