
library(gTrack)
library(gUtils)
library(testthat)

context('gTrack tests')


gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))









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









## formatting

## clear

## dat

## colormap

## show

## .identical.seqinfo()

## plot()

## karyogram()

## track.gencode()

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

## test_that("testing track.gencode() works", {
##
## })



## draw.ranges()


## draw.grl()


## barbs()

## bernsteinp()


## connectors()


## affine.map()

test_that("testing affine.map() works", {

    expect_equal(affine.map(80), 0.5)
    expect_equal(affine.map(c(10, 20))[1], 0)
    expect_equal(affine.map(c(10, 20))[2], 1)

})



## alpha()

test_that("testing alpha() works", {

    expect_match(alpha('blue', 1), '#0000FFFF')
    
})





## blend()


## col.scale()


## brewer.master()


## lighten()



## plot.bank()


## draw.triangle()


## clip_polys()


## diamond()



## triangle()



## .geti()

test_that("testing .geti() works", {

    expect_error(.geti(gr, gr2))
    expect_equal(.geti(5, 10)$x, 7.5)
    expect_equal(.geti(5, 10)$y, 2.5)
    
})


## .getpoints()

test_that("testing .getpoints() works", {

    expect_true(is(.getpoints(c(1, 2, 3), c(5, 6, 7)), 'list'))
    expect_equal(as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$x[1]), 3)
    expect_equal(as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$x[2]), 3.5)
    expect_equal(as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$x[3]), 4)
    expect_equal(as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$x[4]), 3.5)
    expect_equal(as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$y[1]), 2) 
    expect_equal(as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$y[2]), 2.5) 
    expect_equal(as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$y[3]), 2) 
    expect_equal( as.numeric(.getpoints(c(1, 2, 3), c(5, 6, 7))$y[4]), 1.5) 

})




## .all.xpairs

test_that("testing .all.xpairs() works", {

    expect_equal(.all.xpairs(c(3, 5), c(3, 10))[1], 3)
    expect_equal(.all.xpairs(c(3, 5), c(3, 10))[2], 3)
    expect_equal(.all.xpairs(c(3, 5), c(3, 10))[3], 5)
    expect_equal(.all.xpairs(c(3, 5), c(3, 10))[4], 10)
    expect_equal(.all.xpairs(c(3, 5), c(3, 10))[5], 1)
    expect_equal(.all.xpairs(c(3, 5), c(3, 10))[6], 2)

})






## color.bar



## readData



test_that("check gr.flatmap()", {

    expect_true(is(gr.flatmap(example_genes, window=GRanges("1:1000000-1100000")), 'list'))
    foo = gr.flatmap(example_genes, window=GRanges("1:1000000-1100000"))
    expect_equal((foo$window.segs)$start, 1)
    expect_equal((foo$window.segs)$end, 100001)
    expect_equal((foo$grl.segs)$exonCount[1], 2)
    expect_equal((foo$grl.segs)$exonCount[2], 10)
    foobar = gr.flatmap(example_genes, window=GRanges("1:1000000-1100000"), gap=1, squeeze=TRUE) 
    expect_equal((foobar$window.segs)$start, 0)
    expect_equal((foobar$window.segs)$end, 1)
})





test_that("check gr.stripstrand()", {

    expect_match(as.character(strand(gr.stripstrand(gr2)))[1], '*')

})




## format_windows 


## prep_defaults_for_plotting


## extract_data_from_tmp_dat


## get_seqinfo


## enforce_max_ranges


## smooth_yfield


## format_yfield_limits



## draw_x_ticks




## clip.seg

test_that('clip.seg() works', {

    expect_equal(as.numeric(clip.seg(gr2, region = GRanges('1:10-15'))$chr), 1)
    expect_equal(as.numeric(clip.seg(gr2, region = GRanges('1:10-15'))$pos1), 10)
    expect_equal(as.numeric(clip.seg(gr2, region = GRanges('1:10-15'))$pos2), 14)

})



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



### track.straw



## listify()

test_that("testing listify() works", {

    expect_true(is(listify(gr2, length()), 'list'))
    
})




HI=300
WI=600
##gr <- GRanges(c(1,1), IRanges(c(5,6), width=1), strand=c("+","-"), A=c(5,8), B=c(3,2), seqinfo= gUtils::si)

grl <- GRangesList(list(GRanges(c(1,2), IRanges(c(6,10), width=1), strand=c("+","-"),seqinfo=gUtils::si)))
set.seed(137)
mdat <- matrix(sample(10, 4, replace=TRUE), ncol=2, nrow=2)
mdat[upper.tri(mdat)] <- mdat[lower.tri(mdat)]




test_that('show setMethod', {
    
    gt = gTrack(GRanges(1, IRanges(1,100)))

    yfield_showed_gt = show(gt)[[1]]

    expect_equal(yfield_showed_gt, NA)

})


test_that("karyogram method", {
    
    expect_equal(karyogram(TRUE)$y.field, NA)

    expect_equal(karyogram(hg19 = TRUE, bands = TRUE)$y.field, NA)

})


#test_that("multipleTracks", {
#    ## It is also possible to add multiple plots to the same window. Use the concatenation operator.
#    plot(c(gTrack(gr, edges = graph, stack.gap = 5), gTrack(gr, mdata = heatMap, stack.gap = 5)))
#})


##if (FALSE) {
##test_that("draw.paths", {
##    gene1 = sort(sample(gUtils::gr.tile(gUtils::parse.gr('1:1-5e3+'), 50), 5))
##    gene2 = rev(sort(sample(gUtils::gr.tile(gUtils::parse.gr('2:1-5e3-'), 50), 12)))
##    gene3 = sort(sample(gUtils::gr.tile(gUtils::parse.gr('3:1-5e3+'), 50), 8))
##
##    ##Create a column that keeps a counter of the exon number.
##
##    gene1$exon = 1:length(gene1)
##    gene2$exon = 1:length(gene2)
##    gene3$exon = 1:length(gene3)
##
##    ## Combine into GRangesList
##    grl = GRangesList(gene1 = gene1, gene2 = gene2, gene3 = gene3)
##
##    gt.genes = gTrack(grl)
##
##    ## Plot two graphs, one with and one without the draw.paths parameter.
##    fusion = GRangesList(c(grl$gene1[1:3], grl$gene2[5:9], grl$gene3[7:8]))
##    gt.fusion = gTrack(fusion, draw.paths = FALSE)
##    gt.fusion.o = gTrack(fusion, draw.paths = TRUE)
##
##    ## separating the windows for the graph.
##    win = gUtils::parse.gr(c('1:1-1e4', '2:1-1e4', '3:1-1e4'))
##    plot(c(gt.genes, gt.fusion, gt.fusion.o), win +1e3)
##})
##}

## need to fix d variable. 

##if (FALSE) {
##test_that("name", {
##    ## create sequences from chromosomes 1-3.
##    fake.genome = c('1'=1e4, '2'=1e3, '3'=5e3)
##    tiles = gr.tile(fake.genome, 1)
##
##    ## Choose 5 random indices. These indices will store the variants.
##    hotspots = sample(length(tiles), 5)
##
##    ## for each sequence, calculate the shortest distance to one of the hotspots.
##    d = pmin(Inf, values(distanceToNearest(tiles, tiles[hotspots]))$distance, na.rm = TRUE)
##    ## for sequences near the hotspots, the "prob" will be a higher positive number. It becomes smaller as it moves farther from the hotspot.
##    prob = .05 + exp(-d^2/10000)
##
##    ## sample 2000 of the sequences. the one nearer to the hotspots will "probably" be selected.
##    mut = sample(tiles, 2000, prob = prob, replace = TRUE)
##
##    ## graph with different degrees of stack.gap. The higher numeric supplied to stack.gap helps separate the data, visually.
##    gt.mut0 = gTrack(mut, circle = TRUE, stack.gap = 0, name = "Track 0")
##    gt.mut2 = gTrack(mut, circle = TRUE, stack.gap = 2, name = "Track 2")
##    gt.mut10 = gTrack(mut, circle = TRUE, stack.gap = 10, name = "Track 10")
##    gt.mut50 = gTrack(mut, circle = TRUE, stack.gap = 50, name = "Track 50")
##
##    win = si2gr(fake.genome)
##    plot(c(gt.mut0, gt.mut2, gt.mut10, gt.mut50), win)
##})
##}
