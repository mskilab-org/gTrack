library(gTrack)
library(gUtils)
context("gTrack tests")

HI=300
WI=600
gr <- GRanges(c(1,1), IRanges(c(5,6), width=1), strand=c("+","-"), A=c(5,8), B=c(3,2), seqinfo= gUtils::si)

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
gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3"), c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19), end = c(2,4,6,8,10,12,14,16,18,20),  names = head(letters,10)), GC=seq(1,10,length=10), name=seq(5,10,length=10))

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
    
test_that("edges",  {   
    plot(gTrack(gr , edges = graph , stack.gap = 5))
})

test_that("stack.gap", {
    plot(gTrack(gr , stack.gap = 2))
})

test_that("y.field", {
    plot(gTrack(gr , y.field = 'GC')) 
})

test_that("bars", {
    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = 'light blue'))
})

test_that("lines", {
    plot(gTrack(gr , y.field = 'GC' , lines = TRUE , col = 'purple'))
})

test_that("cirlces", {
    plot(gTrack(gr , y.field = 'GC' , circles = TRUE , col = 'magenta' , border = '60'))
})

test_that("colorfield", {
    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , colormaps = list(GC = c("1"="red" , "2" = "blue" , "3"="magenta", "4"="light blue" ,"5"="black" , "6"="green", "7"="brown" , "8"="pink", "9"="yellow", "10" = "orange")) ))
})

test_that("gr.colorfield", {
    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC'))
})

test_that("gr.labelfield", {
    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC' , gr.labelfield = 'name'))
})

test_that("col", {
    graph = data.frame(from = 1:9, to = c(6,9,7,2,4,10,8,5,3) , col = c('red', 'blue', 'green'))
    plot(gTrack(gr , edges = graph , stack.gap = 5))  
})

if(FALSE){

test_that("lwd", {
    graph$lwd = 1.844941
    plot(gTrack(gr, edges = graph, stack.gap = 5))
})

test_that("lty", {
    graph$lty = c(1,2,3)
    plot(gTrack(gr , edges = graph , stack.gap = 5))
})

test_that("h", {
    graph$h = 10 
    plot(gTrack(gr , edges = graph , stack.gap = 5))
})

test_that("mdata",{
    heatMap = matrix(runif(length(gr)^2), nrow = 10, ncol = 10)
    plot(gTrack(gr, mdata = heatMap, stack.gap = 5))   
})

test_that("multipleTracks", {
    ## It is also possible to add multiple plots to the same window. Use the concatenation operator.
    plot(c(gTrack(gr, edges = graph, stack.gap = 5), gTrack(gr, mdata = heatMap, stack.gap = 5)))
})


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

## test_that("track.gencode", {
##     gt.ge = track.gencode()
## })
