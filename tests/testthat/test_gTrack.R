library(gTrack)
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

test_that("edges",  {   
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

