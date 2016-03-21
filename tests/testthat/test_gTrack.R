library(gTrack)
library(gUtils)
context("gTrack tests")

HI=300
WI=600
gr <- GRanges(c(1,1), IRanges(c(5,6), width=1), strand=c("+","-"), A=c(5,8), B=c(3,2), seqinfo= si)
gr <- GRanges(c(1), IRanges(c(5), width=1), strand=c("+"), A=c(5), B=c(3), seqinfo= si)

grl <- GRangesList(list(GRanges(c(1,2), IRanges(c(6,10), width=1), strand=c("+","-"),seqinfo=si)))
mdat <- matrix(sample(10, 4, replace=TRUE), ncol=2, nrow=2)
mdat[upper.tri(mdat)] <- mdat[lower.tri(mdat)]

test_that("gTrack receives input vals into format", {
  g <- gTrack(gr, height=7, xaxis.chronly = TRUE, sep.bg.col='red')
  #png("rtdocs/figures/basic_background_col.png", height=HI, width=WI);
  #gTrack::plot(g)
  #dev.off()

  g <- gTrack(gr, xaxis.chronly = TRUE, mdata=mdat)
  g2 <- gTrack(gr, height=7, xaxis.chronly = TRUE)

  #png("rtdocs/figures/basic_background_col.png", height=HI, width=WI);
  ##plot(g2, GRanges(1, IRanges(4,6)))
  #dev.off()
})
