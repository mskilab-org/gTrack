library(gTrack)
context("gTrack tests")

gr <- GRanges(c(1,2), IRanges(c(5,7), width=5), strand=c("+","-"))

test_that("gTrack receives input vals into format", {

  g <- gTrack(gr, height=7, xaxis.chronly = TRUE)
  expect_equal(g$legend, FALSE)
  expect_equal(g$height, 7)

})
