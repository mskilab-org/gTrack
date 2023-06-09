
#library(covr)
#report()
#roject_path = "~/projects/gTrack" ## this should be the path to your gTrack clone
#devtools::load_all(project_path)

#devtools::install_github('mskilab/skidb')
library(gTrack)
library(gUtils)
library(testthat)
library(bamUtils)


context("gTrack")

create_plot_file <- function(function_call) {
  temp_file <- tempfile(fileext = ".png")
  png(temp_file)
  function_call
  dev.off()
  file.size(temp_file)
}

create_test_data <- function() {
  coverage_gr_path <- system.file("extdata", "ovcar.subgraph.coverage.rds", package = "gTrack")
  coverage_gr <- readRDS(coverage_gr_path)
  
  plus_coverage_gr <- coverage_gr
  strand(plus_coverage_gr) <- "+"
  
  minus_coverage_gr <- coverage_gr
  strand(minus_coverage_gr) <- "-"
  
  # Test gTrack objects
  coverage_gt <- gTrack(coverage_gr)
  coverage_gt_saved <- readRDS(system.file("extdata", "test_data", "coverage_gt.rds", package = "gTrack"))
  expect_equal(coverage_gt, coverage_gt_saved)
  
  plus_coverage_gt <- gTrack(plus_coverage_gr)
  minus_coverage_gt <- gTrack(minus_coverage_gr)
  coverage_gt_rd <- gTrack(coverage_gr, y.field = "cn", circles = TRUE, lwd.border = 0.2, y0 = 0, y1 = 12)

  gg <- readRDS(system.file("extdata", "ovcar.subgraph.rds", package = "gTrack"))
  reads <- readRDS(system.file("extdata", "ovcar.subgraph.reads.rds", package = "gTrack"))
  
  fp <- parse.gr("1:6850000-7050000", seqlengths=coverage_gr$seqinfo)
  
  list(
    coverage_gr = coverage_gr,
    plus_coverage_gr = plus_coverage_gr,
    minus_coverage_gr = minus_coverage_gr,
    coverage_gt = coverage_gt,
    plus_coverage_gt = plus_coverage_gt,
    minus_coverage_gt = minus_coverage_gt,
    coverage_gt_rd = coverage_gt_rd,
    gg = gg,
    fp = fp,
    reads = reads
  )
}

test_that("GRanges and gTrack constructors work as expected", {
  test_data <- create_test_data()
  
  # Check that the GRanges and gTrack objects have the correct class
  expect_is(test_data$coverage_gr, "GRanges")
  expect_is(test_data$coverage_gt, "gTrack")
  expect_is(test_data$plus_coverage_gt, "gTrack")
  expect_is(test_data$minus_coverage_gt, "gTrack")
  
  # Test gTrack objects
  coverage_gt_saved <- readRDS(system.file("extdata", "test_data", "coverage_gt.rds", package = "gTrack"))
  expect_equal(test_data$coverage_gt, coverage_gt_saved)
  
  plus_coverage_gt_saved <- readRDS(system.file("extdata", "test_data", "plus_coverage_gt.rds", package = "gTrack"))
  expect_equal(test_data$plus_coverage_gt, plus_coverage_gt_saved)
  
  minus_coverage_gt_saved <- readRDS(system.file("extdata", "test_data", "minus_coverage_gt.rds", package = "gTrack"))
  expect_equal(test_data$minus_coverage_gt, minus_coverage_gt_saved)
  
  coverage_gt_saved <- readRDS(system.file("extdata", "test_data", "coverage_gt.rds", package = "gTrack"))
  expect_equal(test_data$coverage_gt, coverage_gt_saved)
  
  # Check that the strand information is correctly set
  expect_true(all(as.character(strand(test_data$plus_coverage_gr)) == "+"))
  expect_true(all(as.character(strand(test_data$minus_coverage_gr)) == "-"))
  expect_true(all(as.character(strand(test_data$coverage_gr)) == "*"))
})

test_that("gTrack plot function works as expected", {
  test_data <- create_test_data()
  
  # Test unstranded intervals
  expect_error(plot(test_data$coverage_gt, test_data$fp), NA)
  expect_gt(create_plot_file({plot(test_data$coverage_gt, test_data$fp)}), 0)
  
  # Test stranded intervals with "+"
  expect_error(plot(test_data$plus_coverage_gt, test_data$fp), NA)
  expect_gt(create_plot_file({plot(test_data$plus_coverage_gt, test_data$fp)}), 0)
  
  
  # Test stranded intervals with "-"
  expect_error(plot(test_data$minus_coverage_gt, test_data$fp), NA)
  expect_gt(create_plot_file({plot(test_data$minus_coverage_gt, test_data$fp)}), 0)
  
})

test_that("gTrack scatter plot function works as expected", {
  test_data <- create_test_data()
  
  # Test scatter plot
  expect_error(plot(test_data$coverage_gt_rd, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(test_data$coverage_gt_rd, test_data$fp + 1e5)}), 0)  
})

test_that("gTrack bar plot function works as expected", {
  test_data <- create_test_data()
  
  # Test gTrack object
  coverage_bars_gt <- gTrack(test_data$coverage_gr, y.field = "cn", bars = TRUE, y0 = 0, y1 = 12)
  coverage_bars_gt_saved <- readRDS(system.file("extdata", "test_data", "coverage_bars_gt.rds", package = "gTrack"))
  expect_equal(coverage_bars_gt, coverage_bars_gt_saved)
  
  # Test bar plot
  expect_error(plot(coverage_bars_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(coverage_bars_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("gTrack line plot function works as expected", {
  test_data <- create_test_data()
  
  # Test gTrack object
  coverage_lines_gt <- gTrack(test_data$coverage_gr, y.field = "cn", lines = TRUE, y0 = 0, y1 = 12)
  coverage_lines_gt_saved <- readRDS(system.file("extdata", "test_data", "coverage_lines_gt.rds", package = "gTrack"))
  expect_equal(coverage_lines_gt, coverage_lines_gt_saved)
  
  # Test line plot
  expect_error(plot(coverage_lines_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(coverage_lines_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("gTrack multiple plots function works as expected", {
  test_data <- create_test_data()
  
  coverage_bars_gt <- gTrack(test_data$coverage_gr, y.field = "cn", bars = TRUE, y0 = 0, y1 = 12)
  coverage_lines_gt <- gTrack(test_data$coverage_gr, y.field = "cn", lines = TRUE, y0 = 0, y1 = 12)

  concatenated_gt <- c(test_data$coverage_gt_rd, coverage_bars_gt, coverage_lines_gt)
  
  # Test multiple plots
  expect_error(plot(concatenated_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("gTrack unordered GRangesList (default) function works as expected", {
  test_data <- create_test_data()
  
  # Test gTrack object
  reads_gt <- gTrack(test_data$reads)
  reads_gt_saved <- readRDS(system.file("extdata", "test_data", "reads_gt.rds", package = "gTrack"))
  expect_equal(reads_gt, reads_gt_saved)
  
  concatenated_gt <- c(test_data$coverage_gt_rd, reads_gt)
  
  # Test unordered GRangesList (default) plot
  expect_error(plot(concatenated_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that('mdata() function', {
    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    heatMap = matrix(runif(length(gr)^2), nrow = 3, ncol = 3)
    gTrack_heatMap = gTrack(gr, mdata = heatMap)
    gTrack_heatMap_matrices = mdata(gTrack_heatMap)

    expect_is(gTrack_heatMap_matrices[[1]], "matrix")
    expect_equal(mdata(gTrack_heatMap, GRanges()), NULL)
    expect_is(mdata(gTrack_heatMap, GRanges(1, IRanges(1,100))), "matrix")
})

test_that("gTrack heatmap (mdata) function works as expected", {
  test_data <- create_test_data()
  
  mdata_mat <- readRDS(system.file("extdata", "ovcar.subgraph.mdata.mat.rds", package = "gTrack"))
  mdata_gr <- readRDS(system.file("extdata", "ovcar.subgraph.mdata.gr.rds", package = "gTrack"))
  
  # Test gTrack object
  heatmap_gt <- gTrack(mdata_gr, mdata = mdata_mat, cmap.max = 10)
  heatmap_gt_saved <- readRDS(system.file("extdata", "test_data", "heatmap_gt.rds", package = "gTrack"))
  expect_equal(heatmap_gt, heatmap_gt_saved)

  # Test mdata
  mdata_saved <- matrix(data = 19443, nrow = 1, ncol = 1,
                        dimnames = list(c("1:6043576-6053575"), c("1:6043576-6053575")))
  expect_equal(mdata(heatmap_gt, GRanges()), NULL)
  expect_equal(mdata(heatmap_gt, GRanges(1, IRanges(6043576,6045576))), mdata_saved)

  concatenated_gt <- c(test_data$coverage_gt_rd, heatmap_gt)
  
  # Test heatmap (mdata) plot
  expect_error(plot(concatenated_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("gTrack connections (edges) function works as expected", {
  test_data <- create_test_data()
  
  edges_dat <- readRDS(system.file("extdata", "ovcar.subgraph.edges.dat.rds", package = "gTrack"))
  edges_gr <- readRDS(system.file("extdata", "ovcar.subgraph.edges.gr.rds", package = "gTrack"))

  # Test gTrack object
  edges_gt <- gTrack(edges_gr, edges = edges_dat)
  edges_gt_saved <- readRDS(system.file("extdata", "test_data", "edges_gt.rds", package = "gTrack"))
  expect_equal(edges_gt, edges_gt_saved)
  concatenated_gt <- c(test_data$coverage_gt_rd, edges_gt)
  
  # Test connections (edges) plot
  expect_error(plot(concatenated_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("gGraph function works as expected", {
  test_data <- create_test_data()
  
  concatenated_gt <- c(test_data$coverage_gt_rd, test_data$gg$gt)
  
  # Test gGraph plot
  expect_error(plot(concatenated_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("gWalk function works as expected", {
  test_data <- create_test_data()
  
  wks <- readRDS(system.file("extdata", "ovcar.subgraph.walks.rds", package = "gTrack"))
  concatenated_gt <- c(test_data$coverage_gt_rd, test_data$gg$gt, wks$gt)
  
  # Test gWalk plot
  expect_error(plot(concatenated_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("gMatrix function works as expected", {
  test_data <- create_test_data()
  
  # Test gTrack object
  gm <- readRDS(system.file("extdata", "ovcar.subgraph.hic.rds", package = "gTrack"))
  gm_gt_saved <- readRDS(system.file("extdata", "test_data", "gm_gt.rds", package = "gTrack"))
  expect_equal(gm$gtrack(cmap.max = 1000), gm_gt_saved)
  
  concatenated_gt <- c(test_data$coverage_gt_rd, test_data$gg$gt, gm$gtrack(cmap.max = 1000))
  
  # Test gMatrix plot
  expect_error(plot(concatenated_gt, test_data$fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, test_data$fp + 1e5)}), 0)  
  
})

test_that("plot function handles links parameter correctly", {
  test_data <- create_test_data()

  # Create GRangesList corresponding to ALT edges
  grl = test_data$gg$junctions[type == "ALT"]$grl

  # Test the links parameter
  expect_error(plot(test_data$gg$gt, test_data$fp + 1e5, links = grl), NA)
  expect_gt(create_plot_file({plot(test_data$gg$gt, test_data$fp + 1e5, links = grl)}), 0)  
})

test_that("karyogram method works as expected", {
  fp <- parse.gr("1:1-200000000")
  karyogram_gt_hg18 = karyogram(hg19 = FALSE, bands = TRUE)
  karyogram_gt_hg19 = karyogram(hg19 = TRUE, bands = TRUE)
  karyogram_gt_no_bands = karyogram(hg19 = TRUE, bands = FALSE, arms = TRUE)
  karyogram_gt_no_arms = karyogram(hg19 = TRUE, bands = FALSE, arms = FALSE)
  
  # Test karyogram plot
  expect_error(plot(karyogram_gt_hg18, fp), NA)
  expect_gt(create_plot_file({plot(karyogram_gt_hg18, fp)}), 0)  
})

test_that("gencode constructor works as expected", {
  test_data <- create_test_data()

  fp <- parse.gr('1:6000000-6000100', seqlengths=test_data$coverage_gr$seqinfo$seqlengths)
  gencode_gt <- track.gencode()
  expect_error(plot(gencode_gt, fp + 1e4), NA)
  
  # Test uncached
  Sys.setenv(GENCODE_DIR = system.file('extdata', 'test_data', package = 'gTrack'))
  gencode_gt_uncached <- track.gencode(cached=FALSE)
  expect_error(plot(gencode_gt_uncached, fp + 1e4), NA)
  
  expect_gt(create_plot_file({plot(gencode_gt, fp + 1e4)}), 0)  
})

test_that('reduce method works as expected', {
    gt = gTrack(GRanges(1, IRanges(1,100)))
    expect_identical(reduce(gt), GRanges(1, IRanges(1, 100)))
})

test_that('show setMethod', {
    gt = gTrack(GRanges(1, IRanges(1,100)))
    yfield_shown_gt = show(gt)[[1]]
    expect_equal(yfield_shown_gt, NA)
})

test_that('brewer.master method works as expected', {
  brewer_colors <- brewer.master(3)
  expected_value <- c("#7FC97F", "#BEAED4", "#FDC086")
  expect_equal(brewer_colors, expected_value)
})

test_that('draw.ranges label route works as expected', {
  Sys.setenv(DEFAULT_GENOME = "")
  test_data <- create_test_data()
  plot(test_data$coverage_gt)
  
  coverage_gr_copy <- test_data$coverage_gr
  coverage_gr_copy$label <- 'test'
  
  expect_error(gTrack:::draw.ranges(coverage_gr_copy, angle=0), regexp = NA)
})

# test_that('draw.grl draw.var route works as expected', {
#   project_path = "~/projects/gTrack" ## this should be the path to your gTrack clone
  # devtools::load_all(project_path)
#   test_data <- create_test_data()
#   plot(test_data$coverage_gt, draw.var=TRUE)
#   
#   expect_error(gTrack:::draw.grl(test_data$reads, var=test_data$reads, draw.var=TRUE, windows=GRanges(1, 10), new.plot=FALSE, new.axis=FALSE), regexp = NA)
# })

test_that('col.scale internal function works as expected', {
  col_scale_ref = c("#000000", "#000000", "#000000")
  col_scale <- gTrack:::col.scale(c(1, 2, 3))
  expect_equal(col_scale_ref, col_scale)

  col_scale_inverted_ref = c("#FFFFFF", "#FFFFFF", "#FFFFFF")
  col_scale_inverted <- gTrack:::col.scale(c(1, 2, 3), invert=TRUE)
  expect_equal(col_scale_inverted_ref, col_scale_inverted)
})

test_that('clipping works as expected', {
  test_data <- create_test_data()
  plot(test_data$coverage_gt)
  
  expect_error(gTrack:::draw.ranges(test_data$coverage_gr, angle=0, clip=GRanges(1, IRanges(6043576,6044576))), regexp = NA)
})
