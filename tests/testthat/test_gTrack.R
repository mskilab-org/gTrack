library(gTrack)
library(gUtils)

library(testthat)

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
  
  
  list(
    coverage_gr = coverage_gr,
    plus_coverage_gr = plus_coverage_gr,
    minus_coverage_gr = minus_coverage_gr,
    coverage_gt = coverage_gt,
    plus_coverage_gt = plus_coverage_gt,
    minus_coverage_gt = minus_coverage_gt,
    coverage_gt_rd = coverage_gt_rd,
    gg = gg
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
  fp <- parse.gr("1:6850000-7050000")
  
  # Test unstranded intervals
  expect_error(plot(test_data$coverage_gt, fp), NA)
  expect_gt(create_plot_file({plot(test_data$coverage_gt, fp)}), 0)
  
  # Test stranded intervals with "+"
  expect_error(plot(test_data$plus_coverage_gt, fp), NA)
  expect_gt(create_plot_file({plot(test_data$plus_coverage_gt, fp)}), 0)
  
  
  # Test stranded intervals with "-"
  expect_error(plot(test_data$minus_coverage_gt, fp), NA)
  expect_gt(create_plot_file({plot(test_data$minus_coverage_gt, fp)}), 0)
  
})

test_that("gTrack scatter plot function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  # Test scatter plot
  expect_error(plot(test_data$coverage_gt_rd, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(test_data$coverage_gt_rd, fp + 1e5)}), 0)  
})

test_that("gTrack bar plot function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  # Test gTrack object
  coverage_bars_gt <- gTrack(test_data$coverage_gr, y.field = "cn", bars = TRUE, y0 = 0, y1 = 12)
  coverage_bars_gt_saved <- readRDS(system.file("extdata", "test_data", "coverage_bars_gt.rds", package = "gTrack"))
  expect_equal(coverage_bars_gt, coverage_bars_gt_saved)
  
  # Test bar plot
  expect_error(plot(coverage_bars_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(coverage_bars_gt, fp + 1e5)}), 0)  
  
})

test_that("gTrack line plot function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  # Test gTrack object
  coverage_lines_gt <- gTrack(test_data$coverage_gr, y.field = "cn", lines = TRUE, y0 = 0, y1 = 12)
  coverage_lines_gt_saved <- readRDS(system.file("extdata", "test_data", "coverage_lines_gt.rds", package = "gTrack"))
  expect_equal(coverage_lines_gt, coverage_lines_gt_saved)
  
  # Test line plot
  expect_error(plot(coverage_lines_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(coverage_lines_gt, fp + 1e5)}), 0)  
  
})

test_that("gTrack multiple plots function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  coverage_bars_gt <- gTrack(test_data$coverage_gr, y.field = "cn", bars = TRUE, y0 = 0, y1 = 12)
  coverage_lines_gt <- gTrack(test_data$coverage_gr, y.field = "cn", lines = TRUE, y0 = 0, y1 = 12)

  concatenated_gt <- c(test_data$coverage_gt_rd, coverage_bars_gt, coverage_lines_gt)
  
  # Test multiple plots
  expect_error(plot(concatenated_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, fp + 1e5)}), 0)  
  
})

test_that("gTrack unordered GRangesList (default) function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  # Test gTrack object
  reads <- readRDS(system.file("extdata", "ovcar.subgraph.reads.rds", package = "gTrack"))
  reads_gt <- gTrack(reads)
  reads_gt_saved <- readRDS(system.file("extdata", "test_data", "reads_gt.rds", package = "gTrack"))
  expect_equal(reads_gt, reads_gt_saved)
  
  concatenated_gt <- c(test_data$coverage_gt_rd, reads_gt)
  
  # Test unordered GRangesList (default) plot
  expect_error(plot(concatenated_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, fp + 1e5)}), 0)  
  
})

test_that("gTrack heatmap (mdata) function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  mdata_mat <- readRDS(system.file("extdata", "ovcar.subgraph.mdata.mat.rds", package = "gTrack"))
  mdata_gr <- readRDS(system.file("extdata", "ovcar.subgraph.mdata.gr.rds", package = "gTrack"))
  
  # Test gTrack object
  heatmap_gt <- gTrack(mdata_gr, mdata = mdata_mat, cmap.max = 10)
  heatmap_gt_saved <- readRDS(system.file("extdata", "test_data", "heatmap_gt.rds", package = "gTrack"))
  expect_equal(heatmap_gt, heatmap_gt_saved)
  
  concatenated_gt <- c(test_data$coverage_gt_rd, heatmap_gt)
  
  # Test heatmap (mdata) plot
  expect_error(plot(concatenated_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, fp + 1e5)}), 0)  
  
})

test_that("gTrack connections (edges) function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  edges_dat <- readRDS(system.file("extdata", "ovcar.subgraph.edges.dat.rds", package = "gTrack"))
  edges_gr <- readRDS(system.file("extdata", "ovcar.subgraph.edges.gr.rds", package = "gTrack"))

  # Test gTrack object
  edges_gt <- gTrack(edges_gr, edges = edges_dat)
  edges_gt_saved <- readRDS(system.file("extdata", "test_data", "edges_gt.rds", package = "gTrack"))
  expect_equal(edges_gt, edges_gt_saved)
  concatenated_gt <- c(test_data$coverage_gt_rd, edges_gt)
  
  # Test connections (edges) plot
  expect_error(plot(concatenated_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, fp + 1e5)}), 0)  
  
})

test_that("gGraph function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  concatenated_gt <- c(test_data$coverage_gt_rd, test_data$gg$gt)
  
  # Test gGraph plot
  expect_error(plot(concatenated_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, fp + 1e5)}), 0)  
  
})

test_that("gWalk function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  wks <- readRDS(system.file("extdata", "ovcar.subgraph.walks.rds", package = "gTrack"))
  concatenated_gt <- c(test_data$coverage_gt_rd, test_data$gg$gt, wks$gt)
  
  # Test gWalk plot
  expect_error(plot(concatenated_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, fp + 1e5)}), 0)  
  
})

test_that("gMatrix function works as expected", {
  test_data <- create_test_data()
  fp <- parse.gr("1:6043576-7172800")
  
  # Test gTrack object
  gm <- readRDS(system.file("extdata", "ovcar.subgraph.hic.rds", package = "gTrack"))
  gm_gt_saved <- readRDS(system.file("extdata", "test_data", "gm_gt.rds", package = "gTrack"))
  expect_equal(gm$gtrack(cmap.max = 1000), gm_gt_saved)
  
  concatenated_gt <- c(test_data$coverage_gt_rd, test_data$gg$gt, gm$gtrack(cmap.max = 1000))
  
  # Test gMatrix plot
  expect_error(plot(concatenated_gt, fp + 1e5), NA)
  expect_gt(create_plot_file({plot(concatenated_gt, fp + 1e5)}), 0)  
  
})
