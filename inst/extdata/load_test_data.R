library(gTrack)
library(gUtils)

create_test_data <- function() {
  coverage_gr_path <- system.file("extdata", "ovcar.subgraph.coverage.rds", package = "gTrack")
  coverage_gr <- readRDS(coverage_gr_path)
  
  plus_coverage_gr <- coverage_gr
  strand(plus_coverage_gr) <- "+"
  
  minus_coverage_gr <- coverage_gr
  strand(minus_coverage_gr) <- "-"
  
  coverage_gt <- gTrack(coverage_gr)
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

test_data <- create_test_data()
fp <- parse.gr("1:6043576-7172800")

objects_to_save <- list(
  coverage_gt = test_data$coverage_gt,
  plus_coverage_gt = test_data$plus_coverage_gt,
  minus_coverage_gt = test_data$minus_coverage_gt,
  coverage_gt_rd = test_data$coverage_gt_rd,
  coverage_bars_gt = gTrack(test_data$coverage_gr, y.field = "cn", bars = TRUE, y0 = 0, y1 = 12),
  coverage_lines_gt = gTrack(test_data$coverage_gr, y.field = "cn", lines = TRUE, y0 = 0, y1 = 12),
  reads_gt = gTrack(readRDS(system.file("extdata", "ovcar.subgraph.reads.rds", package = "gTrack"))),
  heatmap_gt = gTrack(
    readRDS(system.file("extdata", "ovcar.subgraph.mdata.gr.rds", package = "gTrack")),
    mdata = readRDS(system.file("extdata", "ovcar.subgraph.mdata.mat.rds", package = "gTrack")),
    cmap.max = 10
  ),
  edges_gt = gTrack(
    readRDS(system.file("extdata", "ovcar.subgraph.edges.gr.rds", package = "gTrack")),
    edges = readRDS(system.file("extdata", "ovcar.subgraph.edges.dat.rds", package = "gTrack"))
  ),
  gm_gt = readRDS(system.file("extdata", "ovcar.subgraph.hic.rds", package = "gTrack"))$gtrack(cmap.max = 1000)
)

# Set dir_path to the extdata folder of your local gTrack repository
local_clone_path <- "~/projects/gTrack"
dir_path <- file.path(local_clone_path, "inst", "extdata", "test_data")

# Check if the directory exists and create one if it doesn't
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

# Save objects to .rds files
for (object_name in names(objects_to_save)) {
  saveRDS(objects_to_save[[object_name]], file = file.path(dir_path, paste0(object_name, ".rds")))
}
