#!/bin/bash

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomicRanges"); install.packages("devtools"); library(devtools); devtools::install_github("mskilab/gUtils")'
