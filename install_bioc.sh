#!/bin/bash

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomicRanges"); install.packages("devtools", repos = "https://cran.rstudio.com/"); install_github("mskilab/gUtils"); devtools::install_github("jimhester/covr");'
