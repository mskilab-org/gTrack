#!/bin/bash

Rscript -e 'install.packages("devtools"); source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomicRanges"); install_github("mskilab/gUtils"); devtools::install_github("jimhester/covr");'
