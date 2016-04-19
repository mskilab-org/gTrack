Vignette Using CNV Data
=======================

** Operations on CNV data of one sample ** 


.. sourcecode:: r
    

    ##load data from TCGA
    
    tcgaData <- read.delim("inst/extdata/BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg18.seg.txt")


::

    ## Warning in file(file, "rt"): cannot open file 'inst/extdata/
    ## BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg18.seg.txt': No such
    ## file or directory



::

    ## Error in file(file, "rt"): cannot open the connection


.. sourcecode:: r
    

    ## convert data.frame to GRanges object
    
    tcgagr <- GRanges(tcgaData)


::

    ## Error in newGRanges("GRanges", seqnames = seqnames, ranges = ranges, strand = strand, : object 'tcgaData' not found


.. sourcecode:: r
    

    # wrap gTrack around TCGA GRanges object
    
    tcgagt <- gTrack(tcgagr)


::

    ## Error in listify(data, GRanges): object 'tcgagr' not found


