How to Graph Super Enhancers
============================

To illustrate gTrack's capability in exploring data sets such as ChiP-Seq data



.. sourcecode:: r
    

    ## The methods section of http://www.nature.com/ng/journal/v48/n2/full/ng.3470.html stated
    ## that GISTIC analyses were performed on the TCGA data available on the TCGA copy number
    ## portal which is created by the Broad Institute of MIT and Harvard.
    
    ## After finding version 3.0 of the SNP pipeline on 22 October 2014, clicking on the IGV
    ## session and an XML document (http://portals.broadinstitute.org/tcga/gisticIgv/session.xml?analysisId=21&tissueId=548&type=.xml)
    ## was returned which stored the web path to the *.seg.gz file. I downloaded that and found
    ## that it stored the log2 ratios (tumor coverage / normal coverage).
    
    ## wget http://www.broadinstitute.org/igvdata/tcga/tcgascape/141024_update/all_cancers.seg.gz
    ## wget http://www.broadinstitute.org/igvdata/tcga/tcgascape/141024_update/sample_info.txt
    
    ## gzip -d all_cancers.seg.gz



.. sourcecode:: r
    

    ##############################                         ##############################
    
    ##############################    Starting Analysis    ##############################
    
    ##############################                         ##############################
    
    library(gTrack) ## main sauce. 
    library(gUtils) ## for dt2gr 
    
    ## Load coverage data into data.table. Very fast, thanks data.table.
    ## seg_data <- fread('../../inst/extdata/files/all_cancers.seg')
    
    ## Coerce into GRanges from data.table because gTrack operates on GRanges.
    seg_ranges <- dt2gr(seg_data)


::

    ## Warning in dt2gr(seg_data): coercing to GRanges via non-standard columns


.. sourcecode:: r
    

    ## Subset to MYC enhancer amplification regions.
    seg_data_chrom8 <- seg_data[ Chromosome == 8]
    
    ## coerce into GRanges from data.table because gTrack operates on GRanges.
    seg_ranges_chrom8 <- dt2gr(seg_data_chrom8)


::

    ## Warning in dt2gr(seg_data_chrom8): coercing to GRanges via non-standard
    ## columns




.. sourcecode:: r
    

    ## first MYC(myc) (s)uper-(e)nhancer.
    myc_se <- parse.gr(c('8:129543949-129554294'))
    
    ## zoom into that region to view CNA.
    win <- myc_se

