How to Graph Super Enhancers
============================

To illustrate gTrack's capability in exploring data sets such as ChiP-Seq data



.. sourcecode:: r
    

    ## The methods section of http://www.nature.com/ng/journal/v48/n2/full/ng.3470.html stated
    ## that GISTIC analyses were performed on the TCGA data available on the TCGA copy number
    ## portal which is created by the Broad Institute of MIT and Harvard.
    
    ## After finding version 3.0 of the SNP pipeline on 22 October 2014, clicking on the IGV
    ## session and an XML document
    (http://portals.broadinstitute.org/tcga/gisticIgv/session.xml?analysisId=21&tissueId=548&type=.xml)
    ## was returned which stored the web path to the *.seg.gz file. I downloaded that and found
    ## that it stored the log2 ratios (tumor coverage / normal coverage).
    
    wget http://www.broadinstitute.org/igvdata/tcga/tcgascape/141024_update/all_cancers.seg.gz
    wget http://www.broadinstitute.org/igvdata/tcga/tcgascape/141024_update/sample_info.txt
    
    gzip -d all_cancers.seg.gz


::

    ## Error: <text>:7:7: unexpected '/'
    ## 6: ## session and an XML document
    ## 7: (http:/
    ##          ^




.. sourcecode:: r
    

    ##############################                         ##############################
    
    ##############################    Starting Analysis    ##############################
    
    ##############################                         ##############################
    
    library(gTrack) ## main sauce. 
    library(gUtils) ## for dt2gr 
    
    ## Load coverage data into data.table. Very fast, thanks data.table.
    seg_data <- fread('../../inst/extdata/files/all_cancers.seg')


::

    ## Error in fread("../../inst/extdata/files/all_cancers.seg"): embedded nul in string: '4\x8a\xbe~i̮\xb1\xa9v\0054\x92\xde\r-\xd2R\xa0M\xaf8h{\xeb\030\0343\xab\001:v\x9b\a\xe7\037\xbc\x8f{\xe8\xf9`E-6\xdf`\xee\xe6C\xec;Cd\xa7\x80g\x83~\x82g\x91k\xa0q\xc9@\x98k4\033\xdeɘ\005\xef\034\x9c\xc3\xcft\xe9\xc5\xdfy\x92\xc6\xdf\035X\xc1\xd7\xcd?\xc0\xd7UF\xf0\x8d\t\177\xf9\xc6r\017\xb4\027-C{G\x8cG\acߡC\x81'ѡ\xc2;\xe8\b_\x80\x8e\xb89\xa1c\xc6\xfbб\xe5\xee\xe8\xf8\xf2Ft\xfc\xee}TY\x9b\x86\xce\037\xf1A\027\xa6䣋/\xbb\xa1*W\016\xf7L5\xc3=\017\xefƽϴ`\xa5\xe2i\xd8=\xd7\037{PK\xeco\x93\x83\003TRp@\xb0\001\x8e\t4\xc51\xfb\x8e\xe0\xd8;\xa58N\xee\006\xbe\023\xebA<*\xae\x91\xb1\xca\034\xf1z3\x87\x8c\xd3\xd0!\xde!]\x88\xf7mm2\xa98\x85\024\xf8ǒ\x82\xb2٤\xfe\xb4\aiP\xd2'\r\xa7\b\xb9\017\020\xf9z\xfc\033\xf9\xbe\xfc\017\xf9\xb9s9\035\xd5\U0001d3b6\x94\xa6\xa3\017Υc\xec\xe4鶺\xbbt\xbb\xf6\va\xc8\032\037aHk\x92`\xb0\xccX0\xb8\xa7&X6\034\025\xac\xa2\x8d\005\xab5\x89\xc2\xe4\xbe\017\x84ɓ\016K\xfc$\xfd:\026\x8f\004\xabl\xa6\033\xea\026\xd1}\ru\xf4\xc0\xf1\xfb\x


.. sourcecode:: r
    

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


