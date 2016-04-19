Vignette Using CNV Data
=======================

**Operations on CNV data of one sample** 


.. sourcecode:: r
    

    getwd()


::

    ## [1] "/Users/knagdimov/Documents/gTrack/gTrack/rtdocs/vignetteCNV"


.. sourcecode:: r
    

    list.files()


::

    ## [1] "vignetteCNV.Rrst"  "vignetteCNV.Rrst~" "vignetteCNV.rst"


.. sourcecode:: r
    

    ##load data from TCGA
    
    tcgaData <- read.delim("gTrack/inst/extdata/BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg18.seg.txt")


::

    ## Warning in file(file, "rt"): cannot open file 'gTrack/inst/extdata/
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




.. sourcecode:: r
    

    # plot gTrack object
    
    plot(tcgagt)


::

    ## Error in plot(tcgagt): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'tcgagt' not found




.. sourcecode:: r
    

    # Changed default number of windows (from entire genome to the first five windows).
    # Windows argument requires a subset of a GRanges Object. Check documentation for more details.
    plot(tcgagt , window = tcgagr[1:5])


::

    ## Error in plot(tcgagt, window = tcgagr[1:5]): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'tcgagt' not found




.. sourcecode:: r
    

    # use amplifications/deletions as y-values
    tcgagt <- gTrack(tcgagr , y.field="Segment_Mean")


::

    ## Error in listify(data, GRanges): object 'tcgagr' not found


.. sourcecode:: r
    

    plot(tcgagt , windows = tcgagrr[1:5] , col = "red")


::

    ## Error in plot(tcgagt, windows = tcgagrr[1:5], col = "red"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'tcgagt' not found




.. sourcecode:: r
    

    # add a second sample to the graph
    # create gTrack object for sample
    tcgaData2 <- read.delim("inst/extdata/BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg19.seg.txt")


::

    ## Warning in file(file, "rt"): cannot open file 'inst/extdata/
    ## BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg19.seg.txt': No such
    ## file or directory



::

    ## Error in file(file, "rt"): cannot open the connection


.. sourcecode:: r
    

    tcgagr2 <- GRanges(tcgaData2)


::

    ## Error in newGRanges("GRanges", seqnames = seqnames, ranges = ranges, strand = strand, : object 'tcgaData2' not found


.. sourcecode:: r
    

    tcgagt2 <- gTrack(tcgagr2 , y.field="Segment_Mean")


::

    ## Error in listify(data, GRanges): object 'tcgagr2' not found




.. sourcecode:: r
    

    # plot the two samples
    plot(c(tcgagt2 , tcgagt), windows = tcgagr2[1:5] , col = "red")


::

    ## Error in plot(c(tcgagt2, tcgagt), windows = tcgagr2[1:5], col = "red"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'tcgagt2' not found




.. sourcecode:: r
    

    # physically separate gaps between tracks
    plot(c(tcgagt2 , tcgagt), windows = tcgagr2[1:5] , col = "red" , ygap = 5)


::

    ## Error in plot(c(tcgagt2, tcgagt), windows = tcgagr2[1:5], col = "red", : error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'tcgagt2' not found




.. sourcecode:: r
    

    # study of the CNVs in breast cancer
    fn = list.files("Level_3/")
    
    # create data.tables for each patient but, combine them into one HUGE data.table using rbindlist
    dt = rbindlist(lapply(fn , function(x) fread(x , colClasses = "character")[ , file:=x]))
    
    # certain arguments (window) of gTrack require numeric vectors. Thus, "character" vectors need
    # to be converted into "numeric" vectors.
    
    dt$Start = type.convert(dt$Start)


::

    ## Error in type.convert(dt$Start): the first argument must be of mode character


.. sourcecode:: r
    

    dt$End = type.convert(dt$End)


::

    ## Error in type.convert(dt$End): the first argument must be of mode character


.. sourcecode:: r
    

    # because we are graphing segment mean, that column also needs to be "numeric"
    dt$Segment_Mean = type.convert(dt$Segment_Mean)


::

    ## Error in type.convert(dt$Segment_Mean): the first argument must be of mode character


.. sourcecode:: r
    

    # convert data.table into GRanges object
    dtgr = GRanges(dt)


::

    ## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'Rle' for signature '"data.table", "missing"'


.. sourcecode:: r
    

    # wrap a gTrack object around it and plot
    dtgt <- gTrack(dtgr , y.field = "Segment_Mean")


::

    ## Error in listify(data, GRanges): object 'dtgr' not found




.. sourcecode:: r
    

    plot(dtgt , window = dtgr[1:5])


::

    ## Error in plot(dtgt, window = dtgr[1:5]): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'dtgt' not found




.. sourcecode:: r
    

    # show amplifications only (use gUtils operators!)
    dtgr = dtgr %Q% (Segment_Mean > 0)


::

    ## Error in dtgr %Q% (Segment_Mean > 0): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error: object 'dtgr' not found


.. sourcecode:: r
    

    dtgt <- gTrack(dtgr , y.field = "Segment_Mean")


::

    ## Error in listify(data, GRanges): object 'dtgr' not found




.. sourcecode:: r
    

    plot(dtgt , window = dtgr[1:5])


::

    ## Error in plot(dtgt, window = dtgr[1:5]): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'dtgt' not found




.. sourcecode:: r
    

    # show deletions only (again, use gUtils operators!)
    
    # recreate the original GRanges object
    dtgr = GRanges(dt)


::

    ## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'Rle' for signature '"data.table", "missing"'


.. sourcecode:: r
    

    # subset properly
    dtgr = dtgr %Q% (Segment_Mean < 0)


::

    ## Error in dtgr %Q% (Segment_Mean < 0): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error: object 'dtgr' not found


.. sourcecode:: r
    

    dtgt <- gTrack(dtgr , y.field = "Segment_Mean")


::

    ## Error in listify(data, GRanges): object 'dtgr' not found




.. sourcecode:: r
    

    plot(dtgt , window = dtgr[1:5])


::

    ## Error in plot(dtgt, window = dtgr[1:5]): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'dtgt' not found


