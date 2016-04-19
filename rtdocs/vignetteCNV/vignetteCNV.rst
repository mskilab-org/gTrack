Vignette Using CNV Data
=======================

**Operations on CNV data of one sample** 


.. sourcecode:: r
    

    setwd("~/Documents/gTrack/gTrack/")
    ok = list.files()
    ##load data from TCGA
    
    tcgaData <- read.delim("inst/extdata/BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg18.seg.txt")
    
    ## convert data.frame to GRanges object
    
    tcgagr <- GRanges(tcgaData)
    
    # wrap gTrack around TCGA GRanges object
    
    tcgagt <- gTrack(tcgagr)



.. sourcecode:: r
    

    # plot gTrack object
    
    plot(tcgagt)

.. figure:: figure/plot-entireData-1.png
    :alt: plot of chunk plot-entireData

    plot of chunk plot-entireData


.. sourcecode:: r
    

    # Changed default number of windows (from entire genome to the first five windows).
    # Windows argument requires a subset of a GRanges Object. Check documentation for more details.
    plot(tcgagt , window = tcgagr[1:5])

.. figure:: figure/plot-subsetofData-1.png
    :alt: plot of chunk plot-subsetofData

    plot of chunk plot-subsetofData


.. sourcecode:: r
    

    # use amplifications/deletions as y-values
    tcgagt <- gTrack(tcgagr , y.field="Segment_Mean")
    plot(tcgagt , windows = tcgagr[1:5] , col = "red")

.. figure:: figure/plot-amplifcations-1.png
    :alt: plot of chunk plot-amplifcations

    plot of chunk plot-amplifcations


.. sourcecode:: r
    

    # add a second sample to the graph
    # create gTrack object for sample
    
    setwd("~/Documents/gTrack/gTrack/")
    tcgaData2 <- read.delim("inst/extdata/BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg19.seg.txt")
    tcgagr2 <- GRanges(tcgaData2)
    tcgagt2 <- gTrack(tcgagr2 , y.field="Segment_Mean")



.. sourcecode:: r
    

    # plot the two samples
    plot(c(tcgagt2 , tcgagt), windows = tcgagr2[1:5] , col = "red")


::

    ## Error in plot(c(tcgagt2, tcgagt), windows = tcgagr2[1:5], col = "red"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error in do.call("rrbind", lapply(args, formatting)) : 
    ##   could not find function "rrbind"




.. sourcecode:: r
    

    # physically separate gaps between tracks
    plot(c(tcgagt2 , tcgagt), windows = tcgagr2[1:5] , col = "red" , ygap = 5)


::

    ## Error in plot(c(tcgagt2, tcgagt), windows = tcgagr2[1:5], col = "red", : error in evaluating the argument 'x' in selecting a method for function 'plot': Error in do.call("rrbind", lapply(args, formatting)) : 
    ##   could not find function "rrbind"




.. sourcecode:: r
    

    # study of the CNVs in breast cancer
    fn = list.files("~/Documents/gTrack/gTrack/inst/extdata/Level_3/")
    
    # create data.tables for each patient but, combine them into one HUGE data.table using rbindlist
    dt = rbindlist(lapply(fn , function(x) fread(x , colClasses = "character")[ , file:=x]))


::

    ## Error in eval(expr, envir, enclos): could not find function "rbindlist"


.. sourcecode:: r
    

    # certain arguments (window) of gTrack require numeric vectors. Thus, "character" vectors need
    # to be converted into "numeric" vectors.
    
    dt$Start = type.convert(dt$Start)


::

    ## Error in dt$Start: object of type 'closure' is not subsettable


.. sourcecode:: r
    

    dt$End = type.convert(dt$End)


::

    ## Error in dt$End: object of type 'closure' is not subsettable


.. sourcecode:: r
    

    # because we are graphing segment mean, that column also needs to be "numeric"
    dt$Segment_Mean = type.convert(dt$Segment_Mean)


::

    ## Error in dt$Segment_Mean: object of type 'closure' is not subsettable


.. sourcecode:: r
    

    # convert data.table into GRanges object
    dtgr = GRanges(dt)


::

    ## Error in as(seqnames, class): no method or default for coercing "function" to "GRanges"


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

    ## Error in eval(expr, envir, enclos): could not find function "%Q%"


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

    ## Error in as(seqnames, class): no method or default for coercing "function" to "GRanges"


.. sourcecode:: r
    

    # subset properly
    dtgr = dtgr %Q% (Segment_Mean < 0)


::

    ## Error in eval(expr, envir, enclos): could not find function "%Q%"


.. sourcecode:: r
    

    dtgt <- gTrack(dtgr , y.field = "Segment_Mean")


::

    ## Error in listify(data, GRanges): object 'dtgr' not found




.. sourcecode:: r
    

    plot(dtgt , window = dtgr[1:5])


::

    ## Error in plot(dtgt, window = dtgr[1:5]): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'dtgt' not found


