How to Plot CNV Data
=======================

gTrack is capable of graphing CNV data, quickly. In this vignette, the y.field parameter will be used.

y.field Parameter
~~~~~~~~~~~~~~~~~

.. sourcecode:: r
    # load breast cancer CNV data from tcgaData.
    setwd("~/gTrack/inst/extdata/Level_3")
    fn = list.files()
    
    # create data.tables for each patient but, combine them into one data.table.
    dt = rbindlist(lapply(fn , function(x) fread(x , colClasses = "character")[ , file:=x]))
    
    # certain arguments (window) of gTrack require numeric vectors. Thus, "character" vectors need
    # to be converted into "numeric" vectors.
    
    dt$Start = type.convert(dt$Start)
    dt$End = type.convert(dt$End)
    
    # because we are graphing segment mean, column also needs to be "numeric".
    dt$Segment_Mean = type.convert(dt$Segment_Mean)
    
    # convert data.table into GRanges object.
    dtgr = GRanges(dt)
    
    # wrap a gTrack object around it and plot
    dtgt <- gTrack(dtgr , y.field = "Segment_Mean")

.. sourcecode:: r
    
    plot(dtgt , window = dtgr[1:5])

.. figure:: figure/plot-allSamples-1.png
    :alt: plot of chunk plot-allSamples

    plot of chunk plot-allSamples


.. sourcecode:: r
    

    # show amplifications only (use gUtils operators!)
    dtgr = dtgr %Q% (Segment_Mean > 0)
    dtgt <- gTrack(dtgr , y.field = "Segment_Mean")



.. sourcecode:: r
    

    plot(dtgt , window = dtgr[1:5])

.. figure:: figure/plot-amplificationsAll-1.png
    :alt: plot of chunk plot-amplificationsAll

    plot of chunk plot-amplificationsAll


.. sourcecode:: r
    

    # show deletions only (again, use gUtils operators!)
    
    # recreate the original GRanges object
    dtgr = GRanges(dt)
    # subset properly
    dtgr = dtgr %Q% (Segment_Mean < 0)
    dtgt <- gTrack(dtgr , y.field = "Segment_Mean")



.. sourcecode:: r
    

    plot(dtgt , window = dtgr[1:5])

.. figure:: figure/plot-deletionsAll-1.png
    :alt: plot of chunk plot-deletionsAll

    plot of chunk plot-deletionsAll
