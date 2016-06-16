Customizing a Small Data Set
===============================

In this vignette, examples of how to segment a data set such as a single GRanges object, how to specify the y-axis of a graph, how to color that same graph, how to add a color to each unique value will be shown. 

gr.tile(gr , w) - Divide GRanges into tiles of length "w"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. sourcecode:: r
    

    ## DO NOT FORGET TO LOAD gUtils library.
    library(gUtils)
    
    #The only interval in this GRanges object has a range of length 100, it'll be divided by 5 and thus, 20 tiles of length 5 will be returned.
    gr <- gr.tile(GRanges(1, IRanges(1,100)), w=5)



.. sourcecode:: r
    

    ## Plot tiles 
    plot(gTrack(gr))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-tiles-1.png
    :alt: plot of chunk plot-tiles

    plot of chunk plot-tiles

gTrack(gr + n) - Extend each range by "n" base pairs 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

.. sourcecode:: r
    

    plot(gTrack(gr+5))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated



::

    ## Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 1 out-of-bound range located on sequence
    ##   1. Note that only ranges located on a non-circular sequence whose
    ##   length is not NA can be considered out-of-bound (use seqlengths()
    ##   and isCircular() to get the lengths and circularity flags of the
    ##   underlying sequences). You can use trim() to trim these ranges.
    ##   See ?`trim,GenomicRanges-method` for more information.


.. figure:: figure/plot-overlappingtiles-1.png
    :alt: plot of chunk plot-overlappingtiles

    plot of chunk plot-overlappingtiles

stack.gap - Specify degree of spacing(in x-direction) between **ADJACENT** tiles. 
~~~~~~~~~


.. sourcecode:: r
    

    gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3") ,
      c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) , end =
        c(2,4,6,8,10,12,14,16,18,20), names = head(letters,10)), GC=seq(1,10,length=10), name=seq(5,10,length=10))



.. sourcecode:: r
    

    plot(gTrack(gr))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-gr-1.png
    :alt: plot of chunk plot-gr

    plot of chunk plot-gr


.. sourcecode:: r
    

    plot(gTrack(gr , stack.gap = 2))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-stack.gap2-1.png
    :alt: plot of chunk plot-stack.gap2

    plot of chunk plot-stack.gap2


.. sourcecode:: r
    

    plot(gTrack(gr , stack.gap = 3))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-stack.gap3-1.png
    :alt: plot of chunk plot-stack.gap3

    plot of chunk plot-stack.gap3

y.field - Specify y-axis of graph 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC'))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-y.fieldGC-1.png
    :alt: plot of chunk plot-y.fieldGC

    plot of chunk plot-y.fieldGC

bars - Plot data points as vertical bars 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**gTrack(gr , bars = TRUE/FALSE)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = 'light blue'))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-bars-1.png
    :alt: plot of chunk plot-bars

    plot of chunk plot-bars

lines - Plot data points as lines.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gTrack(gr , lines = TRUE/FALSE)


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , lines = TRUE , col = 'purple'))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-lines-1.png
    :alt: plot of chunk plot-lines

    plot of chunk plot-lines

circles - Plot data points as circles. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gTrack(gr , circles = TRUE/FALSE)


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , circles = TRUE , col = 'magenta' , border = '60'))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-circles-1.png
    :alt: plot of chunk plot-circles

    plot of chunk plot-circles

colorfield - Specify mapping of colors to values.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , colormaps = list(GC = c("1"="red" , "2" = "blue" , "3"="magenta", "4"="light blue" ,"5"="black" , "6"="green", "7"="brown" , "8"="pink", "9"="yellow", "10" = "orange")) ))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-colorfield-1.png
    :alt: plot of chunk plot-colorfield

    plot of chunk plot-colorfield

gr.colorfield - Automatically specify mapping of colors to values. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC'))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-gr.colorfield-1.png
    :alt: plot of chunk plot-gr.colorfield

    plot of chunk plot-gr.colorfield

gr.labelfield - Plot values for each data point. 
~~~~~~~~~~~~~~


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC' , gr.labelfield = 'name'))


::

    ## Warning in `[<-`(`*tmp*`, null.ix, value = list(<S4 object of class
    ## structure("GRanges", package = "GenomicRanges")>)): implicit list embedding
    ## of S4 objects is deprecated


.. figure:: figure/plot-labelfield-1.png
    :alt: plot of chunk plot-labelfield

    plot of chunk plot-labelfield
