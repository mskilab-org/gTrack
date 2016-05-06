Customizing a Small Data Set
===============================

In this vignette, examples of how to segment a data set such as a single GRanges object, how to specify the y-axis of a graph, how to color that same graph, how to add a color to each unique value will be shown. 

gr.tile(gr , w) - Divide GRanges into tiles of length "w"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. sourcecode:: r
    

    #The only interval in this GRanges object has a range of length 100, it'll be divided by 5 and thus, 20 tiles of length 5 will be returned.
    gr <- gr.tile(GRanges(1, IRanges(1,100)), w=5)



.. sourcecode:: r
    

    ## Plot tiles 
    plot(gTrack(gr))

.. figure:: figure/plot-tiles-1.png
    :alt: plot of chunk plot-tiles

    plot of chunk plot-tiles

gTrack(gr + n) - Extend each range by "n" base pairs 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

.. sourcecode:: r
    

    plot(gTrack(gr+5))

.. figure:: figure/plot-overlappingtiles-1.png
    :alt: plot of chunk plot-overlappingtiles

    plot of chunk plot-overlappingtiles

stack.gap - Specify degree of spacing(in y-direction) between tiles. 
~~~~~~~~~


.. sourcecode:: r
    

    gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3") ,
      c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) , end =
        c(2,4,6,8,10,12,14,16,18,20), names = head(letters,10)), GC=seq(1,10,length=10), name=seq(5,10,length=10))



.. sourcecode:: r
    

    plot(gTrack(gr))

.. figure:: figure/plot-gr-1.png
    :alt: plot of chunk plot-gr

    plot of chunk plot-gr


.. sourcecode:: r
    

    plot(gTrack(gr , stack.gap = 2))

.. figure:: figure/plot-stack.gap2-1.png
    :alt: plot of chunk plot-stack.gap2

    plot of chunk plot-stack.gap2


.. sourcecode:: r
    

    plot(gTrack(gr , stack.gap = 3))

.. figure:: figure/plot-stack.gap3-1.png
    :alt: plot of chunk plot-stack.gap3

    plot of chunk plot-stack.gap3

y.field - Specify y-axis of graph 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC'))

.. figure:: figure/plot-y.fieldGC-1.png
    :alt: plot of chunk plot-y.fieldGC

    plot of chunk plot-y.fieldGC

**gTrack(gr , bars = TRUE/FALSE)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = 'light blue'))

.. figure:: figure/plot-bars-1.png
    :alt: plot of chunk plot-bars

    plot of chunk plot-bars

**gTrack(gr , lines = TRUE/FALSE)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , lines = TRUE , col = 'purple'))

.. figure:: figure/plot-lines-1.png
    :alt: plot of chunk plot-lines

    plot of chunk plot-lines

**gTrack(gr , circles = TRUE/FALSE)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , circles = TRUE , col = 'magenta' , border = '60'))

.. figure:: figure/plot-circles-1.png
    :alt: plot of chunk plot-circles

    plot of chunk plot-circles

**colorfield**

**map values to colors! Legend is automatically added**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , colormaps = list(GC = c("1"="red" , "2" = "blue" , "3"="magenta", "4"="light blue" ,"5"="black" , "6"="green", "7"="brown" , "8"="pink", "9"="yellow", "10" = "orange")) ))

.. figure:: figure/plot-colorfield-1.png
    :alt: plot of chunk plot-colorfield

    plot of chunk plot-colorfield

**gr.colorfield**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC'))

.. figure:: figure/plot-gr.colorfield-1.png
    :alt: plot of chunk plot-gr.colorfield

    plot of chunk plot-gr.colorfield

**gr.labelfield**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC' , gr.labelfield = 'name'))

.. figure:: figure/plot-labelfield-1.png
    :alt: plot of chunk plot-labelfield

    plot of chunk plot-labelfield

**GRangesList**


.. sourcecode:: r
    

    # first, create GRanges object for chroms 1-3. Each chrom stores regions of exons
    chrom1 <- GRanges(seqnames=Rle(rep(1,5)) , ranges = IRanges(c(13214448,13377047,17190862,17284920,30741950) , end=c(13376489,17190004,17283075,30741656,30745210)))
    chrom2 <- GRanges(seqnames=Rle(rep(2,5)) , ranges = IRanges(c(34675467,34737163,50880025,50882016,51098931) , end = c(34737057,50879519,50880979,51089715,51099793)))
    chrom3 <- GRanges(seqnames=Rle(rep(3,5)) , ranges = IRanges(c(5883026,5888521,6651128,6655078,10251906) , end = c(5887648,6646543,6653332,10245198,10254797)))
    chroms <- GRangesList("chrom1" = chrom1 , "chrom2" = chrom2 , "chrom3" = chrom3)

