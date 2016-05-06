Customizing a Small Data Set
===============================

In this vignette, examples of how to segment a data set such as a single GRanges object, how to specify the y-axis of a graph, how to color that same graph, how to add a color to each unique value will be shown. 

gr.tile(gr , w) - Divide GRanges into tiles of length "w"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. sourcecode:: r
    

    #The only interval in this GRanges object has a range of length 100, it'll be divided by 5 and thus, 20 tiles of length 5 will be returned.
    gr <- gr.tile(GRanges(1, IRanges(1,100)), w=5)


::

    ## Error in eval(expr, envir, enclos): could not find function "gr.tile"




.. sourcecode:: r
    

    ## Plot tiles 
    plot(gTrack(gr))


::

    ## Error in plot(gTrack(gr)): could not find function "gTrack"



gTrack(gr + n) - specify degree of overlap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

.. sourcecode:: r
    

    plot(gTrack(gr+7))


::

    ## Error in plot(gTrack(gr + 7)): could not find function "gTrack"



**stack.gap**

**vector or scaler numeric specifiying x gap between stacking non-numeric GRanges or GRangesLists items in track(s).**


.. sourcecode:: r
    

    gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3") ,
      c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) , end =
        c(2,4,6,8,10,12,14,16,18,20), names = head(letters,10)), GC=seq(1,10,length=10), name=seq(5,10,length=10))


::

    ## Error in eval(expr, envir, enclos): could not find function "GRanges"


.. sourcecode:: r
    

    print(gr)


::

    ## Error in print(gr): object 'gr' not found




.. sourcecode:: r
    

    plot(gTrack(gr))


::

    ## Error in plot(gTrack(gr)): could not find function "gTrack"




.. sourcecode:: r
    

    plot(gTrack(gr , stack.gap = 2))


::

    ## Error in plot(gTrack(gr, stack.gap = 2)): could not find function "gTrack"




.. sourcecode:: r
    

    plot(gTrack(gr , stack.gap = 3))


::

    ## Error in plot(gTrack(gr, stack.gap = 3)): could not find function "gTrack"



**gTrack(gr , y.field = 'GC')**

**vector or scalar numeric specifiying gap between tracks (add a dimension to the data)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC'))


::

    ## Error in plot(gTrack(gr, y.field = "GC")): could not find function "gTrack"



**gTrack(gr , bars = TRUE/FALSE)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = 'light blue'))


::

    ## Error in plot(gTrack(gr, y.field = "GC", bars = TRUE, col = "light blue")): could not find function "gTrack"



**gTrack(gr , lines = TRUE/FALSE)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , lines = TRUE , col = 'purple'))


::

    ## Error in plot(gTrack(gr, y.field = "GC", lines = TRUE, col = "purple")): could not find function "gTrack"



**gTrack(gr , circles = TRUE/FALSE)**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , circles = TRUE , col = 'magenta' , border = '60'))


::

    ## Error in plot(gTrack(gr, y.field = "GC", circles = TRUE, col = "magenta", : could not find function "gTrack"



**colorfield**

**map values to colors! Legend is automatically added**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , colormaps = list(GC = c("1"="red" , "2" = "blue" , "3"="magenta", "4"="light blue" ,"5"="black" , "6"="green", "7"="brown" , "8"="pink", "9"="yellow", "10" = "orange")) ))


::

    ## Error in plot(gTrack(gr, y.field = "GC", bars = TRUE, col = NA, colormaps = list(GC = c(`1` = "red", : could not find function "gTrack"



**gr.colorfield**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC'))


::

    ## Error in plot(gTrack(gr, y.field = "GC", bars = TRUE, col = NA, gr.colorfield = "GC")): could not find function "gTrack"



**gr.labelfield**


.. sourcecode:: r
    

    plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = 'GC' , gr.labelfield = 'name'))


::

    ## Error in plot(gTrack(gr, y.field = "GC", bars = TRUE, col = NA, gr.colorfield = "GC", : could not find function "gTrack"



**GRangesList**


.. sourcecode:: r
    

    # first, create GRanges object for chroms 1-3. Each chrom stores regions of exons
    chrom1 <- GRanges(seqnames=Rle(rep(1,5)) , ranges = IRanges(c(13214448,13377047,17190862,17284920,30741950) , end=c(13376489,17190004,17283075,30741656,30745210)))


::

    ## Error in eval(expr, envir, enclos): could not find function "GRanges"


.. sourcecode:: r
    

    chrom2 <- GRanges(seqnames=Rle(rep(2,5)) , ranges = IRanges(c(34675467,34737163,50880025,50882016,51098931) , end = c(34737057,50879519,50880979,51089715,51099793)))


::

    ## Error in eval(expr, envir, enclos): could not find function "GRanges"


.. sourcecode:: r
    

    chrom3 <- GRanges(seqnames=Rle(rep(3,5)) , ranges = IRanges(c(5883026,5888521,6651128,6655078,10251906) , end = c(5887648,6646543,6653332,10245198,10254797)))


::

    ## Error in eval(expr, envir, enclos): could not find function "GRanges"


.. sourcecode:: r
    

    chroms <- GRangesList("chrom1" = chrom1 , "chrom2" = chrom2 , "chrom3" = chrom3)


::

    ## Error in eval(expr, envir, enclos): could not find function "GRangesList"


