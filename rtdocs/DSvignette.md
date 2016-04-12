---
output: html_document
---
## gTrack Tutorial On An Artifically Small Data Set 



# Overlapping Segments 
### gr.tile(gr , w)
##### tile ranges across GRanges 
  


::

      ## Given a GRanges object, separate the interval by "w". 
      ## In this example, the output will be 20 equal sized tiles
    
      gr <- gr.tile(GRanges(1, IRanges(1,100)), w=5)
        
      ## Plot how the tiles would like 
      plot(gTrack(gr))


.. figure:: figure/unnamed-chunk-2-1.png
    :alt: plot of chunk unnamed-chunk-2

    plot of chunk unnamed-chunk-2

### Overlapping Tiles - gTrack(gr + n)
##### indicate the degree of overlap


::

      ## "1" in this example adds one more title to each interval 
      plot(gTrack(gr+1))


.. figure:: figure/unnamed-chunk-3-1.png
    :alt: plot of chunk unnamed-chunk-3

    plot of chunk unnamed-chunk-3

### gTrack(gr , stack.gap = n) 
##### vector or scaler numeric specifiying x gap between stacking non-numeric GRanges or GRangesLists items in track(s). 


::

      # first, create GRanges object.
      gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3") ,
      c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) , end =
      c(2,4,6,8,10,12,14,16,18,20), names = head(letters,10)), GC=seq(1,10,length=10), name=seq(5,10,length=10))
    
      print(gr)



::

    ## GRanges object with 10 ranges and 2 metadata columns:
    ##     seqnames    ranges strand |        GC             name
    ##        <Rle> <IRanges>  <Rle> | <numeric>        <numeric>
    ##   a     chr1  [ 1,  2]      * |         1                5
    ##   b     chr2  [ 3,  4]      * |         2 5.55555555555556
    ##   c     chr2  [ 5,  6]      * |         3 6.11111111111111
    ##   d     chr2  [ 7,  8]      * |         4 6.66666666666667
    ##   e     chr1  [ 9, 10]      * |         5 7.22222222222222
    ##   f     chr1  [11, 12]      * |         6 7.77777777777778
    ##   g     chr3  [13, 14]      * |         7 8.33333333333333
    ##   h     chr3  [15, 16]      * |         8 8.88888888888889
    ##   i     chr3  [17, 18]      * |         9 9.44444444444444
    ##   j     chr3  [19, 20]      * |        10               10
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome; no seqlengths



::

      # instead of plotting GRanges in a linear fashion, like this example 
      plot(gTrack(gr))



::

    ## Note: method with signature 'vectorORfactor#RangesNSBS' chosen for function 'extractROWS',
    ##  target signature 'matrix#RangesNSBS'.
    ##  "matrix#ANY" would also be valid


.. figure:: figure/unnamed-chunk-4-1.png
    :alt: plot of chunk unnamed-chunk-4

    plot of chunk unnamed-chunk-4

::

      # use the stack.gap parameter! 1by default but, if changed to any numeric(n), n levels will be created where 
      # the data will be plotted. 
      
      # 2 levels 
      plot(gTrack(gr , stack.gap=2))


.. figure:: figure/unnamed-chunk-4-2.png
    :alt: plot of chunk unnamed-chunk-4

    plot of chunk unnamed-chunk-4

::

      # 3 levels
      plot(gTrack(gr , stack.gap=3))


.. figure:: figure/unnamed-chunk-4-3.png
    :alt: plot of chunk unnamed-chunk-4

    plot of chunk unnamed-chunk-4
  
### gTrack(gr , y.field='GC')
##### vector or scalar numeric specifying gap between tracks (add a dimension to the data) 


::

      plot(gTrack(gr , y.field='GC'))


.. figure:: figure/unnamed-chunk-5-1.png
    :alt: plot of chunk unnamed-chunk-5

    plot of chunk unnamed-chunk-5

### gTrack(gr, bars=TRUE/FALSE)


::

      plot(gTrack(gr , y.field='GC', bars = TRUE, col = "light blue"))


.. figure:: figure/unnamed-chunk-6-1.png
    :alt: plot of chunk unnamed-chunk-6

    plot of chunk unnamed-chunk-6

### gTrack(gr , lines=TRUE/FALSE)


::

      plot(gTrack(gr , y.field='GC' , lines = TRUE , col = "purple"))


.. figure:: figure/unnamed-chunk-7-1.png
    :alt: plot of chunk unnamed-chunk-7

    plot of chunk unnamed-chunk-7

### gTrack(gr, circles = TRUE/FALSE)


::

      plot(gTrack(gr , y.field='GC' , circles = TRUE , col = "magenta" , border='60'))


.. figure:: figure/unnamed-chunk-8-1.png
    :alt: plot of chunk unnamed-chunk-8

    plot of chunk unnamed-chunk-8
 
#### colorfield 
##### map values to colors! Legend is automatically added 


::

      plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , colormaps = list(GC = c("1"="red" , "2" = "blue" , "3"="magenta", "4"="light blue" ,"5"="black" , "6"="green", "7"="brown" , "8"="pink", "9"="yellow", "10" = "orange")) ))


.. figure:: figure/unnamed-chunk-9-1.png
    :alt: plot of chunk unnamed-chunk-9

    plot of chunk unnamed-chunk-9

#### gr.colorfield 
##### same as colorfield but, gTrack will choose the mapping values for you! Good for quick plot generations of large data sets. 


::

      plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = "GC"))


.. figure:: figure/unnamed-chunk-10-1.png
    :alt: plot of chunk unnamed-chunk-10

    plot of chunk unnamed-chunk-10


#### gr.labelfield 
##### add certain meta data values to data points on graph 


::

      plot(gTrack(gr , y.field = 'GC' , bars = TRUE , col = NA , gr.colorfield = "GC", gr.labelfield="name" ))


.. figure:: figure/unnamed-chunk-11-1.png
    :alt: plot of chunk unnamed-chunk-11

    plot of chunk unnamed-chunk-11

# GRangesList



::

      # first, create GRanges object for chroms 1-3. Each chrom stores regions of exons
      chrom1 <- GRanges(seqnames=Rle(rep(1,5)) , ranges = IRanges(c(13214448,13377047,17190862,17284920,30741950) , end=c(13376489,17190004,17283075,30741656,30745210)))
    
      chrom2 <- GRanges(seqnames=Rle(rep(2,5)) , ranges = IRanges(c(34675467,34737163,50880025,50882016,51098931) , end = c(34737057,50879519,50880979,51089715,51099793)))
      
      chrom3 <- GRanges(seqnames=Rle(rep(3,5)) , ranges = IRanges(c(5883026,5888521,6651128,6655078,10251906) , end = c(5887648,6646543,6653332,10245198,10254797)))
      
      chroms <- GRangesList("chrom1" = chrom1 , chrom2 = "chrom2" , chrom3 = "chrom3")



::

    ## Error in GRangesList(chrom1 = chrom1, chrom2 = "chrom2", chrom3 = "chrom3"): all elements in '...' must be GRanges objects



