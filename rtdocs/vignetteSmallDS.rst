Vignette Using a Small Data Set
===============================

Overlapping Segments

gr.tile(gr , w)

tile ranges across GRanges


.. sourcecode:: r
    

    ## Given a GRanges object, separate the interval by "w".
    ## In this example, the output will be 20 equal sized tiles
    
    print(GRanges(1 , IRanges(1,100)))


::

    ## GRanges object with 1 range and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]        1  [1, 100]      *
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths


.. sourcecode:: r
    

    gr <- gr.tile(GRanges(1, IRanges(1,100)), w=5)



.. sourcecode:: r
    

    ## Plot how the tiles would like
    plot(gTrack(gr))

.. figure:: figure/rst-plot-1.png
    :alt: plot of chunk rst-plot

    plot of chunk rst-plot

Overlapping Tiles - gTrack(gr + n)
indicate the degree of overlap


.. sourcecode:: r
    

    plot(gTrack(gr+10))

.. figure:: figure/rst-plot2-1.png
    :alt: plot of chunk rst-plot2

    plot of chunk rst-plot2




