Graphing Point Mutations
========================

To illustrate gTrack's functionality in graphing point mutations, a data set of sequences is created and a few of them will are picked as variants. This data will be graphed and because there are outliers (variants), they will be easily visable. This vignette also exemplifies how/when to use the gTrack **title** parameter.  

title Parameter
~~~~~~~~~~~~~~


.. sourcecode:: r
    

    ## create sequences from chromosomes 1-3. 
    fake.genome = c('1'=1e4, '2'=1e3, '3'=5e3)
    tiles = gr.tile(fake.genome, 1)
    
    ## Choose 5 random indices. These indices will store the variants. 
    hotspots = sample(length(tiles), 5)
    
    ## for each sequence, calculate the shortest distance to one of the hotspots.
    d = pmin(Inf, values(distanceToNearest(tiles, tiles[hotspots]))$distance, na.rm = TRUE)
    ## for sequences near the hotspots, the "prob" will be a higher positive number. It becomes smaller as it moves farther from the hotspot. 
    prob = .05 + exp(-d^2/10000)


.. sourcecode:: r
    

    ## sample 2000 of the sequences. the one nearer to the hotspots will "probably" be selected.
    mut = sample(tiles, 2000, prob = prob, replace = TRUE) 
    
    ## graph with different degrees of stack.gap. The higher numeric supplied to stack.gap helps separate the data, visually. 
    gt.mut0 = gTrack(mut, circle = TRUE, stack.gap = 0, title = "Track 0")
    gt.mut2 = gTrack(mut, circle = TRUE, stack.gap = 2, title = "Track 2")
    gt.mut10 = gTrack(mut, circle = TRUE, stack.gap = 10, title = "Track 10")
    gt.mut50 = gTrack(mut, circle = TRUE, stack.gap = 50, title = "Track 50")



.. sourcecode:: r
    

    win = si2gr(fake.genome)
    plot(c(gt.mut0, gt.mut2, gt.mut10, gt.mut50), win)


::

    ## Error in `[.data.frame`(x@formatting, , name): undefined columns selected


