Vignette Showing How to Graph Mutations 
=======================================

**Draw.paths**



.. sourcecode:: r
    

    ## To create fake genes, create two GRanges objects each filled with random ranges from chromosomes 1 and 2. Ranges fall within the 1-5e3 sequence
    
    gene1 = sort(sample(gr.tile(parse.gr('1:1-5e3+'), 50), 5))
    gene2 = rev(sort(sample(gr.tile(parse.gr('2:1-5e3-'), 50), 12)))
    gene3 = sort(sample(gr.tile(parse.gr('3:1-5e3+'), 50), 8))
    
    ##Create a column that keeps track of the exons
    
    gene1$exon = 1:length(gene1)
    gene2$exon = 1:length(gene2)
    gene3$exon = 1:length(gene3)
    
    ## Combine into GRangesList
    grl = GRangesList(gene1 = gene1, gene2 = gene2, gene3 = gene3)
    
    gt.genes = gTrack(grl)
    
    ## Plot but, connect ranges using draw.paths
    fusion = GRangesList(c(grl$gene1[1:3], grl$gene2[5:9], grl$gene3[7:8]))
    gt.fusion = gTrack(fusion, draw.paths = FALSE, gr.labelfield = 'exon')
    gt.fusion.o = gTrack(fusion, draw.paths = TRUE, gr.labelfield = 'exon')
    
    
    win = parse.gr(c('1:1-1e4', '2:1-1e4', '3:1-1e4'))


.. sourcecode:: r
    

    plot(c(gt.genes, gt.fusion, gt.fusion.o), win +1e3)

.. figure:: figure/-plotList-1.png
    :alt: plot of chunk -plotList

    plot of chunk -plotList

**Example of How to Graph Fake Variants** 


.. sourcecode:: r
    

    ## Create a GRanges
    fake.genome = c('1'=1e4, '2'=1e3, '3'=5e3)
    tiles = gr.tile(fake.genome, 1)
    
    ## Choose 5 random indices 
    hotspots = sample(length(tiles), 5)
    
    d = values(distanceToNearest(tiles, tiles[hotspots]))$distance
    prob = .05 + exp(-d^2/10000)


.. sourcecode:: r
    

    mut = sample(tiles, 2000, prob = prob, replace = TRUE) 


::

    ## Error in sample.int(length(x), size, replace, prob): incorrect number of probabilities


.. sourcecode:: r
    

    win = si2gr(fake.genome)
    
    gt.mut0 = gTrack(mut, circle = TRUE, stack.gap = 0)
    gt.mut2 = gTrack(mut, circle = TRUE, stack.gap = 2)
    gt.mut10 = gTrack(mut, circle = TRUE, stack.gap = 10)
    gt.mut50 = gTrack(mut, circle = TRUE, stack.gap = 50)



.. sourcecode:: r
    

    plot(c(gt.mut0, gt.mut2, gt.mut10, gt.mut50), win)

.. figure:: figure/mutations2-plot-1.png
    :alt: plot of chunk mutations2-plot

    plot of chunk mutations2-plot


