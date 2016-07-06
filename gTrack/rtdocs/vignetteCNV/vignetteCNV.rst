Vignette Using CNV Data
=======================


.. sourcecode:: r
    

    #### make sure to load gUtils
    library(gTrack)
    #### load in SCNA data
    #### these are level 3 segments from over 1600 TCGA breast cancer cases
    #### downloaded from the TCGA and converted to a single GRanges object 
    scna = seg2gr(system.file(c('extdata/files'), package = 'gTrack'))


::

    ## Error in if (nrow(segs) == 0) return(GRanges(seqlengths = seqlengths(seqinfo))): argument is of length zero


.. sourcecode:: r
    

    #### we define amplification events and deletion events
    #### and then do a simple recurrence analysis on amplifications
    
    #### amplification events are defined using %Q% gUtils operator to filter
    #### on scna metadata for segments with seg.mean greater than 1
    #### greater than 50 markers and width less than 1MB
     amps = scna %Q% (seg.mean<(-1) & num.mark > 50 & width < 1e7)


::

    ## Error in scna %Q% (seg.mean < (-1) & num.mark > 50 & width < 1e+07): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error: object 'scna' not found


.. sourcecode:: r
    

    #### apply a similar filter to define deletions
     dels = scna %Q% (seg.mean<(-1) & num.mark > 50 & width < 1e7)


::

    ## Error in scna %Q% (seg.mean < (-1) & num.mark > 50 & width < 1e+07): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error: object 'scna' not found


.. sourcecode:: r
    

    #### compute the amplification "score" as the total number of amplification
    #### events in a given region
     amp.score = as(coverage(amps), 'GRanges')


::

    ## Error in coverage(amps): error in evaluating the argument 'x' in selecting a method for function 'coverage': Error: object 'amps' not found


.. sourcecode:: r
    

    #### define the peaks of amplification as the 100 regions with
    #### the highest amplification score
     amp.peaks = amp.score %Q% (rev(order(score))) %Q% (1:100)


::

    ## Error in amp.score %Q% (rev(order(score))) %Q% (1:100): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error in amp.score %Q% (rev(order(score))) : 
    ##   error in evaluating the argument 'x' in selecting a method for function '%Q%': Error: object 'amp.score' not found


.. sourcecode:: r
    

    #### "reduce" or merge the top peaks to find areas of recurrent amplification
    amp.peaks = reduce(amp.peaks+1e5) %$% amp.peaks %Q% (rev(order(score)))


::

    ## Error in reduce(amp.peaks + 1e+05) %$% amp.peaks %Q% (rev(order(score))): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error in reduce(amp.peaks + 1e+05) %$% amp.peaks : 
    ##   error in evaluating the argument 'x' in selecting a method for function '%$%': Error in reduce(amp.peaks + 1e+05) : 
    ##   error in evaluating the argument 'x' in selecting a method for function 'reduce': Error: object 'amp.peaks' not found


.. sourcecode:: r
    

    ### do a similar analysis for dels
    del.score = as(coverage(dels), 'GRanges')


::

    ## Error in coverage(dels): error in evaluating the argument 'x' in selecting a method for function 'coverage': Error: object 'dels' not found


.. sourcecode:: r
    

    del.peaks = del.score %Q% (rev(order(score))) %Q% (1:100)


::

    ## Error in del.score %Q% (rev(order(score))) %Q% (1:100): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error in del.score %Q% (rev(order(score))) : 
    ##   error in evaluating the argument 'x' in selecting a method for function '%Q%': Error: object 'del.score' not found


.. sourcecode:: r
    

    del.peaks = reduce(del.peaks+1e5) %$% del.peaks %Q% (rev(order(score)))


::

    ## Error in reduce(del.peaks + 1e+05) %$% del.peaks %Q% (rev(order(score))): error in evaluating the argument 'x' in selecting a method for function '%Q%': Error in reduce(del.peaks + 1e+05) %$% del.peaks : 
    ##   error in evaluating the argument 'x' in selecting a method for function '%$%': Error in reduce(del.peaks + 1e+05) : 
    ##   error in evaluating the argument 'x' in selecting a method for function 'reduce': Error: object 'del.peaks' not found


.. sourcecode:: r
    

    #### load in GRanges of GENCODE genes
    genes = readRDS(system.file("extdata", 'genes.rds', package = "gTrack"))


::

    ## Warning in gzfile(file, "rb"): cannot open compressed file '', probable
    ## reason 'No such file or directory'



::

    ## Error in gzfile(file, "rb"): cannot open the connection


.. sourcecode:: r
    

    #### use %$% operator to annotate merged amp and del peaks with "gene name" metadata
    amp.peaks = amp.peaks %$% genes[, 'gene_name']


::

    ## Error in amp.peaks %$% genes[, "gene_name"]: error in evaluating the argument 'x' in selecting a method for function '%$%': Error: object 'amp.peaks' not found


.. sourcecode:: r
    

    del.peaks = del.peaks %$% genes[, 'gene_name']


::

    ## Error in del.peaks %$% genes[, "gene_name"]: error in evaluating the argument 'x' in selecting a method for function '%$%': Error: object 'del.peaks' not found


.. sourcecode:: r
    

    ### now that we've computed scores and annotated peaks
    ### we want to inspect these peaks and plot them with gTrack
    
    ### load in precomputed gTrack of hg19 GENCODE annotation
    ### (note this is different from the GENCODE genes which is a GRanges
    ### we loaded in a previous line .. this is purely for visualization)
    ge = track.gencode()


::

    ## Pulling gencode annotations from /data/research/mski_lab/Software/R/gTrack/extdata/gencode.composite.collapsed.rds


.. sourcecode:: r
    

    #### build a gTrack of amps colored in red with black border
    #### and one of dels colored in blue 
    gt.amps = gTrack(amps,  col = 'red', name = 'Amps')


::

    ## Error in listify(data, GRanges): object 'amps' not found


.. sourcecode:: r
    

    gt.dels = gTrack(dels, col = 'blue', name = 'Dels')


::

    ## Error in listify(data, GRanges): object 'dels' not found


.. sourcecode:: r
    

    #### build a gTrack of amp and del score as a line plot
    gt.amp.score = gTrack(amp.score, y.field = 'score',
        col = 'red', name = 'Amp score', line = TRUE)


::

    ## Error in listify(data, GRanges): object 'amp.score' not found


.. sourcecode:: r
    

    gt.del.score = gTrack(del.score, y.field = 'score',
        col = 'blue', name = 'Amp score', line = TRUE)


::

    ## Error in listify(data, GRanges): object 'del.score' not found


.. sourcecode:: r
    

    #### build a gTrack of peaks of amp and del peaks
    gt.amp.peaks = gTrack(amp.peaks, gr.labelfield = 'gene_name', 
        col = 'pink', border = 'black', name = 'Amp peaks', height = 5)


::

    ## Error in listify(data, GRanges): object 'amp.peaks' not found


.. sourcecode:: r
    

    gt.del.peaks = gTrack(del.peaks, gr.labelfield = 'gene_name',
        col = 'lightblue', border = 'black', name = 'Amp peaks', height = 5)


::

    ## Error in listify(data, GRanges): object 'del.peaks' not found


.. sourcecode:: r
    

    ### let's look at the top amplification peak
    amp.peaks[1]


::

    ## Error in eval(expr, envir, enclos): object 'amp.peaks' not found


.. sourcecode:: r
    

    ### interesting! this looks like a novel peak with genes that have
    ### not previously been associated with breast cancer
    ### ("RP11-122L4.1, AC123767.1, CTD-2024D23.1, ADAM18, ADAM2")
    
    ### let's look at the data supporting this peak - including
    ### the underlying amp events, amp score, and peak region boundary

