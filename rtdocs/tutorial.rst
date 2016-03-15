
Tutorial 
--------

.. code-block:: bash
   
   # make sure to load gTrack
   library(gTrack)

.. code-block:: bash 

   # load data from TCGA 
   tcgaData <- read.delim("inst/extdata/BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg18.seg.txt")

   # convert data.frame to GRanges object
   tcgagr <- GRanges(tcgaData)
   
   # wrap gTrack around TCGA GRanges object 
   tcgagt <- gTrack(tcgagr)
   
   # plot gTrack object 
   plot(tcgagt)

.. figure:: figures/tcgatoomanywindows.png
   :alt:
   :scale: 75%

.. code-block:: bash
   
   # Changed default number of windows (from entire genome to the first five windows). 
   # Windows argument requires a subset of a GRanges Object. Check documentation for more details. 
   plot(tcgagt , window = tcgagr[1:5])   

.. figure:: figures/tcga5windows.png 
   :alt:
   :scale: 75%  

.. code-block:: bash
   
   # use amplifications/deletions as y-values 
   tcgagt <- gTrack(tcgagr , y.field="Segment_Mean")
   plot(tcgagt , windows = tcgagrr[1:5] , col = "red")   
 
.. figure:: figures/deletions.png
   :alt:
   :scale: 75% 

.. code-block:: bash

   # add a second sample to the graph
   # create gTrack object for sample
   tcgaData2 <- read.delim("inst/extdata/BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg19.seg.txt")
   tcgagr2 <- GRanges(tcgaData2)
   tcgagt2 <- gTrack(tcgagr2 , y.field="Segment_Mean")
   
   # plot the two samples 
   plot(c(tcgagt2 , tcgagt), windows = tcgagr2[1:5] , col = "red")

.. figure:: figures/deletions2.png
   :alt:
   :scale: 75%
      
.. code-block:: bash
   
   # physically separate gaps between tracks
   plot(c(tcgagt2 , tcgagt), windows = tcgagr2[1:5] , col = "red" , ygap = 5)

.. figure:: figures/ygap.png
   :alt:
   :scale: 75%

