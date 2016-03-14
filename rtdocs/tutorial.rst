
Tutorial 
--------

.. code-block:: bash
   
   # make sure to load gTrack
   library(gTrack)

.. code-block:: bash 

   # load data from TCGA 
   tcgaData <- read.delim("BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082.hg18.seg.txt")

   # convert data.frame to GRanges object
   tcgagr <- makeGRangesFromDataFrame(tcgaData)
   
   # wrap gTrack around TCGA GRanges object 
   tcgagt <- gTrack(tcgagr)
   
   # plot gTrack object 
   plot(tcgagt)

.. figure:: figures/tcgatoomanywindows.png
   :alt:
   :scale: 75%

.. code-block:: bash
   
   # Changed default number of windows (from entire genome to the first five windows 
   plot(tcgagt , window = tcgagr[1:5])   

.. figure:: figures/tcga5windows.png 
   :alt:
   :scale: 75%  