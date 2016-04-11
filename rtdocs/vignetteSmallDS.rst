===============================
Vignette Using a Small Data Set
===============================

.. code-block:: bash 

   library(gTrack)
   library(gUtils)


.. code-block:: bash

   ## Given a GRanges object, separate the interval by "w".
   ## In this example, the output will be 20 equal sized tiles

   gr <- gr.tile(GRanges(1, IRanges(1,100)), w=5)

   ## Plot how the tiles would like
   plot(gTrack(gr))

.. figure:: figures/tilePlot.png

   :alt:
   :scale: 75% 

.. code-block:: bash

  ## "1" in this example adds one more title to each interval
  plot(gTrack(gr+1))

