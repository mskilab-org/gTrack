A reST document for knitr
=========================

Khagay - testing testing this document with automatically generated images

**knitr**:



::

    options(width = 75)
    library(knitr)
    # do not use the sourcecode directive
    render_rst(strict=TRUE)
    
    # global chunk options
    opts_chunk$set(fig.width=5, fig.height=5)



More examples
-------------

A code chunk begins with ``.. {r label, options}``, and ends with ``.. ..``
(note the space in between). Optionally you can precede all R code with two
dots, e.g.



::

     1+1



::

    ## [1] 2



::

     set.seed(123); rnorm(10)



::

    ##  [1] -0.56047565 -0.23017749  1.55870831  0.07050839  0.12928774
    ##  [6]  1.71506499  0.46091621 -1.26506123 -0.68685285 -0.44566197



::

     warning('do not forget the space after ..!')



::

    ## Warning: do not forget the space after ..!



Here is a plot:



::

    library(ggplot2)
    qplot(hp, mpg, data=mtcars) + geom_smooth()


.. figure:: figure/rst-plot-1.png
    :alt: A ggplot2 example
    :width: 360px

    A ggplot2 example

Inline R code is like this: the value of pi is 3.1415927.
