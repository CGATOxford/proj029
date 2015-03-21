.. highlight:: rst

==============
Introduction
==============

The purpose of this documentation is to allow those that were interested enough to 
read our metagenomics paper to follow the analyses that we used step by step and
recreate results and figures from our analyses. As with all high-throughput
sequencing experiments a number of tasks take a long time to run and so where
possible we have provided the output files so that these tasks can be skipped. We
will however show you the code that generated the output.


Overview
=========

The code for the analyses is distributed across multiple locations. This is due
to the fact that some of it is R code, some python, some is used in scripts and
some is part of module files for pipelines. The directory structure should enable 
people to easily locate the source code i.e. modules are in modules/ scripts are
in scripts/ etc. etc. For ouur purposes and to enable re-running of analyses and
facilitate reproducibity, scripts and module functions are wrapped up in `ruffus`_
pipelines. We also extensively use functions and scripts that are part of the 
`CGAT code collection`_. 

However, this documentation will take you through the running of 3rd party tools,
custom scripts and functions for analyses without having to specifically set up
your environment for using `CGAT pipelines`_


.. _ruffus: http://www.ruffus.org.uk/

.. _CGAT code collection: https://github.com/CGATOxford/cgat

.. _CGAT pipelines: https://www.cgat.org/downloads/public/cgat/documentation/UsingPipelines.html


Dependencies
=============

To run the analyses you will have to install a bit of software. Shown are the versions
used in the analysis although newer versions may now be available and could be substituted




+------------+--------------+
| Software   | Version used |
+============+==============+
|   R        | 3.1.0        |
+------------+--------------+
|   Python   | 2.7.1        |
+------------+--------------+
|   DIAMOND  | 0.3.9        |
+------------+--------------+
|   mtools   | `here`_      |
+------------+--------------+
|  fastx     |              |
|  toolkit   | 0.0.13       |
+------------+--------------+


R libraries

+---------------+
| R libraries   |
+===============+
| ggplot2       |
+---------------+
| metagenomeSeq |
+---------------+
| gplots        |
+---------------+
| gtools        |
+---------------+
| pheatmap      |
+---------------+
| vegan         |
+---------------+


Python libraries


+--------------------+
|  python libraries  |
+====================+
|      pandas        |
+--------------------+
|      sqlite3       |
+--------------------+
|      numpy         |
+--------------------+
|      rpy2          |
+--------------------+


.. _here: http://ab.inf.uni-tuebingen.de/data/software/megan5/download/mtools.zip


