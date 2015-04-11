.. highlight:: rst

==============
Introduction
==============

The aims of the project and analysis are explained in the paper. Briefly, we were 
interested to discover community structure and functional changes in the microbiota associated with 
inflammation in a mouse model of Helicobacter hepaticus + aIL10R colitis. To this
end we used shotgun metagenomic and metatranscriptomic sequencing. Using these two
methods we were interested in comparing RNA- and DNA-based methods for microbiome
profiling and integrating the two to indentify functional groups that are the most
likely to represent colitis-induced functions.


The purpose of this documentation is to allow those that were interested enough to 
read our metagenomics paper to follow the analyses that we used step by step and
recreate results and figures from our analyses. As with all high-throughput
sequencing experiments, a number of tasks take a long time to run and created files are
large. While we expect that the analysis can be run through from start to finish it may be
more apropriate to only run certain sections. Therefore we have provided the output files 
so that some tasks can be skipped.


Overview
=========

The code for the analyses is distributed across multiple locations. This is due
to the fact that some of it is R code, some python, some is used in scripts and
some is part of module files for pipelines. Modules and scripts are also used
from various repositories. The directory structure should enable people to easily 
locate the source code i.e. modules are in modules/ scripts are
in scripts/ etc. For our purposes and to enable re-running of analyses and
facilitate reproducibity, scripts and module functions are wrapped up in `ruffus`_
pipelines. We also extensively use functions and scripts that are part of the 
`CGAT code collection`_. This documentation will take you through the running of 3rd party tools,
custom scripts and functions for analyses without having to specifically set up
your environment for using `CGAT pipelines`_. Nevertheless, pipelines are publically
available at `CGATOxford/CGATPiplines`_ and `CGATOxford/proj029/Proj029Pipelines`_ 

For all analyses it is recommended that you follow the naming of directories, files etc 
as they are given. This is because there are a number of steps where these are
hardcoded. 


Disclaimer
----------

The code used in the analysis was not designed for running on all systems. As such 
we know that it works on our system (Red Hat Enterprise Linux Server release 6.6 (Santiago)).

.. _ruffus: http://www.ruffus.org.uk/

.. _CGAT code collection: https://github.com/CGATOxford/cgat

.. _CGAT pipelines: https://www.cgat.org/downloads/public/cgat/documentation/UsingPipelines.html

.. _CGATOxford/CGATPiplines: https://github.com/CGATOxford/CGATPipelines

.. _CGATOxford/proj029/Proj029Pipelines: https://github.com/CGATOxford/proj029/


Dependencies
=============

To run the analyses you will have to install a bit of software. Shown are the versions
used in the analysis although newer versions may now be available and could be substituted.

First of all make sure that you have the CGAT code collection installed - this includes cgat and
CGATPipelines. It is recommended that with third party python modules that you use a virtual environment::
  
    # install virtualenv
    pip install virtualenv
    
    # setup virtual environment 
    virtualenv --no-site-packages metagenomics

    # start using the virtual environment
    source metagenomics/bin/activate 

    # clone the CGAT repository and install (scripts and modules)
    git clone https://github.com/CGATOxford/cgat.git
    cd cgat
    python setup.py install

    # clone CGAT pipelines (we use Pipeline.py module extensively)
    git clone https://github.com/CGATOxford/CGATPipelines.git
    cd CGATPipelines
    python setup.py install
    
You will also have to install the proj029 respository::

    git clone https://github.com/CGATOxford/proj029.git
    cd proj029
    python setup.py install



3rd party software
-------------------

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
|   fastx    |              |
|   toolkit  | 0.0.13       |
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


