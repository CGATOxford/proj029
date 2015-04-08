

=======================================================
Comparing metagenomic and metatranscriptomic data sets
=======================================================


So we are going to compare datasets. First we could make a separate working directory for the comparison::

    $ mkdir comparison.dir
    $ cd comparison.dir

The first analysis that we do is to do some simple comparisons between the DNA and RNA data sets (Fig. 1 in the paper). 
This involves looking at the number of genera/NOGs that were detected in each method (and the overlap) and how abundance
as measured by DNA compares with abundance of RNA. Here we will use the PipelineMetaomics.py module
to look at the overlaps. For genera we can type::


    >> import PipelineMetaomics
    >> PipelineMetaomics.buildDetectionOverlap("<path>/RNA/genus.diamond.aggregated.counts.tsv.gz", 
                                               "<path>/DNA/genus.diamond.aggregated.counts.tsv.gz", 
                                               "genus_overlap.tsv")


This will produce an outfile named genus_overlap.tsv that contains the number of genera detected with > 1 read in >= 1 
sample in DNA and RNA analyses as well as the overlap. We run this also for the NOG count tables. 

We then compared the abundances of features (e.g. genera) that were detected in DNA, RNA or both data sets (baseed on RPM values from 
counts tables). To do this we run the follwong functions::

    >> PipelineMetaomics.plotAbundanceLevelsOfOverlap("<path>/RNA/genus.diamond.aggregated.counts.tsv.gz",
                                                      "<path>/DNA/genus.diamond.aggregated.counts.tsv.gz",
                                                      "genus_abundance_distributions.png")

This produces the following plot.

.. image:: ../images/genus_abundance_distributions.png
    :align: center
    :width: 300pt
    :height: 300pt

This function also produces the file "genus_abundance_distributions.sig" which contains the significance of the differences
between sets of features based on the wilkoxon rank sum test. 


We can then compare abundance estimates for those commonly detected genera (or NOGs) based on metagenomeSeq normalised
abundances.:: 
 

    >> PipelineMetaomics.scatterplotAbundanceEstimates("<path>DNA/genus.diamond.aggregated.counts.norm.matrix",
                                                       "<path>RNA/genus.diamond.aggregated.counts.norm.matrix",
                                                       "genus_abundance_correlation.pdf")


This produces the plot below in the file genus_abundance_correlation.pdf and prints out our correlation coefficient of 0.98.


.. image:: ../images/genus_abundance_correlations.png
    :align: center
    :width: 300pt
    :height: 300pt


Again, this was repeated using the same functions for NOG-based counts and normalised counts (Fig. 2).





