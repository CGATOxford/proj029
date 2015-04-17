
====================================
MetagenomeSeq initial normalisation
====================================

The count data we generated in the previous section needs to be normalised to account for differences in library size between samples. 
We have used `metagenomeSeq`_  to normalise our data as well as to assess differential abundance of genera and NOGs in both DNA and RNA 
data sets (do this in RNA/ and DNA/ directories). 


In this step we will remove any features (genera/NOGs) whose abundance is < 0.1 reads per million (RPM). Because we are running the
RNA and DNA analyses separately there is the issue that not all features are present in both data sets. Therefore, in this initial
normalisation step we will produce the normalised data for each data set separately, compare the two data sets in the next section 
and then run metagenomeSeq a second time using RNA and DNA counts tables that are restricted to those features that are present in both 
at an abundance > 0.1RPM.

To create normalised counts we use the run_metagenomeseq.R script in the R/ directory of `CGATOxford/cgat`_ . This script produces two files. 
The first is the normalised counts table (ends with .norm.matrix) and the second is the results of the differential abundance analysis (ends
with .diff.tsv).

Do this for the relevant genera and gene counts tables::

    $ gunzip -c genus.diamond.aggregated.counts.tsv.gz > genus.diamond.aggregated.counts.tsv
    $ <path_to_R_install>/Rscript <path_to_cgat>/cgat/R/run_metagenomeseq.R --k 4 --a 0.1 -c genus.diamond.aggregated.counts.tsv -p genus.diamond.aggregated.counts


The options specified control the number of samples (--k) that must contain the feature above a certain reads per million (RPM) (--a)
threshold. The input counts table is specified with the -c option and the -p option specifies the output file prefix. The
command above will produce the two files - genus.diamond.aggregated.counts.norm.matrix and genus.diamond.aggregated.counts.diff.tsv.
We will use the .norm.matrix files in the next section for looking at overlaps between DNA and RNA data sets and correlating abundance estimates.

We load the differential abundance tables into our csvdb databases using the csv2db.py script in the cgat/ repository. 
This is necessary for intersection queries bewteen data sets in the next section.

In both RNA/ and DNA/ directories do::


    $ cat genus.diamond.aggregated.counts.diff.tsv | 
      python <path_to_cgat>/cgat/scripts/csv2db.py --backend=sqlite 
                                                   --retry                              
                                                   --table=genus_diamond_aggregated_counts_diff    
                                                   > genus.diamond.aggregated.counts.diff.tsv.load


    $ cat gene_counts.diff.tsv | 
      python <path_to_cgat>/cgat/scripts/csv2db.py --backend=sqlite 
                                                   --retry                              
                                                   --table=gene_counts_diff    
                                                   > gene_counts.diff.tsv.load


.. _metagenomeSeq: http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html

.. _CGATOxford/cgat: https://github.com/CGATOxford/cgat
