
=============================
MetagenomeSeq Normalisation
=============================

The count data we generated in the previous section needs to be normalised to account for differences in library size between samples. 
We have used `metagenomeSeq`_  to normalise our data as well as to assess differential abundance of genera and NOGs in both DNA and RNA 
data sets (do this in RNA/ and DNA/ directories). The normalised counts are used for downstream analysis. In order to gain a matrix of 
normalised counts we use the run_metagenomeseq.R script in the R/ directory of `CGATOxford/proj029`_ . This script produces two files. 
The first is the normalised counts table and the second is the results of the differential abundance analysis. For the time
being we will ignore the differential abundance results and use the normalised counts for some analysis::

    $ gunzip -c genus.diamond.aggregated.counts.tsv.gz > genus.diamond.aggregated.counts.tsv
    $ Rscript  cgat/R/run_metagenomeseq.R --k 4 --a 0.1 -c genus.diamond.aggregated.counts.tsv -p genus.diamond.aggregated.counts


The options specified control the number of samples (--k) that must contain the feature above a certain reads per million (RPM) (--a)
threshold. The input counts table is specified with the -c option and the -p option specifies the output file prefix. The
command above will produce the two files - genus_counts.norm.matrix and genus_counts.diff.tsv. We will use the .norm.matrix files
in the next sections for looking at overlaps between DNA and RNA data sets as well as performing principle components analysis (PCA).

First we load the counts tables into our csvdb databases using the csv2db.py script in the cgat/ repository. 
This is necessary for intersection queries bewteen data sets later on::

    $ zcat genus.diamond.aggregated.counts.tsv.gz | 
      python cgat/scripts/csv2db.py --backend=sqlite 
                                    --retry                              
                                    --table=genus_diamond_aggregated_counts    
                                    > genus.diamond.aggregated.counts.tsv.load

and the differential abundance tables::

    $ cat genus.diamond.aggregated.counts.diff.tsv | 
      python cgat/scripts/csv2db.py --backend=sqlite 
                                    --retry                              
                                    --table=genus_diamond_aggregated_counts_diff    
                                    > genus.diamond.aggregated.counts.diff.tsv.load


    $ cat gene_counts.diff.tsv | 
      python cgat/scripts/csv2db.py --backend=sqlite 
                                    --retry                              
                                    --table=gene_counts_diff    
                                    > gene_counts.diff.tsv.load


.. _metagenomeSeq: http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html

.. _CGATOxford/proj029: https://github.com/CGATOxford/proj029
