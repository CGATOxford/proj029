
=================================
MetagenomeSeq filtered analysis
=================================


Now we have compared the RNA and DNA data sets and have sets of features that we want to analyse,
we can run metagenomeSeq on a filtered set of features (using common_genera.tsv and common_genes.tsv).

For example to produce new counts tables with filtered features for genera we can do::
 
    $ cd RNA
    $ zcat genus.diamond.aggregated.counts.tsv.gz | python proj029/scripts/counts2restrictedcounts.py 
                                                  --restrict-to=../compare_datasets/common_genera.tsv
                                                  --log=restrict_genera.log
                                                  > genus.diamond.aggregated.counts.restricted.tsv


We can then use this filtered counts table for metagenomeSeq analysis::

    $ Rscript cgat/R/run_metagenomeseq.R --k 4 --a 0.1 -c genus.diamond.aggregated.counts.restricted.tsv -p genus.diamond.aggregated.counts

::
    NOTE: here we have overwritten the previous results tables from metagenomeSeq.


The resulting .norm.matrix files are used in subsequent principle components analysis and the .diff.tsv files are the final
results of differential abundance testing.

Load the differential abundance table into the database::

    $ cat genus.diamond.aggregated.counts.diff.tsv | 
      python <path_to_cgat>/cgat/scripts/csv2db.py --backend=sqlite 
                                                   --retry                              
                                                   --table=genus_diamond_aggregated_counts_diff    
                                                   > genus.diamond.aggregated.counts.diff.tsv.load



::
    NOTE: here we have overwritten the database tables






