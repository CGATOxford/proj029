
=================================
MetagenomeSeq filtered analysis
=================================


Now we have analysed the RNA and DNA data sets in isolation we can start pulling them together. The first
thing to do is to re-run the metagenomeSeq analysis on a filtered dataset. This means restricting the
analysis to those features that were present at > 0.1RPM in both DNA and RNA data sets.

For genera we can do::
 
    $ cd RNA
    $ zcat genus.diamond.aggregated.counts.tsv.gz | python proj029/scripts/counts2restrictedcounts.py 
                                                  --restrict-to=../compare_datasets/common_genera.tsv
                                                  --log=restrict_genera.log
                                                  > genus.diamond.aggregated.counts.restricted.tsv


We can then use this filtered counts table for metagenomeSeq analysis::

    $ Rscript cgat/R/run_metagenomeseq.R --k 4 --a 0.1 -c genus.diamond.aggregated.counts.restricted.tsv -p genus.diamond.aggregated.counts.restricted







