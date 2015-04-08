

=============================
MetagenomeSeq Normalisation
=============================


Count data need to be normalised to account for differences in library depth between samples. We have used `metagenomeSeq`_ 
to normalise our data as well as to assess differential abundance of genera and NOGs in both DNA and RNA data sets (do this in RNA/ and DNA/ directories). 
The normalised counts are used for downstream analysis. In order to gain a matrix of normalised counts we use the run_metagenomeseq. R script
in the R/ directory. This script produces two files - the first is the normalised counts table and the second is the 
results of the differential abundance analysis::

    $ gunzip genus.diamond.aggregated.counts.tsv.gz

    $ Rscript < R/run_metagenomeseq.R --k 4 --a 0.1 -c genus.diamond.aggregated.counts.tsv -p genus.diamond.aggregated.counts


The options specified control the number of samples (--k) that must contain the feature above a certain reads per million (RPM) (--a)
threshold. The input counts table is specified with the -c option and the -p option specifies the output file prefix. The
command above will produce the two files - genus_counts.norm.matrix and genus_counts.diff.tsv. We will use the .norm.matrix files
in the next sections for looking at overlaps between DNA and RNA data sets as well as performing Principle components analysis.

However first we load the tables into our csvdb databases::

   >> import CGAT.Pipeline as P
   >> P.load("genus.diamond.aggregated.counts.tsv", "genus.diamond.aggregated.counts.tsv.load")
   >> P.load("gene_counts.tsv", "gene_counts.tsv.load")




.. _metagenomeSeq: http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html

