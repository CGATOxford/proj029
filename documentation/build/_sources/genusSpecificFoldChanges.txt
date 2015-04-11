
=====================================
Estimating genus-level fold changes
=====================================

In the previous section we looked at the proportion of reads from each genus that
were assigned to each NOG. This was an average across samples and we were further
interested in defining which genera contributed the most to our observed changes in
NOG expression. To define genus-level fold changes between colitic mice and steady-state
mice took two steps. In step 1 we used the script nogs2genera.py as we did previously to 
produce a per genus-NOG count (as opposed to percentage). 


so for example::

    $ python <path_to_proj029>/scripts/nogs2genera.py
             -m gene2cog.tsv.gz
             -d common_genes.tsv
             --level=genus
             --counts
             --alignment-taxa=stool-HhaIL10R-R1.lca
             --alignment-genes=stool-HhaIL10R-R1.igc.tsv.gz                    
             --log=stool-HhaIL10R-R1.diamond.ptaxa.tsv.gz.log
             | gzip > associate_taxa.dir/stool-HhaIL10R-R1.diamond.ctaxa.tsv.gz


The only difference here is that we use the --counts switch so that percentages aren't returned. Once again
each sample for both DNA and RNA are run through this process and count tables are combined with
cgat/combine_tables.py. This is a time-consuming and memory intensive step. The files are therefore
provided in the data/RNA and data/DNA directories - called associated_taxa_counts.tsv.gz

The second step is straightforward and draws on another script that we used earlier. Each count table (DNA and RNA)
is run through metagenomeSeq.R to normalise the data and to estimate fold changes. The .diff.tsv files are
then used to plot DNA fold changes vs. RNA fold changes at the genus level. We include the gene_counts.diff files
so that we can include the overall NOG fold changes in the plots. 

Create the genus-level .diff.tsv files for RNA and DNA. E.g.::

    $ cd RNA
    $ ln -s <path_to_data>/data/RNA/associated_taxa_counts.tsv.gz .
    $ gunzip -c associated_taxa_counts.tsv.gz > associated_taxa_counts.tsv
    $ Rscript <path_to_cgat>/cgat/R/run_metagenomeseq.R --k 4 --a 0.1 -c associated_taxa_counts.tsv -p associated_taxa_counts

Once this is done for RNA and DNA we can plot the fold changes. Apologies but you will have to create a subfolder called
scatterplot_genus_cog_fold.dir in compare_datasets::

    $ mkdir scatterplot_genus_cog_fold.dir

    >> import Proj029Pipelines.PipelineMetaomics as PipelineMetaomics
    >> PipelineMetaomics.scatterplotPerCogTaxaDNAFoldRNAFold(<path_to_rna>/RNA/associated_taxa_counts.diff.tsv,
                                                             <path_to_dna>/DNA/associated_taxa_counts.diff.tsv,
                                                             <path_to_rna>/RNA/gene_counts.diff.tsv,
                                                             <path_to_rna>/RNA/gene_counts.diff.tsv)



This will produce plots of our genes of interest i.e. those involved in oxidative stress
resistance and glycan utilisation.

As an example here is the plot for COG0783: Dps/Ferritin

.. image:: ../images/COG0783.scatters.png
    :align: center
    :width: 400pt
    :height: 400pt

This is where we end the analysis of the microbiota. The next analyses are on the host response to colitis.
We are not going into the anlayses of these data as they are fairly standard i.e. LIMMA to perform microarray
analysis and DESeq for RNA-seq analysis.

We hope that this documentation was useful!

