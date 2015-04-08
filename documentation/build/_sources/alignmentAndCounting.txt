===============================
Alignment and counting
===============================

We are picking up the analysis at the point at which raw fastq files have
been filtered for contaminating adapters, pairs have been flashed and reads
mapping to rRNA and mouse have been removed. Details of these steps are in the 
paper. The RNA and DNA analyses were run separately using the same pipeline to avoid 
filename conflicts. Therefore if you want to follow some of the analyses then
it is probably useful to create separate working directories e.g::

    $ mkdir RNA
    $ mkdir DNA 


Alignment
==========

The first step in the analysis is to assign reads to taxa and functional groups.
In order to do this we need some databases. We have used the NCBI non-redundant
protein database (`nr`_) for taxanomic purposes and the integrated gene catalogue
(`IGC`_) for assessing functions.

.. _nr: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

.. _IGC: ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.pep.gz


We used DIAMOND to align sequences to the reference databases so first we built
database indexes::

    $ diamond makedb --db nr --in nr --threads 16

and::

    $ diamond makedb --db igc --in IGC.pep --threads 16



NOTE: DIAMOND takes uncompressed files as input


Next we ran the alignment for each fastq file (reads) to each database 
(fastq files need to be converted to fasta format first e.g. using `fastx toolkit`_). As 
an example:

conversion::

    $ fastq-to-fasta stool-HhaIL10R-R1.fastq.gz > stool-HhaIL10R-R1.fasta 

alignment::

    $ diamond blastx --db nr --query stool-HhaIL10R-R1.fasta -o stool-HhaIL10R-R1.diamond.tsv


This was done for each sample (DNA-seq and RNA-seq) against each data base resulting in

16x DNA-seq to nr

16x RNA-seq to nr

16x DNA-seq to IGC

16x RNA-seq to IGC


Counting
=========

Our main objective was to identify community structure and functional changes between mice
with and without colitis. In order to do that we first needed to count the number of
reads that mapped to each genus and each functional category. This is where the analysis
of taxonomy and function diverge. 


For taxnomic profiling We used the lowest common ancestor approach (LCA) to 
assign reads to genera implemented using lcamapper.sh from mtools (see dependencies). This 
requires the mapping file of gi number to taxonomy id that is distributed with mtools so for each 
sample for both DNA and RNA data sets we do for example::

    $ lcamapper.sh -i stool-HhaIL10R-R1.diamond.tsv -f Detect -ms 50 -me 0.01 -tp 50 -gt gi_taxid_prot.bin -o stool-HhaIL10R-R1.lca

Again each file that was aligned to the ncbi nr database is used as input to lcamapper.sh. To obtain counts per
genus we use the lca2table.py script that is in the scripts/ directory of this repository::

    $ cat stool-HhaIL10R-R1.lca | python scripts/lca2table.py --summarise=taxa-counts --log=stool-HhaIL10R-R1.counts.tsv.log > stool-HhaIL10R-R1.counts.tsv

 
At this point each file that contains taxa counts is loaded into an sqlite database. This makes subsetting etc easier
downstream. We use the Pipeline.py module from the CGAT code collection to load the tables. As I mentioned
these analyses are wrapped up in ruffus pipelines and are therefore in python scripts. However to 
load the table from python we simply do::


    >> import CGAT.Pipeline as P
    >> P.load("stool-HhaIL10R-R1.lca", "stool-HhaIL10R-R1.lca.load")

This will create a database called "csvdb" in the working directory and will have loaded the table
stool_HhaIL10R_R1_lca. In our analyses we were interested in genus abundances and so we create
flat files from the database simply by::

   $ sqlite3 csvdb 'SELECT taxa, count FROM stool_HhaIL10R_R1_lca WHERE taxa == "genus"' > stool-HhaIL10R-R1.lca.counts.tsv


Next we combine counts for each sample into a single table with rows as genera and samples as columns. We use a convenient
script from CGAT code collection to do this - combine_tables.py.

e.g for example if we are in the RNA analysis working directory we run combine_tables.py by specifying that we want missing
genera from a sample to be given a 0 count, we want to merge on column1 (genera), we want to take the "count"
column and we want to combine all tables that end in .lca.counts.tsv. We also specify the prefixes to be used
for each column in the resulting combined table::


    $ python cgat/scripts/combine_tables.py                    
             --missing=0                    
             --columns=1                    
             --take=count                    
             --glob=*.lca.counts.tsv 
             --prefixes=stool-HhaIL10R-R4.lca,
                        stool-HhaIL10R-R3.lca, 
                        stool-Hh-R4.lca,
                        stool-Hh-R3.lca,
                        stool-WT-R4.lca,
                        stool-aIL10R-R1.lca,
                        stool-WT-R3.lca,
                        stool-WT-R2.lca,
                        stool-aIL10R-R4.lca,
                        stool-Hh-R2.lca,
                        stool-Hh-R1.lca,
                        stool-aIL10R-R2.lca,
                        stool-WT-R1.lca,
                        stool-HhaIL10R-R1.lca,
                        stool-HhaIL10R-R2.lca,
                        stool-aIL10R-R3.lca 
     | gzip > genus.diamond.aggregated.counts.tsv.gz


Again we do this for both RNA and DNA data sets. 


The next task is to produce a counts table similar to the one above but for functions. We have aligned to the IGC and we use
their annotations of eggNOG functions (NOGs) as the features to be counted. For counting per NOG, we extract the best hit DIAMOMD alignment
for eac read, map the gene to NOG using an additional mapping file and count. To produce a count table for a single sample
we use the diamond2counts.py script in the scripts/ directory. The input is the DIAMOND alignment file. ::


    $ zcat stool-HhaIL10R-R1.igc.tsv.gz | python scripts/diamond2counts.py                     
                                                 --method=best 
                                                 --cog-map=data/gene2cog.tsv.gz 
                                                 --sum-cog                    
                                                 --log=stool-HhaIL10R-R1.igc.counts.tsv.gz.log                    
                                        | gzip > stool-HhaIL10R-R1.nogs.counts.tsv.gz


Again, we combine tables for each sample into a final counts table using combine_tables.py. Now we have count tables for genera and
NOGs for both metagenomic and metatranscriptomic data we can start doing some analysis.


So to recap, in our working directories, DNA/ and RNA/ we now have count tables for alignments to genera and NOGs. Given the large 
sizes of raw and alignment files these have been deposited at the EBI ENA (ADD LINK). However for reproducing our downstream analysis
we have provided the count tables in the data/ directory.



.. _fastx toolkit: http://hannonlab.cshl.edu/fastx_toolkit/






