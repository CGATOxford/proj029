===============================
Alignment and counting
===============================

We are picking up the analysis at the point at which raw fastq files have
been filtered for contaminating adapters, pairs have been flashed and reads
mapping to rRNA and mouse have been removed. Details of these steps are in the 
paper. Both raw fastq files and processed fastq files are available at the EBI ENA
under accession number EMTAB-XXXX so the analysis can be run without having to 
perform the pre-processing steps yourself::

    NOTE: The alignment and counting steps can take a long time and large amount of memory. If
    you do not want to run these steps then we have provided the counts tables for genera
    and functional groups (NOGs) - see bottom of this page.


The RNA and DNA analyses were run separately using the same pipeline. For the our
purposes we will run the analyses separately in two separate directories. Create these
directories::

    $ mkdir RNA
    $ mkdir DNA 


Now download the .fastq files from the EBI ENA into their respective directories.



Alignment
==========

The first step in the analysis is to assign reads to taxa and functional groups.
In order to do this we need some databases. We have used the NCBI non-redundant
protein database (`nr`_) for taxonomic purposes and the integrated gene catalogue
(`IGC`_) for assessing functions. Download these into a directory of your choice::

    $ mkdir databases
    $ cd databases
    $ wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
    $ wget ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.pep.gz

.. _nr: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

.. _IGC: ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.pep.gz


We used DIAMOND to align sequences to the reference databases so first we build
database indexes::

    $ diamond makedb --db nr --in nr --threads 16

and::

    $ diamond makedb --db igc --in IGC.pep --threads 16


NOTE: DIAMOND takes uncompressed files as input so you will need to unzip first

Next we ran the alignment for each fastq file (reads) to each database 
(fastq files need to be converted to fasta format first e.g. using `fastx toolkit`_). As 
an example:

    # change into the RNA directory (contains raw data) for example
    $ cd RNA

conversion::

    $ fastq-to-fasta stool-HhaIL10R-R1.fastq.gz > stool-HhaIL10R-R1.fasta 

alignment::

    $ diamond blastx --db ../databases/nr --query stool-HhaIL10R-R1.fasta -o stool-HhaIL10R-R1.diamond.tsv


This was done for each sample (DNA-seq and RNA-seq) against each database, resulting in

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

For taxonomic profiling We used the lowest common ancestor approach (LCA) to 
assign reads to genera implemented using lcamapper.sh from mtools (see dependencies). This 
requires the mapping file of gi number to taxonomy id that is distributed with mtools so for each 
sample for both DNA and RNA data sets we did for example::

    $ lcamapper.sh -i stool-HhaIL10R-R1.diamond.tsv -f Detect -ms 50 -me 0.01 -tp 50 -gt gi_taxid_prot.bin -o stool-HhaIL10R-R1.lca

Again each file that was aligned to the ncbi nr database is used as input to lcamapper.sh. To obtain counts per
genus (we used genus level analysis throughout) we use the lca2table.py script that is in the scripts/ directory 
of the `CGATOxford/cgat`_ repository::

    $ cat stool-HhaIL10R-R1.lca | python cgat/scripts/lca2table.py --summarise=taxa-counts --log=stool-HhaIL10R-R1.lca.counts.tsv.log > stool-HhaIL10R-R1.lca.counts.tsv


.. _CGATOxford/cgat: https://github.com/CGATOxford/cgat 

 
At this point each file that contains taxa counts is loaded into an sqlite database. This makes subsetting etc easier
downstream. We use the c2v2db.py script from cgat/ respository to do this::

    $ cat stool-HhaIL10R-R1.counts.tsv | python cgat/scripts/csv2db.py
                                    --backend=sqlite 
                                    --retry                              
                                    --table=stool_HhaIL10R_R1_lca_counts
                                    > stool-HhaIL10R-R1.lca.counts.tsv.load


.. _CGATOxford/CGATPipelines: https://github.com/CGATOxford/CGATPipelines 

This will create a database called "csvdb" in the working directory and will have loaded the table
stool_HhaIL10R_R1_lca. In our analyses we were interested in genus abundances and so we create
flat files from the database simply by::

   $ sqlite3 csvdb 'SELECT taxa, count FROM stool_HhaIL10R_R1_lca_counts WHERE taxa == "genus"' > stool-HhaIL10R-R1.lca.counts.tsv


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
             --prefixes=stool-HhaIL10R-R4,
                        stool-HhaIL10R-R3, 
                        stool-Hh-R4,
                        stool-Hh-R3,
                        stool-WT-R4,
                        stool-aIL10R-R1,
                        stool-WT-R3,
                        stool-WT-R2,
                        stool-aIL10R-R4,
                        stool-Hh-R2,
                        stool-Hh-R1,
                        stool-aIL10R-R2,
                        stool-WT-R1,
                        stool-HhaIL10R-R1,
                        stool-HhaIL10R-R2,
                        stool-aIL10R-R3 
     | gzip > genus.diamond.aggregated.counts.tsv.gz


Again we do this for both RNA and DNA data sets. This produces a table (truncated for visual reasons)
    
    +-----------------+-----------------+-----------------+---------------------+-----------------+-----------------+-----------------------+
    |taxa             |stool-WT-R1_count|stool-WT-R3_count|stool-aIL10R-R1_count|stool-Hh-R2_count|stool-Hh-R1_count|stool-HhaIL10R-R4_count|
    +=================+=================+=================+=====================+=================+=================+=======================+
    |Methylobacillus  |228              |517              |560                  |406              |201              |353                    |
    +-----------------+-----------------+-----------------+---------------------+-----------------+-----------------+-----------------------+
    |Methanosphaera   |98               |224              |194                  |175              |65               |132                    |
    +-----------------+-----------------+-----------------+---------------------+-----------------+-----------------+-----------------------+
    |Desulfarculus    |3                |6                |12                   |2                |5                |2                      |
    +-----------------+-----------------+-----------------+---------------------+-----------------+-----------------+-----------------------+
    |Polaromonas      |859              |2021             |2034                 |1111             |616              |1806                   |
    +-----------------+-----------------+-----------------+---------------------+-----------------+-----------------+-----------------------+
    |Caldanaerobacter |3330             |5367             |5847                 |3645             |2072             |6571                   |
    +-----------------+-----------------+-----------------+---------------------+-----------------+-----------------+-----------------------+



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


Again, we combine tables for each sample into a final counts table using combine_tables.py to give

    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |ref         |stool-HhaIL10R-R4.diamond_count|stool-HhaIL10R-R3.diamond_count|stool-Hh-R4.diamond_count|stool-Hh-R3.diamond_count|stool-WT-R4.diamond_count|stool-aIL10R-R1.diamond_count|
    +============+===============================+===============================+=========================+=========================+=========================+=============================+
    |NOG243840   |2                              |4                              |6                        |10                       |1                        |0                            |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |NOG281778   |1                              |5                              |4                        |4                        |2                        |1                            |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |NOG41625    |113                            |744                            |414                      |1273                     |404                      |567                          |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |bactNOG18808|1                              |2                              |8                        |8                        |15                       |21                           |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |COG3010     |2118                           |2395                           |1061                     |1738                     |2483                     |1043                         |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |proNOG56664 |7                              |2                              |3                        |2                        |4                        |5                            |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |COG3012     |15                             |60                             |34                       |101                      |38                       |22                           |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |COG3014     |27                             |18                             |32                       |92                       |41                       |17                           |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+
    |COG3015     |13                             |298                            |91                       |206                      |73                       |1148                         |
    +------------+-------------------------------+-------------------------------+-------------------------+-------------------------+-------------------------+-----------------------------+




Now we have count tables for genera and NOGs for both metagenomic and metatranscriptomic data we can start doing some analysis.

So to recap, in our working directories, DNA/ and RNA/ we now have count tables for alignments to genera and NOGs. Given the large 
sizes of raw and alignment files these have been deposited at the EBI ENA (ADD LINK). However for reproducing our downstream analysis
we have provided the count tables in the data/DNA/ and data/RNA/ directories (genus.diamond.aggregated.counts.tsv.gz and gene_counts.tsv.gz).
These can therefore be used for downstream analysis. 

If you have not run the above steps yourself then link to the counts tables we have provided and load into the csvdb database.

For RNA::

    $ cd RNA
    $ ln -s ../data/RNA/genus.diamond.aggregated.counts.tsv.gz .
    $ ln -s ../data/RNA/gene_counts.tsv.gz .
    $ zcat genus.diamond.aggregated.counts.tsv.gz | python cgat/scripts/csv2db.py
                                                    --backend=sqlite 
                                                    --retry                              
                                                    --table=genus_diamond_aggregated_counts
                                                    > genus.diamond.aggregated.counts.tsv.gz.load

    $ zcat gene_counts.tsv.gz | python cgat/scripts/csv2db.py
                                --backend=sqlite 
                                --retry                              
                                --table=gene_counts
                                > gene_counts.tsv.gz.load


For DNA::

    $ cd ../DNA
    $ ln -s ../data/DNA/genus.diamond.aggregated.counts.tsv.gz .
    $ ln -s ../data/DNA/gene_counts.tsv.gz .
    $ zcat genus.diamond.aggregated.counts.tsv.gz | python cgat/scripts/csv2db.py
                                                    --backend=sqlite 
                                                    --retry                              
                                                    --table=genus_diamond_aggregated_counts
                                                    > genus.diamond.aggregated.counts.tsv.gz.load


    $ zcat gene_counts.tsv.gz | python cgat/scripts/csv2db.py
                                --backend=sqlite 
                                --retry                              
                                --table=gene_counts
                                > gene_counts.tsv.gz.load


.. _fastx toolkit: http://hannonlab.cshl.edu/fastx_toolkit/






