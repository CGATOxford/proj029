import os, sys
import sqlite3
import gzip 

database = sys.argv[1]

dbh = sqlite3.connect(database)
cc = dbh.cursor()

######################
# essential genes data
######################
sys.stdout.write("contig\torf\thmm_id\tlength\tgc\tcontrol_cov\thh_cov\ttaxa\n")
statement = """select length.scaffold_name
               , gene.orf
               , gene.hmm_profile as hmm_id
               , length.length
               , gc.pGC
               , control.cov_mean as control_cov
               , hh.cov_mean as hh_cov
               , taxa.phylum 
               from idba_agg_agg_agg_filtered_contigs_lengths as length
               , agg_agg_agg_filtered_contigs_genes_essential_hmm_contigs as gene
               , idba_agg_agg_agg_filtered_contigs_gc as gc
               , idba_filteredmm10_gut_control_R1_filtered_contigs_coverage_stats as control
               , idba_filteredmm10_gut_hh_il10_R1_filtered_contigs_coverage_stats as hh
               , agg_agg_agg_filtered_contigs_taxa as taxa 
               where length.scaffold_name = gc.id 
               and control.contig = length.scaffold_name 
               and hh.contig = length.scaffold_name 
               and taxa.id = length.scaffold_name 
               and gene.contig = length.scaffold_name"""
for data in cc.execute(statement).fetchall():
    sys.stdout.write("\t".join(map(str, data)) + "\n")
