# load modules
from ruffus import *

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV

import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random

import numpy as np
import sqlite3
# import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
# import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import pandas


###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
    ["pipeline.ini"])


PARAMS = P.PARAMS


###################################################
###################################################
###################################################


@follows(mkdir("associate_species.dir"))
@transform("alignment_genes.dir/*.genes.tsv.gz",
           regex("(\S+)/(\S+).genes.tsv.gz"),
           add_inputs("alignment_taxa.dir/*.lca.gz"),
           r"associate_species.dir/\2.species.counts.tsv.gz")
def buildCountTaxaInCogsRNA(infiles, outfile):
    '''
    reads that map to all NOGs are cross-referenced
    with their genera assignment - provides per taxa counts
    '''

    job_options="-l mem_free=20G"
    alignment_genes = infiles[0]
    track = P.snip(
        os.path.basename(alignment_genes), ".diamond.genes.tsv.gz")
    alignment_taxa = [
        x for x in infiles[1:] if os.path.basename(x).startswith(track)][0]

    statement = '''python ../../../src/scripts/diff2genera.py 
                   -m /ifs/projects/proj029/data/IGC/gene2cog.tsv.gz
                   -d gene_list.tsv
                   --alignment-taxa=%(alignment_taxa)s
                   --alignment-genes=%(alignment_genes)s
                   --counts
                   --level=species
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()



#########################################
#########################################
#########################################


@merge(buildCountTaxaInCogsRNA, "associate_species.dir/species_counts.tsv.gz")
def mergeSpeciesCounts(infiles, outfile):
    '''
    merge associated species counts
    '''
        # USE THE SAME GLOB AS IN THE COMBINING TABLES SCRIPT - maintain correct order
    prefixes = [
        P.snip(os.path.basename(x), ".species.counts.tsv.gz") 
        for x in glob.glob("associate_species.dir/*.species.counts.tsv.gz")]
    prefixes = ",".join(prefixes)

    statement = '''python %(scriptsdir)s/combine_tables.py
                   --missing=0
                   --columns=1,2
                   --take=preads
                   --glob=associate_species.dir/*.species.counts.tsv.gz
                   --prefixes=%(prefixes)s
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()


#########################################
#########################################
#########################################

@follows(mkdir("diff.dir"))
@transform(mergeSpeciesCounts, 
           regex("(\S+)/(\S+).tsv.gz"),
           r"diff.dir/\2.diff.tsv")
def runMetagenomeSeq(infile, outfile):
    '''
    run metagenomeSeq per species
    '''
    rscriptsdir = PARAMS.get("rscriptsdir")
    rscript = PARAMS.get("metagenomeseq_rscript")
    prefix = os.path.join(os.path.dirname(outfile),
                          P.snip(os.path.basename(infile), ".tsv.gz"))
    k = PARAMS.get("metagenomeseq_k")
    a = PARAMS.get("metagenomeseq_a")
    statement = '''%(rscript)s %(rscriptsdir)s/run_metagenomeseq.R
                   -c %(infile)s 
                   -p %(prefix)s
                   --k %(k)i 
                   --a %(a)f > %(outfile)s.log'''

    P.run()



#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
