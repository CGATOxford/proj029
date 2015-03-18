"""
=======================================
defining a set of housekeeping genes
=======================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

"""

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
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import pandas
import PipelineProj029

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["pipeline.ini"])


PARAMS = P.PARAMS

###################################################################
# connecting to database
###################################################################

def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''
    dbh = sqlite3.connect(PARAMS["database"])
    return dbh

###################################################
###################################################
###################################################

# global variables
COG = PARAMS.get("cog")

###################################################
###################################################
###################################################

@transform(PARAMS.get("igc_file"), 
           regex("(\S+)/(\S+)"),
           r"gene_names.tsv.gz")
def getGenesAssociatedWithCOG(infile, outfile):
    '''
    get all genes from the IGC that are associated with
    a COG
    '''
    cog = PARAMS.get("cog")
    gene2cog = PARAMS.get("gene2cog")

    statement = '''python %(scriptsdir)s/cog2fasta.py 
                  -m %(gene2cog)s
                  -c %(cog)s 
                  -f %(infile)s 
                  --output-list 
                  --log=%(outfile)s.log'''
    P.run()

#########################################
#########################################
#########################################

@follows(mkdir("fastq.dir"), getGenesAssociatedWithCOG)
@transform("*.fastq.gz", 
           regex("(\S+).fastq.gz"), 
           r"fastq.dir/%s_\1.list.gz" % COG)
def buildReadList(infile, outfile):
    '''
    build list of fastq reads mapping to COG
    '''
    afile = P.snip(infile, ".fastq.gz") + ".diamond.genes.tsv.gz"
    alignment = os.path.join(PARAMS.get("rna_communities_dir"), afile)
    statement = '''python %(scriptsdir)s/genes2reads.py
                   -g gene_names.tsv.gz
                   -a %(alignment)s 
                    --log=%(outfile)s.log
                    | gzip
                   > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@transform("*.fastq.gz", 
           regex("(\S+).fastq.gz"), 
           add_inputs(buildReadList),
           r"fastq.dir/%s_\1.fastq.gz" % COG)
def buildFastq(infiles, outfile):
    '''
    filter fastq files based on reads
    mapping to COG
    '''
    fastq = infiles[0]
    reads = [x for x in infiles[1:] if P.snip(x, ".list.gz") in outfile][0]

    statement = '''zcat %(fastq)s 
                   | python %(cgat_scriptsdir)s/fastq2fastq.py 
                   --apply=%(reads)s
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
