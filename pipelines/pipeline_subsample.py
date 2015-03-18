"""
=======================================
Compare RNA and DNA sequencing results
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

###################################################
###################################################
###################################################

@transform("*.counts.tsv.gz", suffix(".tsv.gz"), ".subsampled.tsv.gz")
def subsample(infile, outfile):
    '''
    subsample based on parameters in pipeline.ini
    '''
    sample = PARAMS.get("sample")
    statement = '''zcat %(infile)s 
                   | python %(scriptsdir)s/counts2subsample.py
                   -s %(sample)s
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@merge(subsample, "agg-agg-agg.counts.subsampled.tsv.gz")
def aggregateCounts(infile, outfile):
    '''
    aggregate the susampled data together
    '''
    # USE THE SAME GLOB AS IN THE COMBINING TABLES SCRIPT - maintain correct order
    prefixes = [
        P.snip(os.path.basename(x), ".genes.counts.subsampled.tsv.gz") for x in glob.glob("*.genes.counts.subsampled.tsv.gz")]
    prefixes = ",".join(prefixes)

    statement = '''python %(cgatscriptsdir)s/combine_tables.py
                   --missing=0
                   --columns=1
                   --take=count
                   --glob=*.genes.counts.subsampled.tsv.gz
                   --prefixes=%(prefixes)s
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@transform(aggregateCounts,
           suffix(".tsv.gz"), 
           ".diff.tsv")
def runMetagenomeSeq(infile, outfile):
    '''
    run metagenomeSeq - a tool for calculating significance
    based on gene counts
    '''
    rscriptsdir = PARAMS.get("rscriptsdir")
    rscript = PARAMS.get("metagenomeseq_rscript")
    prefix = P.snip(infile, ".tsv.gz")

    a = PARAMS.get("metagenomeseq_a")
    k = PARAMS.get("metagenomeseq_k")

    statement = '''%(rscript)s %(rscriptsdir)s/run_metagenomeseq.R
                   -c %(infile)s 
                   -p %(prefix)s
                   --k %(k)i 
                   --a %(a)i > %(outfile)s.log'''

    P.run()

#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
