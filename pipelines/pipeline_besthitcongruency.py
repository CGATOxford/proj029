"""
=======================================================
Check congruency of aligned reads in terms of COG
annotations
=======================================================

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

import numpy
import sqlite3
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R

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

@follows(mkdir("congruency.dir"))
@transform("*.tsv.gz", 
           regex("(\S+).tsv.gz"),
           add_inputs(PARAMS.get("annotations_gene2cog")),
           r"congruency.dir/\1.congruency.tsv.gz")
def buildCongruency(infiles, outfile):
    '''
    build a table of reads with the % congruent with
    the best hit. We take the top n alignments as representing
    a random sample
    '''
    n = PARAMS.get("sample_sample")
    alignments, annotations = infiles
    temp = P.getTempFilename(".")
    statement = '''zcat %(alignments)s | cut -f1 | uniq | head -n %(n)i > %(temp)s;
                   python %(scriptsdir)s/best_hit_congruency.py
                   -a %(alignments)s
                   -m %(annotations)s
                   -r %(temp)s
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s;
                   rm -rf %(temp)s'''
    P.run()

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@transform(buildCongruency, suffix(".tsv.gz"), ".pdf")
def plotCongruency(infile, outfile):
    '''
    histogram congruency
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, se = "\t")''' % infile)
    R('''congruent_over_50 <- (length(dat$pcongruent[dat$pcongruent >= 50])/nrow(dat))*100''')

    R('''ggplot(dat, aes(x = pcongruent)) + geom_histogram() + labs(title = congruent_over_50)''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(plotCongruency)
def full():
    pass

#########################################
#########################################
#########################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
