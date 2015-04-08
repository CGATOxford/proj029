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
import CGATPipelines.Pipeline as P
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

@follows(mkdir("subsampled.dir"))
@transform("*.diamond.tsv.gz", regex("(\S+).tsv.gz"), r"subsampled.dir/\1.subsampled.tsv.gz")
def subsampleAlignments(infile, outfile):
    '''
    downsample data to look at properties between
    nanopore and illumina data on a reasonable number
    of alignments
    '''
    sample = PARAMS.get("alignments_sample")
    statement = '''python %(projectscriptsdir)s/diamond2sample.py
                   -a %(infile)s
                   -s %(sample)i
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s
                '''

###################################################
###################################################
###################################################
# look at coverage consistency for B.fragilis
###################################################
###################################################
###################################################

@follows(mkdir("bfragilis.dir"))
@merge(["agg-agg-agg.filtered.contigs.taxa.gz",
        "stool-HhaIL10R-R2.filtered.contigs.coverage.binned.taxa"],
       "bfragilis.dir/bfragilis_coverage.tsv.gz")
def buildCoverageOverBfragilis(infiles, outfile):
    '''
    build data set of coverage over b.fragilis 
    reads
    '''
    tmp = P.getTempFilename(".")
    statement = '''python %(cgatscriptsdir)s/combine_tables.py 
                   --glob=*taxa*
                   --log=%(outfile)s.log
                   | gzip > %(tmp)s;
                   zcat %(tmp)s 
                   | grep Bacteroides_fragilis 
                   | gzip > %(outfile)s;
                   rm -rf %(tmp)s'''
    P.run()
    
###################################################
###################################################
###################################################

@transform(buildCoverageOverBfragilis, suffix(".tsv.gz"), ".reads.tsv")
def buildBfragilisIdsNanopore(infile, outfile):
    '''
    get list of nanopore reads assigned to 
    B.fragilis
    '''
    statement = '''zcat %(infile)s | cut -f1 > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################

@follows(mkdir("fasta.dir"))
@transform("stool-nanopore-R1.fasta.gz", regex("(\S+).gz"), r"fasta.dir/\1.gz")
def formatFastaNanopore(infile, outfile):
    '''
    read names have lost compatibility between 
    files so need to format the original fasta
    file
    '''
    statement = '''zcat %(infile)s 
                   | sed -e 's/ [0-9].*//g' 
                   | gzip > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################

@transform(formatFastaNanopore,
           regex("(\S+)/(\S+)-(\S+)-R1.fasta.gz"),
           add_inputs(buildBfragilisIdsNanopore),
           r"bfragilis.dir/\2-\3bfragilis-R1.fasta.gz")
def buildReadsBfragilisNanopore(infiles, outfile):
    '''
    filter fastq file for reads mapping to 
    b fragilis
    '''
    fasta, ids = infiles
    statement = '''zcat %(fasta)s 
                  | python %(projectscriptsdir)s/fasta2filtered.py
                  -f %(ids)s
                  --log=%(outfile)s.log
                  | gzip > %(outfile)s
                '''
    P.run()

###################################################
###################################################
###################################################

attributes = ["stool-HhaIL10R-R2.filtered.contigs.coverage.stats.gz",
              "agg-agg-agg.filtered.contigs.lengths.tsv",
              "agg-agg-agg.filtered.contigs.taxa.gz",
              "agg-agg-agg.filtered.contigs.gc.tsv"]

@follows(mkdir("attributes.dir"))
@merge(attributes,
       "attributes.dir/contig_attributes.tsv")

def buildContigAttributes(infiles, outfile):
    '''
    merge attributes obtained by mapping
    illumina reads to nanopore contigs
    '''
    R('''read <- function(inf){ return(
                                    read.csv(inf, 
                                             header = T, 
                                             stringsAsFactors = F, 
                                             sep = "\t",
                                             row.names = 1)
                                       )
                               }'''
      )
    R('''cov <- read("%s")''' % infiles[0])
    R('''lengths <- read("%s")''' % infiles[1])
    R('''taxa <- read("%s")''' % infiles[2])
    R('''gc <- read("%s")''' % infiles[3])
    
    # gc and lengths have full read names others dont
    R('''rownames(lengths) <- unlist(strsplit(rownames(lengths), " "))[seq(1,nrow(lengths)*2,2)]''')
    R('''rownames(gc) <- unlist(strsplit(rownames(gc), " "))[seq(1,nrow(gc)*2,2)]''')
    
    # merge based on rownames - will introduce NAs presumably
    R('''dat <- data.frame(
                           cov[rownames(cov),]$cov_mean,
                           lengths[rownames(cov),1],
                           gc[rownames(cov),]$pGC,
                           taxa[rownames(cov),]$species,
                           taxa[rownames(cov),]$genus,
                           taxa[rownames(cov),]$family,
                           taxa[rownames(cov),]$X_order,
                           taxa[rownames(cov),]$class,
                           taxa[rownames(cov),]$phylum)''')
    R('''colnames(dat) <- c("coverage", 
                            "length", 
                            "gc", 
                            "species",
                            "genus",
                            "family",
                            "order",
                            "class",
                            "phylum")''')
    R('''dat$id <- rownames(cov)''')
    
    R('''write.table(dat, file = "%s", row.names = F, sep = "\t")''' % outfile)
    

#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
