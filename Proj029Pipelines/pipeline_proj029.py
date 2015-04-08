"""
=====================================================
Analysis for project 29
=====================================================

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
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
#from pandas import *
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

###################################################################
###################################################################
###################################################################
# This first section deals with collating the information from
# pipeline_metagenomeassembly.py. We produce plots of relative
# abundance correllations etc between different samples
###################################################################
###################################################################
###################################################################

@follows(mkdir("metaphlan.dir"))
@split(os.path.join(PARAMS.get("communities_dir"), PARAMS.get("communities_db")), 
       "metaphlan.dir/relab.*.matrix")
def buildRelativeAbundanceMatricesMetaphlan(infile, outfiles):
    '''
    build a matrix combining the relative abundance estimations
    for all samples
    '''
    # filenames to derive tablenames
    dirname = PARAMS.get("communities_dir")
    exp = "metaphlan.dir/*.relab"
    files = glob.glob(os.path.join(dirname, exp))
    tablenames = [
        os.path.basename(x).replace(".", "_").replace("-", "_") \
            for x in files]
    tablenames.sort()

    for level in ["phylum", "class", "order", "family", "genus", "species"]:
        outfile = os.path.join("communities.dir", "relab." + level + ".matrix")
        PipelineProj029.buildRelativeAbundanceMatrix(infile,
                                                     tablenames,
                                                     outfile,
                                                     level = level)

###################################################################
###################################################################
###################################################################

@follows(mkdir("kraken.dir"))
@split(os.path.join(PARAMS.get("communities_dir"), PARAMS.get("communities_db")), 
       "kraken.dir/*norm*.matrix")
def buildAbundanceMatricesKraken(infile, outfiles):
    '''
    build a matrix combining the rpm estimations
    for all samples
    '''
    # filenames to derive tablenames
    dirname = PARAMS.get("communities_dir")
    exp = "kraken.dir/*.counts.norm.tsv.gz"
    files = glob.glob(os.path.join(dirname, exp))
    tablenames = [
        "kraken_" + os.path.basename(P.snip(x, ".tsv.gz")).replace(".", "_").replace("-", "_") \
            for x in files]
    tablenames.sort()

    for level in ["phylum", "class", "order", "family", "genus", "species"]:
        outfile = os.path.join("kraken.dir", "counts.norm." + level + ".matrix")
        PipelineProj029.buildRelativeAbundanceMatrix(infile,
                                                     tablenames,
                                                     outfile,
                                                     level = level)

        
###################################################################
###################################################################
###################################################################

@follows(mkdir("diamond.dir"))
@split(os.path.join(PARAMS.get("communities_dir"), PARAMS.get("communities_db")), 
       "diamond.dir/*norm*.matrix")
def buildAbundanceMatricesDiamond(infile, outfiles):
    '''
    build a matrix combining the rpm estimations
    for all samples
    '''
    # filenames to derive tablenames
    dirname = PARAMS.get("communities_dir")
    exp = "diamond.dir/*.taxa.count"
    files = glob.glob(os.path.join(dirname, exp))
    tablenames = [
        os.path.basename(x).replace(".", "_").replace("-", "_") \
            for x in files]
    tablenames.sort()

    for level in ["phylum", "class", "order", "family", "genus", "species"]:
        outfile = os.path.join("diamond.dir", "counts.norm." + level + ".matrix")
        PipelineProj029.buildRelativeAbundanceMatrix(infile,
                                                     tablenames,
                                                     outfile,
                                                     level = level)


###################################################################
###################################################################
###################################################################
COMMUNITIES_TARGETS = []
communities_targets = {"kraken": buildAbundanceMatricesKraken,
                       "metaphlan": buildRelativeAbundanceMatricesMetaphlan,
                       "diamond": buildAbundanceMatricesDiamond}
for x in P.asList(PARAMS.get("classifiers")):
    COMMUNITIES_TARGETS.append(communities_targets[x])

@transform(COMMUNITIES_TARGETS,
           suffix(".matrix"), ".barplot.pdf")
def barplotAbundances(infile, outfile):
    '''
    barplot the species relative abundances
    '''
    threshold = PARAMS.get("communities_threshold")
    PipelineProj029.barplotAbundances(infile,
                                      outfile,
                                      threshold)

###################################################################
###################################################################
###################################################################

@transform(COMMUNITIES_TARGETS, suffix(".matrix"), ".ratio.tsv")
def calculateFirmicutesBacteroidetesRatio(infile, outfile):
    '''
    barplot the species relative abundances
    '''
    threshold = PARAMS.get("communities_threshold")
    PipelineProj029.calculateFirmicutesBacteroidetesRatio(infile,
                                                          outfile,
                                                          threshold)

###################################################################
###################################################################
###################################################################

@transform(calculateFirmicutesBacteroidetesRatio, suffix(".tsv"), ".pdf")
def plotFirmicutesBacteroidetesRatio(infile, outfile):
    '''
    produce boxplot of firmicutes/bacteroidetes ratio
    '''
    PipelineProj029.plotFirmicutesBacteroidetesRatio(infile,
                                                     outfile)
    
###################################################################
###################################################################
###################################################################
@transform(calculateFirmicutesBacteroidetesRatio, suffix(".tsv"), ".signif")
def calculateSignificanceOfFirmicutesBacteroidetesRatio(infile, outfile):
    '''
    use tuleyHSD to calculate significance between groups with
    multiple testing correction
    '''
    PipelineProj029.calculateSignificanceOfFirmicutesBacteroidetesRatio(infile, 
                                                                        outfile)


###################################################################
###################################################################
###################################################################

@jobs_limit(1, "R")
@transform(COMMUNITIES_TARGETS, suffix(".matrix"), ".barplot.numbers.pdf")
def plotHowManySpecies(infile, outfile):
    '''
    how many samples have how many species?
    '''
    PipelineProj029.plotHowManySpecies(infile, outfile)

###################################################################
###################################################################
###################################################################

@transform(COMMUNITIES_TARGETS, suffix(".matrix"), ".heatmap.pdf")
def heatmapAbundances(infile, outfile):
    '''
    heatmap the species relative abundances
    '''

    threshold = PARAMS.get("communities_threshold")
    PipelineProj029.heatmapAbundances(infile,
                                      outfile,
                                      threshold,
                                      "covariates.tsv")

###################################################################
###################################################################
###################################################################

@transform(COMMUNITIES_TARGETS, suffix(".matrix"), ".signif")
def testSignificanceOfAbundances(infile, outfile):
    '''
    use an anova to test significance. This is not ideal but serves
    as a quick look
    '''
    PipelineProj029.anovaTest(infile,
                              "covariates.tsv", 
                              outfile,
                              threshold = PARAMS.get("communities_threshold"),
                              over = PARAMS.get("communities_over"))

###################################################################
###################################################################
###################################################################

@transform(testSignificanceOfAbundances, 
           suffix(".signif"), 
           add_inputs(COMMUNITIES_TARGETS),
           ".signif.pdf")
def plotSignificantResults(infiles, outfile):
    '''
    barplot those taxa that are different across groups
    '''
    inf = infiles[0]
    track = P.snip(inf, ".signif")

    abundance_file = [m for m in infiles[1] if m.find(track) != -1][0]
    threshold = PARAMS.get("communities_threshold")
    PipelineProj029.plotSignificantResults(inf, abundance_file, outfile, threshold)

###################################################################
###################################################################
###################################################################
# KEGG analysis
###################################################################
###################################################################
###################################################################

@follows(mkdir("kegg.dir"))
@merge(glob.glob(os.path.join(PARAMS.get("communities_dir"), "kegg.dir/*.kegg.counts")), "kegg.dir/kegg.pathways.matrix")
def combineKeggTables(infiles, outfile):
    '''
    merge counts for kegg pathways
    '''
    headers = ",".join(
        [re.match(".*.dir/(.*).kegg.counts", x).groups()[0]                                                        
         for x in infiles])
    directory = os.path.dirname(infiles[0])

    statement = '''python %(scriptsdir)s/combine_tables.py
                   --glob=%(directory)s/*.kegg.counts
                   --headers=%(headers)s
                   --columns=1
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################
@follows(barplotAbundances,
         plotFirmicutesBacteroidetesRatio,
         calculateSignificanceOfFirmicutesBacteroidetesRatio,
         plotHowManySpecies,
         heatmapAbundances,
         plotSignificantResults)
def full():
    pass

#########################################
#########################################
#########################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
