"""
============================================
Compare metagenomic and metatranscriptomic
data sets
============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
========

This pipeline is an extended analysis after running
pipeline_metagenomecommunities.py on both the RNA-seq
and DNA-seq data sets. Its purpose is to compare the 
results between data sets and to combine them in a
meaningful way in order to define candidate genera and
functions of the microbiome that are to be taken forward
into downstream experiments.

The pipeline is a wrapper in order to parallelise tasks
using various scripts and functions defined in modules/,
scripts/ and R/ directories.

configuration
==============

Parameters for the pipeline are specified in a pipeline.ini
file that should be present in the working directory.

code
====

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
# import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
# import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import pandas
import PipelineMetaomics


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

##########################
# connecting to database
##########################

def connect():
    '''
    connect to database.
    '''
    dbh = sqlite3.connect(PARAMS["database"])
    return dbh

###################################################
###################################################
###################################################
# SECTION 1: Build stats around differentially
# abundant genera and NOGs in metagenomic and
# metatranscriptomic analyses.
###################################################
###################################################
###################################################

@follows(mkdir("diff_stats.dir"))
@transform([os.path.join(
        PARAMS.get("rna_communitiesdir"), 
        "counts.dir/genus.diamond.aggregated.counts.diff.load"),
           os.path.join(
        PARAMS.get("dna_communitiesdir"), 
        "counts.dir/genus.diamond.aggregated.counts.diff.load")],
           regex("(\S+)/(.*NA)/(\S+)/genus.diamond.aggregated.counts.diff.load"),
           r"diff_stats.dir/\2.diff.stats.tsv")
def buildGeneraDiffStats(infile, outfile):
    '''
    build stats on the number of genera that are differentially
    abundant at different levels of significance and fold changes
    '''
    dbh = connect()
    cc = dbh.cursor()
    db = os.path.join(P.snip(os.path.dirname(infile), "counts.dir"), "csvdb")
    PipelineMetaomics.buildDiffStats(infile,
                                     outfile,
                                     db,
                                     cc)

#########################################
#########################################
#########################################

@follows(mkdir("diff_stats.dir"))
@transform([os.path.join(
        PARAMS.get("rna_communitiesdir"), 
        "genes.dir/gene_counts.diff.load"),
           os.path.join(
        PARAMS.get("dna_communitiesdir"), 
        "genes.dir/gene_counts.diff.load")],
           regex("(\S+)/(.*NA)/(\S+)/gene_counts.diff.load"),
           r"diff_stats.dir/\2.genes.diff.stats.tsv")
def buildGeneDiffStats(infile, outfile):
    '''
    build stats on the number of NOGs that are differentially
    abundant at different levels of significance and fold changes
    '''
    dbh = connect()
    cc = dbh.cursor()
    db = os.path.join(P.snip(os.path.dirname(infile), "genes.dir"), "csvdb")
    PipelineMetaomics.buildDiffStats(infile, 
                                     outfile,
                                     db,
                                     cc)

#########################################
#########################################
#########################################

@follows(buildGeneraDiffStats, 
         buildGeneDiffStats)
def diff_stats():
    pass

###################################################
###################################################
###################################################
# SECTION 2: Build lists of overlapping sets
# of genera/NOGs - total and differentially
# abundant. Also scatterplot abundance estimates
# between data sets
###################################################
###################################################
###################################################

@follows(mkdir("diff.dir"))
@merge([os.path.join(
            PARAMS.get("rna_communitiesdir"), "csvdb"),
        os.path.join(
            PARAMS.get("dna_communitiesdir"), "csvdb")],
       "diff.dir/common_genes.tsv")
def buildCommonGeneList(infiles, outfile):
    '''
    get a list of NOGs that were common beteween DNA
    and RNA analysis - used for downstream analysis
    '''
    rnadb, dnadb = infiles
    PipelineMetaomics.buildCommonList(rnadb,
                                      dnadb,
                                      outfile)

#########################################
#########################################
#########################################

@follows(mkdir("diff.dir"))
@merge([os.path.join(
            PARAMS.get("rna_communitiesdir"), "csvdb"),
        os.path.join(
            PARAMS.get("dna_communitiesdir"), "csvdb")],
           "diff.dir/common_genera.tsv")
def buildCommonGeneraList(infiles, outfile):
    '''
    get a list of genera that were common beteween DNA
    and RNA analysis - used for downstream analysis
    '''
    rnadb, dnadb = infiles
    PipelineMetaomics.buildCommonList(rnadb,
                                      dnadb,
                                      outfile)

#########################################
#########################################
#########################################

@follows(mkdir("diff.dir"))
@transform([os.path.join(PARAMS.get("rna_communitiesdir"), "csvdb"),
            os.path.join(PARAMS.get("dna_communitiesdir"), "csvdb")],
           regex("(\S+)/(.*NA).*/csvdb"),
           add_inputs(buildCommonGeneList),
           r"diff.dir/\2_HhaIL10R_vs_WT.diff.genes.tsv")
def buildGeneDiffList(infiles, outfile):
    '''
    build a list of differentially expressed cogs
    between HhaIL10R and WT
    '''
    db, common = infiles
    PipelineMetaomics.buildGeneDiffList(db,
                                        common,
                                        outfile)

#########################################
#########################################
#########################################

@merge(buildGeneDiffList, "diff.dir/genes_overlap.tsv")
def buildDiffGenesOverlap(infiles, outfile):
    '''
    overlap gene lists
    '''
    dnafile, rnafile = infiles
    PipelineMetaomics.buildDiffGeneOverlap(dnafile,
                                           rnafile,
                                           outfile)

#########################################
#########################################
#########################################

@merge([buildCommonGeneList, buildDiffGenesOverlap],
       "diff.dir/genes_overlap.sig")
def testSignificanceOfGenesOverlap(infiles, outfile):
    '''
    test significance of overlapping diff gene
    lists bewteen RNA and DNA using hypergeometric
    test
    '''
    common, overlap = infiles
    PipelineMetaomics.testSignificanceOfOverlap(common,
                                                overlap,
                                                outfile)

###################################################
###################################################
###################################################

@follows(mkdir("compare_detected.dir"))
@split([os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "counts.dir/*diamond*.aggregated.counts.tsv.gz"),
        os.path.join(
            PARAMS.get("dna_communitiesdir"), 
            "counts.dir/*diamond*.aggregated.counts.tsv.gz")],
       "compare_detected.dir/*.overlap.tsv")
def buildTaxaDetectionOverlap(infiles, outfiles):
    '''
    build taxa detection overlap
    '''
    levels = ["phylum", "class", "order", "family", "genus", "species"]
    for level in levels:
        rnacounts, dnacounts = [inf for inf in infiles if inf.find(level) != -1]
        outfile = "compare_detected.dir/%s.overlap.tsv" % level
        PipelineMetaomics.buildDetectionOverlap(rnacounts,
                                                dnacounts,
                                                outfile)

###################################################
###################################################
###################################################
@follows(mkdir("compare_detected.dir"))
@split([os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "genes.dir/gene_counts.tsv.gz"),
       os.path.join(
        PARAMS.get("dna_communitiesdir"), 
        "genes.dir/gene_counts.tsv.gz")],
       "compare_detected.dir/genes.overlap.tsv")
def buildGeneDetectionOverlap(infiles, outfile):
    '''
    build NOG detection overlap
    '''
    rnacounts, dnacounts = infiles
    PipelineMetaomics.buildDetectionOverlap(rnacounts, 
                                            dnacounts,
                                            outfile)


###################################################
###################################################
###################################################

@follows(mkdir("compare_abundance.dir"))
@split([os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "counts.dir/*diamond*.norm.matrix"),
        os.path.join(
            PARAMS.get("dna_communitiesdir"), 
            "counts.dir/*diamond*.norm.matrix")],
       "compare_abundance.dir/*.pdf")
def scatterplotTaxaAbundanceEstimates(infiles, outfiles):
    '''
    scatterplot abundance estimates for each sample
    '''
    levels = ["phylum", "class", "order", "family", "genus", "species"]
    for level in levels:
        print level
        outfile = "compare_abundance.dir/average_abundance.%s.png" % level
        rnamatrix, dnamatrix = [inf for inf in infiles if inf.find(level) != -1]
        PipelineMetaomics.scatterplotAbundanceEstimates(dnamatrix, 
                                                        rnamatrix, 
                                                        outfile)

###################################################
###################################################
###################################################

@follows(mkdir("compare_abundance.dir"))
@split([os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "genes.dir/*.norm.matrix"),
        os.path.join(
            PARAMS.get("dna_communitiesdir"), 
            "genes.dir/*.norm.matrix")],
       "compare_abundance.dir/average_abundance_genes.png")
def scatterplotGeneAbundanceEstimates(infiles, outfile):
    '''
    scatterplot abundance estimates for each sample
    '''
    rnamatrix, dnamatrix = infiles
    PipelineMetaomics.scatterplotAbundanceEstimates(dnamatrix, 
                                                    rnamatrix, 
                                                    outfile)

###################################################
###################################################
###################################################

@follows(mkdir("compare_detected.dir"))
@split([os.path.join(
        PARAMS.get("rna_communitiesdir"), 
        "counts.dir/*diamond*.aggregated.counts.tsv.gz"), 
       os.path.join(
        PARAMS.get("dna_communitiesdir"), 
        "counts.dir/*diamond*.aggregated.counts.tsv.gz")],
       "compare_detected.dir/*.abundance.pdf")
def plotAbundanceLevelsOfTaxaOverlap(infiles, outfiles):
    '''
    plot the RPM values of taxa that do and don't
    overlap between data sets
    '''
    rnacounts, dnacounts = [inf for inf in infiles if inf.find("genus") != -1]
    outfile = "compare_detected.dir/genus.overlap.abundance.pdf"

    PipelineMetaomics.plotAbundanceLevelsOfOverlap(rnacounts, 
                                                   dnacounts,
                                                   outfile,
                                                   of = "genus")

###################################################
###################################################
###################################################

@follows(mkdir("compare_detected.dir"))
@merge([os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "genes.dir/gene_counts.tsv.gz"),
       os.path.join(
        PARAMS.get("dna_communitiesdir"), 
        "genes.dir/gene_counts.tsv.gz")],
       "compare_detected.dir/genes.abundance.pdf")
def plotAbundanceLevelsOfGeneOverlap(infiles, outfile):
    '''
    plot abundances for unique and common genes
    '''
    rnacounts, dnacounts = infiles
    PipelineMetaomics.plotAbundanceLevelsOfOverlap(rnacounts,
                                                   dnacounts,
                                                   outfile,
                                                   of = "genes")

###################################################
###################################################
###################################################

@follows(buildCommonGeneList,
         buildCommonGeneraList,
         buildGeneDiffList,
         testSignificanceOfGenesOverlap,
         buildTaxaDetectionOverlap,
         buildGeneDetectionOverlap,
         scatterplotTaxaAbundanceEstimates,
         scatterplotGeneAbundanceEstimates,
         plotAbundanceLevelsOfTaxaOverlap,
         plotAbundanceLevelsOfGeneOverlap)
def overlaps():
    pass

###################################################
###################################################
###################################################
# SECTION 3: PCA of normalised counts for 
# metagenomics and metatranscriptomics
###################################################
###################################################
###################################################

@follows(mkdir("pca.dir"))
@jobs_limit(1, "R")
@transform([os.path.join(
            PARAMS.get("dna_communitiesdir"), 
            "genes.dir/gene_counts.norm.matrix"),
            os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "genes.dir/gene_counts.norm.matrix"),
            os.path.join(
            PARAMS.get("dna_communitiesdir"), 
            "counts.dir/genus.diamond.aggregated.counts.norm.matrix"),
            os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "counts.dir/genus.diamond.aggregated.counts.norm.matrix")],
           regex("(\S+)/(\S+).matrix"),
           r"pca.dir/\2.loadings.tsv")
def runPCA(infile, outfile):
    '''
    run PCA analysis - output the plot
    and the loadings file
    '''
    PipelineMetaomics.runPCA(infile,
                             outfile)

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@follows(runPCA)
@transform("pca.dir/*na_*loadings.tsv", suffix(".tsv"), ".pdf")
def plotPCALoadings(infile, outfile):
    '''
    plot the loadings associated with differentially
    expressed genes
    '''
    PipelineMetaomics.plotPCALoadings(infile,
                                      outfile)

#########################################
#########################################
#########################################

@follows(plotPCALoadings)
def PCA():
    pass

###################################################
###################################################
###################################################
# SECTION 4: Identify candidate responsive NOGs
# by comparing metagenomic and metatranscriptomic
# data sets
###################################################
###################################################
###################################################

@follows(mkdir("rna_dna_ratio.dir"))
@merge([os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "genes.dir/gene_counts.diff.tsv"),
        os.path.join(
            PARAMS.get("dna_communitiesdir"), 
            "genes.dir/gene_counts.diff.tsv")],
       "rna_dna_ratio.dir/ratio_genes.tsv")
def buildRNADNARatio(infiles, outfile):
    '''
    build the ratio between RNA fold changes
    and DNA fold changes
    '''
    rnadiff, dnadiff = infiles
    PipelineMetaomics.buildRNADNARatio(dnadiff,
                                       rnadiff,
                                       outfile)

###################################################
###################################################
###################################################

@transform(buildRNADNARatio, 
           suffix(".tsv"), 
           add_inputs(buildGeneDiffList),
           ".annotated.tsv")
def annotateRNADNARatio(infiles, outfile):
    '''
    annotate the ratio data with how genes look in terms
    of DNA and RNA regulation
    '''
    infile, dna, rna = infiles
    PipelineMetaomics.annotateRNADNARatio(infile, 
                                          dna,
                                          rna,
                                          outfile)

###################################################
###################################################
###################################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".pdf")
def plotSets(infile, outfile):
    '''
    plot the fold changes in RNA and DNA analyses
    and label by how they are regulated in DNA and
    RNA analyses
    MUST HAVE GOI FILE IN WORKING DIR - not ideal
    '''
    PipelineMetaomics.plotSets(infile,
                               outfile)

###################################################
###################################################
###################################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".outsidepi.tsv")
def buildGenesOutsidePredictionInterval(infile, outfile):
    '''
    annotate genes as being outside prediction
    interval - these are the NOGs that we are
    defining as colitis-responsive
    '''
    PipelineMetaomics.buildGenesOutsidePredictionInterval(infile,
                                                          outfile)

#########################################
#########################################
#########################################

@follows(plotSets,
         buildGenesOutsidePredictionInterval)
def define_colitis_responsive_nogs():
    pass

###################################################
###################################################
###################################################
# SECTION 5: Associate COGs with the genera that
# express them - select reads that map to both
# cog and taxa
###################################################
###################################################
###################################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".list")
def buildGeneListForTaxaAssociations(infile, outfile):
    '''
    build a list of COGs for use with annotating to
    taxa
    '''
    statement = '''cat %(infile)s | 
                   cut -f1 | 
                   tail -n+2 
                   > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@follows(mkdir("associate_taxa.dir"), 
         buildGeneListForTaxaAssociations)
@transform(os.path.join(
        PARAMS.get("rna_communitiesdir"), 
        "genes.dir/*.genes.tsv.gz"),
           regex("(\S+)/(\S+).genes.tsv.gz"),
           add_inputs(os.path.join(
            PARAMS.get(
                "rna_communitiesdir"), 
            "diamond.dir/*.diamond.lca.gz")),
           r"associate_taxa.dir/\2.ctaxa.tsv.gz")
def buildCountTaxaInCogsRNA(infiles, outfile):
    '''
    reads that map to all NOGs are cross-referenced
    with their genera assignment - provides per taxa counts
    '''
    job_options="-l mem_free=20G"
    m = PARAMS.get("genes_map")
    alignment_genes = infiles[0]
    track = P.snip(
        os.path.basename(alignment_genes), ".diamond.genes.tsv.gz")
    alignment_taxa = [
        x for x in infiles[1] if os.path.basename(x).startswith(track)][0]

    statement = '''python %(projscripts)s/diff2genera.py 
                   -m %(m)s
                   -d rna_dna_ratio.dir/ratio_genes.annotated.list
                   --alignment-taxa=%(alignment_taxa)s
                   --alignment-genes=%(alignment_genes)s
                   --counts
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@merge(buildCountTaxaInCogsRNA, 
       "associate_taxa.dir/associated_taxa_counts.tsv.gz")
def mergeCountTaxaInCogsRNA(infiles, outfile):
    '''
    merge tables of taxa associated with genes
    '''
    pattern = os.path.dirname(infiles[0]) + "/*.ctaxa.tsv.gz"
    prefixes = ",".join(
        [P.snip(os.path.basename(x), ".diamond.ctaxa.tsv.gz") for x in glob.glob(pattern)])
    statement = '''python %(scriptsdir)s/combine_tables.py
                   --glob=%(pattern)s
                   --missing=0
                   --prefixes=%(prefixes)s
                   --columns=1,2
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@follows(mkdir("associate_taxa.dir"), 
         buildGeneListForTaxaAssociations)
@transform(os.path.join(
        PARAMS.get("rna_communitiesdir"), 
        "genes.dir/*.genes.tsv.gz"),
           regex("(\S+)/(\S+).genes.tsv.gz"),
           add_inputs(os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "diamond.dir/*.diamond.lca.gz")),
           r"associate_taxa.dir/\2.ptaxa.tsv.gz")
def buildProportionTaxaInCogsRNA(infiles, outfile):
    '''
    reads that map to differentially expressed genes are cross-referenced
    with their genera assignment
    '''
    job_options="-l mem_free=25G"
    m = PARAMS.get("genes_map")
    alignment_genes = infiles[0]
    track = P.snip(os.path.basename(alignment_genes), ".diamond.genes.tsv.gz")
    alignment_taxa = [x for x in infiles[1] if os.path.basename(x).startswith(track)][0]

    statement = '''python %(projscripts)s/diff2genera.py 
                   -m %(m)s
                   -d rna_dna_ratio.dir/ratio_genes.annotated.list
                   --alignment-taxa=%(alignment_taxa)s
                   --alignment-genes=%(alignment_genes)s
                   --level=genus
                   --log=%(outfile)s.log
                  | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@merge(buildProportionTaxaInCogsRNA, 
       "associate_taxa.dir/associated_ptaxa.tsv.gz")
def mergeProportionTaxaInCogsRNA(infiles, outfile):
    '''
    merge tables of taxa associated with genes
    '''
    temp = P.getTempFilename(".")
    pattern = os.path.dirname(infiles[0]) + "/*.ptaxa.tsv.gz"
    prefixes = ",".join([P.snip(os.path.basename(x), ".diamond.ptaxa.tsv.gz") for x in glob.glob(pattern)])
    statement = '''python %(scriptsdir)s/combine_tables.py
                   --glob=%(pattern)s
                   --missing=0
                   --prefixes=%(prefixes)s
                   --columns=1,2
                   --log=%(outfile)s.log
                    > %(temp)s'''
    P.run()

    inf = IOTools.openFile(temp)
    header = inf.readline()[:-1].split("\t")
    header = header[0].split("-") + header[1:]
    outf = IOTools.openFile(outfile, "w")
    outf.write("\t".join(header) + "\n")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        cog_taxa = data[0].split("-")
        new_data = cog_taxa + data[1:]
        outf.write("\t".join(new_data) + "\n")
    outf.close()
    os.unlink(temp)

#########################################
#########################################
#########################################

@follows(mkdir("associate_taxa.dir"), 
         buildGeneListForTaxaAssociations)
@transform(os.path.join(
        PARAMS.get("dna_communitiesdir"), 
        "genes.dir/*.genes.tsv.gz"),
           regex("(\S+)/(\S+).genes.tsv.gz"),
           add_inputs(os.path.join(
            PARAMS.get("dna_communitiesdir"), 
            "diamond.dir/*.diamond.lca.gz")),
           r"associate_taxa.dir/\2.dna_ctaxa.tsv.gz")
def buildCountTaxaInCogsDNA(infiles, outfile):
    '''
    reads that map to all NOGs are cross-referenced
    with their genera assignment - provides per taxa counts
    '''
    job_options="-l mem_free=20G"
    m = PARAMS.get("genes_map")
    alignment_genes = infiles[0]
    track = P.snip(
        os.path.basename(alignment_genes), ".diamond.genes.tsv.gz")
    alignment_taxa = [
        x for x in infiles[1] if os.path.basename(x).startswith(track)][0]

    statement = '''python %(projscripts)s/diff2genera.py 
                   -m %(m)s
                   -d rna_dna_ratio.dir/ratio_genes.annotated.list
                   --alignment-taxa=%(alignment_taxa)s
                   --alignment-genes=%(alignment_genes)s
                   --counts
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@merge(buildCountTaxaInCogsDNA, 
       "associate_taxa.dir/associated_taxa_counts_dna.tsv.gz")
def mergeCountTaxaInCogsDNA(infiles, outfile):
    '''
    merge tables of taxa associated with genes
    '''
    pattern = os.path.dirname(infiles[0]) + "/*.dna_ctaxa.tsv.gz"
    prefixes = ",".join(
        [P.snip(os.path.basename(x), ".diamond.dna_ctaxa.tsv.gz") 
         for x in glob.glob(pattern)])
    statement = '''python %(scriptsdir)s/combine_tables.py
                   --glob=%(pattern)s
                   --missing=0
                   --prefixes=%(prefixes)s
                   --columns=1,2
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@follows(mergeCountTaxaInCogsRNA,
         mergeProportionTaxaInCogsRNA,
         mergeCountTaxaInCogsDNA)
def associate_taxa():
    pass

###################################################
###################################################
###################################################
# SECTION 6: Plot those genera that contribute
# predominantly to NOG expression and run
# metagenomeSeq to get per taxa cog differences
# and fold changes
###################################################
###################################################
###################################################

@transform(mergeProportionTaxaInCogsRNA, 
           suffix(".tsv.gz"), 
           ".average_ptaxa.matrix")
def buildGenusCogCountsMatrix(infile, outfile):
    '''
    build cog x genus proportion
    matrix
    '''
    PipelineMetaomics.buildGenusCogCountsMatrix(infile,
                                                outfile)

###################################################
###################################################
###################################################

@merge([buildGenusCogCountsMatrix, 
        buildGenesOutsidePredictionInterval], 
       "associate_taxa.dir/taxa_cogs_matrix_annotated.pdf")
def plotMaxTaxaContribution(infiles, outfile):
    '''
    plot the distribution of maximum genus
    contribution per gene set
    '''
    matrix, annotations = infiles
    PipelineMetaomics.plotMaxTaxaContribution(matrix,
                                              annotations,
                                              outfile)

###################################################
###################################################
###################################################

@merge([buildGenusCogCountsMatrix, 
        buildGenesOutsidePredictionInterval], 
       "associate_taxa.dir/taxa_cogs_matrix_annotated.sig")
def testSignificanceOfMaxTaxaContribution(infiles, outfile):
    '''
    Test significance of distribution differences. Compared to NS
    group
    '''
    matrix, annotations = infiles
    PipelineMetaomics.testSignificanceOfMaxTaxaContribution(matrix,
                                                            annotations,
                                                            outfile)

###################################################
###################################################
###################################################

@transform(buildGenusCogCountsMatrix, 
           suffix(".matrix"), 
           add_inputs(buildGenesOutsidePredictionInterval),
           ".pdf")
def heatmapTaxaCogProportionMatrix(infiles, outfile):
    '''
    plot the taxa associated with each cog on
    a heatmap
    '''
    matrix, annotations = infiles
    PipelineMetaomics.heatmapTaxaCogProportionMatrix(matrix,
                                                     annotations,
                                                     outfile)

###################################################
###################################################
###################################################

@follows(mkdir("taxa_cogs_diff.dir"))
@transform([mergeCountTaxaInCogsRNA,
            mergeCountTaxaInCogsDNA],
            regex("(\S+)/(\S+).tsv.gz"), r"taxa_cogs_diff.dir/\2.diff.tsv")
def runMetagenomeSeqPerCOGAndTaxa(infile, outfile):
    '''
    run metagenomeSeq on taxa-COG counts
    '''
    rscriptsdir = PARAMS.get("rscriptsdir")
    rscript = PARAMS.get("metagenomeseq_rscript")
    prefix = P.snip(outfile, ".diff.tsv")

    k = PARAMS.get("metagenomeseq_genes_k")
    a = PARAMS.get("metagenomeseq_genes_a")

    statement = '''%(rscript)s %(rscriptsdir)s/run_metagenomeseq.R
                   -c %(infile)s 
                   -p %(prefix)s
                   --k %(k)i 
                   --a %(a)i > %(outfile)s.log'''

    P.run()

#########################################
#########################################
#########################################

@follows(mkdir("scatterplot_genus_cog_fold.dir"))
@split([runMetagenomeSeqPerCOGAndTaxa, 
        os.path.join(
            PARAMS.get("rna_communitiesdir"), 
            "genes.dir/gene_counts.diff.tsv"),
        os.path.join(
            PARAMS.get("dna_communitiesdir"),
            "genes.dir/gene_counts.diff.tsv")],
       "scatterplot_genus_cog_fold.dir/*scatters.pdf")
def scatterplotPerCogTaxaDNAFoldRNAFold(infiles, outfiles):
    '''
    plot fold change of RNA and DNA for NOGs
    of interest
    '''
    taxa_cog_rnadiff, taxa_cog_dnadiff, cog_rnadiff, cog_dnadiff = infiles
    PipelineMetaomics.scatterplotPerCogTaxaDNAFoldRNAFold(taxa_cog_rnadiff,
                                                          taxa_cog_dnadiff,
                                                          cog_rnadiff,
                                                          cog_dnadiff)

#########################################
#########################################
#########################################

@follows(plotMaxTaxaContribution,
         testSignificanceOfMaxTaxaContribution,
         heatmapTaxaCogProportionMatrix,
         runMetagenomeSeqPerCOGAndTaxa,
         scatterplotPerCogTaxaDNAFoldRNAFold)
def genus_contributions_to_cogs():
    pass

#########################################
#########################################
#########################################

@follows(diff_stats,
         overlaps,
         PCA,
         define_colitis_responsive_nogs,
         associate_taxa,
         genus_contributions_to_cogs)
def full():
    pass

#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
