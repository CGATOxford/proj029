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

@follows(mkdir("annotations.dir"))
@merge([PARAMS.get("annotations_rhea_file"), PARAMS.get("annotations_microme_file"), PARAMS.get("annotations_genomes")],
       "annotations.dir/pyruvate_oxidase.tsv")
def buildPyruvateOxidaseAnnotations(infiles, outfile):
    '''
    build a set of annotations as to whether or not genomes
    from genera that we analyse contain pyruvate oxidase 
    reactions
    '''
    rhea2kegg, microme, genomes = infiles
    statement = '''cat %(rhea2kegg)s | python %(projscripts)s/rhea2annotations.py 
                  --microme-file=%(microme)s
                  --genomes=%(genomes)s
                  --kegg-id=R00207
                  --log=%(outfile)s.log
                  > %(outfile)s
               '''
    P.run()

###################################################
###################################################
###################################################

@follows(mkdir("annotations.dir"))
@merge([PARAMS.get("annotations_rhea_file"), PARAMS.get("annotations_microme_file"), PARAMS.get("annotations_genomes")],
       "annotations.dir/catalase.tsv")
def buildCatalaseAnnotations(infiles, outfile):
    '''
    build a set of annotations as to whether or not genomes
    from genera that we analyse contain pyruvate oxidase 
    reactions
    '''
    rhea2kegg, microme, genomes = infiles
    statement = '''cat %(rhea2kegg)s | python %(projscripts)s/rhea2annotations.py 
                  --microme-file=%(microme)s
                  --genomes=%(genomes)s
                  --kegg-id=R00009
                  --log=%(outfile)s.log
                  > %(outfile)s
               '''
    P.run()

###################################################
###################################################
###################################################
# investigating the differences in taxonomic
# abundances. 
###################################################
###################################################
###################################################

@follows(mkdir("diff_stats.dir"))
@transform(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.diff.load"))
           + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.diff.load")),
           regex("(\S+)/(.*NA)/(\S+)/genus.diamond.aggregated.counts.diff.load"),
           r"diff_stats.dir/\2.diff.stats.tsv")
def buildGeneraDiffStats(infile, outfile):
    '''
    build stats on the number of genera that are differentially
    abundant at different levels
    '''
    dbh = connect()
    cc = dbh.cursor()
    db = os.path.join(P.snip(os.path.dirname(infile), "counts.dir"), "csvdb")
    tablename = P.toTable(os.path.basename(infile))
    statement = "ATTACH '%(db)s' as diff;" % locals() 
    cc.execute(statement)
    
    # build table of results at different thresholds
    ps = [0.01, 0.05, 0.1]
    fcs = [0, 0.5, 1, 1.5, 2]

    # build results for each pair
    pairs = [("HhaIL10R", "WT"), ("WT", "aIL10R"), ("Hh", "WT")]

    outf = open(outfile, "w")
    outf.write("group1\tgroup2\tadj_P_Val\tlogFC\tnumber\n")

    for pair in pairs:
        p1, p2 = pair[0], pair[1]
        for p, fc in itertools.product(ps, fcs):
            statement = """SELECT COUNT(*) 
                           FROM diff.%(tablename)s
                           WHERE group1 == "%(p1)s"
                           AND group2 == "%(p2)s"
                           AND adj_P_Val < %(p)f
                           AND abs(logFC) > %(fc)f""" % locals()

            for data in cc.execute(statement).fetchall():
                outf.write("\t".join([p1, p2, str(p), str(fc), str(data[0])]) + "\n") 
    outf.close()
            
###################################################
###################################################
###################################################

@follows(mkdir("diff_stats.dir"))
@transform(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.diff.load"))
           + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.diff.load")),
           regex("(\S+)/(.*NA)/(\S+)/gene_counts.diff.load"),
           r"diff_stats.dir/\2.genes.diff.stats.tsv")
def buildGeneDiffStats(infile, outfile):
    '''
    build stats on the number of genera that are differentially
    abundant at different levels
    '''
    dbh = connect()
    cc = dbh.cursor()
    db = os.path.join(P.snip(os.path.dirname(infile), "genes.dir"), "csvdb")
    tablename = P.toTable(os.path.basename(infile))
    statement = "ATTACH '%(db)s' as diff;" % locals() 
    cc.execute(statement)
    
    # build table of results at different thresholds
    ps = [0.01, 0.05, 0.1]
    fcs = [0, 0.5, 1, 1.5, 2]

    # build results for each pair
    pairs = [("HhaIL10R", "WT"), ("WT", "aIL10R"), ("Hh", "WT")]

    outf = open(outfile, "w")
    outf.write("group1\tgroup2\tadj_P_Val\tlogFC\tnumber\n")

    for pair in pairs:
        p1, p2 = pair[0], pair[1]
        for p, fc in itertools.product(ps, fcs):
            statement = """SELECT COUNT(*) 
                           FROM diff.%(tablename)s
                           WHERE group1 == "%(p1)s"
                           AND group2 == "%(p2)s"
                           AND adj_P_Val < %(p)f
                           AND abs(logFC) > %(fc)f""" % locals()

            for data in cc.execute(statement).fetchall():
                outf.write("\t".join([p1, p2, str(p), str(fc), str(data[0])]) + "\n") 
    outf.close()
 

@follows(buildGeneraDiffStats, 
         buildGeneDiffStats)
def diff_stats():
    pass


#########################################
#########################################
#########################################

@follows(mkdir("diff.dir"))
@transform([os.path.join(PARAMS.get("rna_communitiesdir"), "csvdb"),
            os.path.join(PARAMS.get("dna_communitiesdir"), "csvdb")],
           regex("(\S+)/(.*NA).*/csvdb"),
           r"diff.dir/\2_HhaIL10R_vs_WT.diff.genus.tsv")
def buildGenusDiffList(infile, outfile):
    '''
    build a list of differentially expressed genera
    between HhaIL10R and WT
    '''
    dbh = sqlite3.connect(infile)
    cc = dbh.cursor()
    outf = open(outfile, "w")

    for data in cc.execute("""SELECT taxa, logFC 
                              FROM genus_diamond_aggregated_counts_diff
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              AND P_Value < 0.05""").fetchall():
        if data[1] > 0:
            status = 1
        elif data[1] < 0:
            status = 2

        outf.write("%s\t%s\n" % (data[0], status))
    outf.close()

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@follows(mkdir("enrichment.dir"))
@transform(buildGenusDiffList, regex("(\S+)/(\S+).tsv"),
           add_inputs(buildPyruvateOxidaseAnnotations),
           r"enrichment.dir/\2_diff_pyruvate_oxidase_up.tsv")
def geDiffAndPyruvateOxidaseUp(infiles, outfile):
    '''
    output those genera that are up-regulated and
    have pyruvate oxidase encoded function
    '''
    R('''diff <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''anno <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])

    # get annotated genomes with pyruvate oxidase
    R('''withpo <- anno[,1][anno[,2] == "yes"]''')

    # get upregulated genera
    R('''up <- diff[,1][diff[,2] == "1"]''')

    # get intersection
    R('''upwithpo <- data.frame(intersect(withpo, up))''')

    # write results
    R('''write.table(upwithpo, file = "%s", quote = F, row.names = F, sep = "\t")''' % outfile)

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@follows(mkdir("enrichment.dir"))
@transform(buildGenusDiffList, regex("(\S+)/(\S+).tsv"),
           add_inputs(buildPyruvateOxidaseAnnotations),
           r"enrichment.dir/\2_pyruvate_oxidase_up.sig")
def testEnrichmentOfPyruvateOxidaseUp(infiles, outfile):
    '''
    test whether diff genus are associated with 
    pyruvate oxidase activity - using fishers exact test
    '''
    R('''diff <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''anno <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])

    # get annotated genomes with pyruvate oxidase and without
    R('''withpo <- anno[,1][anno[,2] == "yes"]''')
    R('''nopo <- anno[,1][anno[,2] == "no"]''')

    # get upregulated genera
    R('''up <- diff[,1][diff[,2] == "1"]''')
    R('''down <- diff[,1][diff[,2] == "2"]''')
    
    # build matrix for testing 
    # the use of "po" stands for pyruvate oxidase
    R('''upwithpo <- length(intersect(up, withpo))''')
    R('''upnopo <- length(intersect(up, nopo))''')
    R('''notupwithpo <- length(setdiff(withpo, upwithpo))''')
    R('''notupnopo <- length(setdiff(nopo, upnopo))''')

    R('''dat <- data.frame("up" = c(upwithpo, upnopo),
                           "notup" = c(notupwithpo, notupnopo))''')

    # barplot those with pyruvate oxidase
    R('''toplot <- data.frame("percentage" = c(dat[1,2] / (dat[2,2] + dat[1,2])*100,
                              dat[1,1] / (dat[1,1] + dat[2,1])*100))''')
    R('''toplot$status <- c("background", "increased")''')
    R('''print(toplot)''')

    outname_barplot = P.snip(outfile, ".sig") + ".pdf"
    R('''library(ggplot2)''')
    R('''ggplot(toplot, aes(x = status, y = percentage)) + geom_bar(position = "dodge", stat = "identity")''')
    R('''ggsave("%s")''' % outname_barplot)

    # perform test and write results
    R('''ftest <- fisher.test(dat)''')
    R('''res <- data.frame("estimate" = ftest$estimate, "p-value" = ftest$p.value)''')
    R('''write.table(res, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@follows(mkdir("enrichment.dir"))
@transform(buildGenusDiffList, regex("(\S+)/(\S+).tsv"),
           add_inputs(buildCatalaseAnnotations),
           r"enrichment.dir/\2_catalase_up.sig")
def testEnrichmentOfCatalaseUp(infiles, outfile):
    '''
    test whether diff genus are associated with 
    superoxide dismutase activity - using fishers exact test
    '''
    R('''diff <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''anno <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])

    # get annotated genomes with superoxide dismutase and without
    R('''withpo <- anno[,1][anno[,2] == "yes"]''')
    R('''nopo <- anno[,1][anno[,2] == "no"]''')

    # get upregulated genera
    R('''up <- diff[,1][diff[,2] == "1"]''')
    R('''down <- diff[,1][diff[,2] == "2"]''')
    
    # build matrix for testing 
    # the use of "po" in this case refers to superoxide dismutase
    R('''upwithpo <- length(intersect(up, withpo))''')
    R('''upnopo <- length(intersect(up, nopo))''')
    R('''notupwithpo <- length(setdiff(withpo, upwithpo))''')
    R('''notupnopo <- length(setdiff(nopo, upnopo))''')

    R('''dat <- data.frame("up" = c(upwithpo, upnopo),
                           "notup" = c(notupwithpo, notupnopo))''')

    # barplot those with pyruvate oxidase
    R('''toplot <- data.frame("percentage" = c(dat[1,2] / (dat[2,2] + dat[1,2])*100,
                              dat[1,1] / (dat[1,1] + dat[2,1])*100))''')
    R('''toplot$status <- c("background", "increased")''')
    R('''print(toplot)''')

    outname_barplot = P.snip(outfile, ".sig") + ".pdf"
    R('''library(ggplot2)''')
    R('''ggplot(toplot, aes(x = status, y = percentage)) + geom_bar(position = "dodge", stat = "identity")''')
    R('''ggsave("%s")''' % outname_barplot)

    # perform test and write results
    R('''ftest <- fisher.test(dat)''')
    R('''res <- data.frame("estimate" = ftest$estimate, "p-value" = ftest$p.value)''')
    R('''write.table(res, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@follows(mkdir("enrichment.dir"))
@transform(buildGenusDiffList, regex("(\S+)/(\S+).tsv"),
           add_inputs(buildPyruvateOxidaseAnnotations),
           r"enrichment.dir/\2_pyruvate_oxidase_down.sig")
def testEnrichmentOfPyruvateDown(infiles, outfile):
    '''
    test whether diff genus are associated with 
    pyruvate oxidase activity - using fishers exact test
    '''
    R('''diff <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''anno <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])

    # get annotated genomes with pyruvate oxidase and without
    R('''withpo <- anno[,1][anno[,2] == "yes"]''')
    R('''nopo <- anno[,1][anno[,2] == "no"]''')

    # get downregulated genera
    R('''up <- diff[,1][diff[,2] == "1"]''')
    R('''down <- diff[,1][diff[,2] == "2"]''')
    
    R('''print(down)''')

    # build matrix for testing 
    # the use of "po" stands for pyruvate oxidase
    R('''downwithpo <- length(intersect(down, withpo))''')
    R('''downnopo <- length(intersect(down, nopo))''')
    R('''notdownwithpo <- length(setdiff(withpo, downwithpo))''')
    R('''notdownnopo <- length(setdiff(nopo, downnopo))''')

    R('''dat <- data.frame("down" = c(downwithpo, downnopo),
                           "notdown" = c(notdownwithpo, notdownnopo))''')
    R('''print(dat)''')

    # barplot those with pyruvate oxidase
    R('''toplot <- data.frame("percentage" = c(dat[1,2] / (dat[2,2] + dat[1,2])*100,
                              dat[1,1] / (dat[1,1] + dat[2,1])*100))''')
    R('''toplot$status <- c("background", "decreased")''')

    outname_barplot = P.snip(outfile, ".sig") + ".pdf"
    R('''library(ggplot2)''')
    R('''ggplot(toplot, aes(x = status, y = percentage)) + geom_bar(position = "dodge", stat = "identity")''')
    R('''ggsave("%s")''' % outname_barplot)

    # perform test and write results
    R('''ftest <- fisher.test(dat)''')
    R('''res <- data.frame("estimate" = ftest$estimate, "p-value" = ftest$p.value)''')
    R('''write.table(res, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@follows(mkdir("enrichment.dir"))
@merge([buildGenusDiffList, buildPyruvateOxidaseAnnotations],
           r"enrichment.dir/intersection_pyruvate_oxidase.sig")
def testEnrichmentOfPyruvateOxidaseIntersection(infiles, outfile):
    '''
    test the enrichment fo aerobic bacteria in the intersection
    of DNA and RNA differences - upregulated only
    '''
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])

    # get the intersection of upregulated genera
    R('''diff <- intersect(dna[,1][dna[,2] == "1"], rna[,1][rna[,2] == "1"])''')

    # read the annotation file
    R('''anno <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[2])

    # get annotated genomes with pyruvate oxidase and without
    R('''withpo <- anno[,1][anno[,2] == "yes"]''')
    R('''nopo <- anno[,1][anno[,2] == "no"]''')

    # build matrix for testing 
    # the use of "po" stands for pyruvate oxidase
    R('''diffwithpo <- length(intersect(diff, withpo))''')
    R('''diffnopo <- length(intersect(diff, nopo))''')
    R('''notdiffwithpo <- length(setdiff(withpo, diffwithpo))''')
    R('''notdiffnopo <- length(setdiff(nopo, diffnopo))''')

    R('''dat <- data.frame("diff" = c(diffwithpo, diffnopo),
                           "notdiff" = c(notdiffwithpo, notdiffnopo))''')
    R('''print(dat)''')

    # barplot those with pyruvate oxidase
    R('''toplot <- data.frame("percentage" = c(dat[1,2] / (dat[2,2] + dat[1,2])*100,
                              dat[1,1] / (dat[1,1] + dat[2,1])*100))''')
    R('''toplot$status <- c("background", "intersection_up")''')

    outname_barplot = P.snip(outfile, ".sig") + ".pdf"
    R('''library(ggplot2)''')
    R('''ggplot(toplot, aes(x = status, y = percentage)) + geom_bar(position = "dodge", stat = "identity")''')
    R('''ggsave("%s")''' % outname_barplot)

    # perform test and write results
    R('''ftest <- fisher.test(dat)''')
    R('''res <- data.frame("estimate" = ftest$estimate, "p-value" = ftest$p.value)''')
    R('''write.table(res, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("diff.dir"))
@merge([os.path.join(PARAMS.get("rna_communitiesdir"), "csvdb"),
            os.path.join(PARAMS.get("dna_communitiesdir"), "csvdb")],
           "diff.dir/common_genes.tsv")
def buildCommonGeneList(infiles, outfile):
    '''
    get a list of genes that were common beteween DNA
    and RNA analysis - used for downstream analysis
    '''
    dbh_rna = sqlite3.connect(infiles[0])
    cc_rna = dbh_rna.cursor()
    dbh_dna = sqlite3.connect(infiles[1])
    cc_dna = dbh_dna.cursor()

    outf = open(outfile, "w")
    rna = set()
    dna = set()
    for gene in cc_rna.execute("""SELECT taxa 
                              FROM gene_counts_diff 
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              """).fetchall():
        rna.add(gene[0])

    for gene in cc_dna.execute("""SELECT taxa 
                              FROM gene_counts_diff 
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              """).fetchall():
        dna.add(gene[0])

    for gene in rna.intersection(dna):
        outf.write(gene + "\n")

#########################################
#########################################
#########################################

@follows(mkdir("diff.dir"))
@merge([os.path.join(PARAMS.get("rna_communitiesdir"), "csvdb"),
            os.path.join(PARAMS.get("dna_communitiesdir"), "csvdb")],
           "diff.dir/common_genera.tsv")
def buildCommonGeneraList(infiles, outfile):
    '''
    get a list of genera that were common beteween DNA
    and RNA analysis - used for downstream analysis
    '''
    dbh_rna = sqlite3.connect(infiles[0])
    cc_rna = dbh_rna.cursor()
    dbh_dna = sqlite3.connect(infiles[1])
    cc_dna = dbh_dna.cursor()

    outf = open(outfile, "w")
    rna = set()
    dna = set()
    for gene in cc_rna.execute("""SELECT taxa 
                              FROM genus_diamond_aggregated_counts_diff 
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              """).fetchall():
        rna.add(gene[0])

    for gene in cc_dna.execute("""SELECT taxa 
                              FROM genus_diamond_aggregated_counts_diff
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              """).fetchall():
        dna.add(gene[0])

    for gene in rna.intersection(dna):
        outf.write(gene + "\n")


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
    # list for sql statement
    common = set([x[:-1] for x in open(infiles[1]).readlines()])
    common = "(" + ",".join(['"'+x+'"' for x in common]) + ")"

    # connect to database
    dbh = sqlite3.connect(infiles[0])
    cc = dbh.cursor()

    # remove any genes that are different between Hh and steady state
    # or between aIL10R and steady state
    # also only include genes that are different with between colitis
    # and other control conditions
    # hh = set([x[0] for x in cc.execute("""SELECT taxa 
    #                    FROM gene_counts_diff
    #                    WHERE group1 == "Hh" 
    #                    AND group2 == "WT" 
    #                    AND adj_P_Val < 0.05""").fetchall()])

    # hh = "(" + ",".join(['"'+x+'"' for x in hh]) + ")"

    # ail10r = set([x[0] for x in cc.execute("""SELECT taxa 
    #                    FROM gene_counts_diff
    #                    WHERE group1 == "WT" 
    #                    AND group2 == "aIL10R" 
    #                    AND adj_P_Val < 0.05""").fetchall()])

    # ail10r = "(" + ",".join(['"'+x+'"' for x in ail10r]) + ")"

    outf = open(outfile, "w")
    for gene in cc.execute("""SELECT taxa 
                              FROM gene_counts_diff 
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              AND adj_P_Val < 0.05
                              AND (logFC > 1 OR logFC < -1)
                              AND taxa IN %s
                              ORDER BY logFC DESC""" % common).fetchall():
        outf.write(gene[0] + "\n")
    outf.close()

#                              AND taxa NOT IN %s
#                              AND taxa NOT IN %s

#########################################
#########################################
#########################################

@merge(buildGeneDiffList, "diff.dir/genes_overlap.tsv")
def buildDiffOverlap(infiles, outfile):
    '''
    overlap gene lists
    '''
    dna = [set(open(x).readlines()) for x in infiles if "DNA" in x][0]
    rna = [set(open(x).readlines()) for x in infiles if "RNA" in x][0]
    ndna = len(dna)
    nrna = len(rna)
    overlap = len(dna.intersection(rna))
    outf = open(outfile, "w")
    outf.write("nDNA\tnRNA\tnoverlap\n%(ndna)i\t%(nrna)i\t%(overlap)i\n" % locals())
    outf.close()

#########################################
#########################################
#########################################

@merge([buildCommonGeneList, buildDiffOverlap],
       "diff.dir/genes_overlap.sig")
def testSignificanceOfGenesOverlap(infiles, outfile):
    '''
    test significance of overlapping diff gene
    lists bewteen RNA and DNA using hypergeometric
    test
    '''
    common, overlap = infiles
    R('''pop <- read.csv("%s", header = F, sep = "\t", stringsAsFactors = F)''' % common)
    R('''overlaps <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % overlap)

    # total genes in population
    R('''npop <- nrow(pop)''')

    # x = number of white balls picked = overlap
    R('''x <- overlaps$noverlap''')

    # m = total number of white balls = total diff in RNA analysis
    R('''m <- overlaps$nRNA''')

    # n = total number of black balls = total - diff in RNA analysis
    R('''n <- npop - m''')

    # k = total balls sampled = number of genera different in DNA analysis
    R('''k <- overlaps$nDNA''')

    # hypergeometric test
    R('''p <- 1-phyper(x,m,n,k)''')

    # write result
    R('''res <- matrix(ncol = 2, nrow = 5)''')
    R('''res[1,1] <- "x"''')
    R('''res[2,1] <- "m"''')
    R('''res[3,1] <- "n"''')
    R('''res[4,1] <- "k"''')
    R('''res[5,1] <- "p-value"''')
    R('''res[1,2] <- x''')
    R('''res[2,2] <- m''')
    R('''res[3,2] <- n''')
    R('''res[4,2] <- k''')
    R('''res[5,2] <- p''')
    R('''write.table(as.data.frame(res), file = "%s", quote = F, sep = "\t", row.names = F)''' % outfile)


#########################################
#########################################
#########################################

@merge([os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.diff.tsv"),
        os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.diff.tsv")],
       "diff.dir/genus_overlap.sig")
def testSignificanceOfDiffGenusOverlap(infiles, outfile):
    '''
    use hypergeometric test to determine significance
    of overlap between differential genera in DNA and
    RNA analyses
    '''
    dna, rna = infiles
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''dna <- dna[dna$group1 == "HhaIL10R" & dna$group2 == "WT",]''')
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''rna <- rna[rna$group1 == "HhaIL10R" & rna$group2 == "WT",]''')

    # must use commonly detected
    R('''common <- intersect(rna$taxa, dna$taxa)''')
    R('''dna <- dna[dna$taxa %in% common,]''')
    R('''rna <- rna[rna$taxa %in% common,]''')

    # total population = common genera
    R('''npop <- length(common)''')

    # x = number of white balls picked = overlap
    R('''x <- length(intersect(dna$taxa[dna$P.Value < 0.05], rna$taxa[rna$P.Value < 0.05]))''')

    # m = total number of white balls = total diff in RNA analysis
    R('''m <- length(rna$taxa[rna$P.Value < 0.05])''')

    # n = total number of black balls = total - diff in RNA analysis
    R('''n <- npop - m''')

    # k = total balls sampled = number of genera different in DNA analysis
    R('''k <- length(dna$taxa[dna$P.Value < 0.05])''')

    # hypergeometric test
    R('''p <- 1-phyper(x,m,n,k)''')

    # write result
    R('''res <- matrix(ncol = 2, nrow = 5)''')
    R('''res[1,1] <- "x"''')
    R('''res[2,1] <- "m"''')
    R('''res[3,1] <- "n"''')
    R('''res[4,1] <- "k"''')
    R('''res[5,1] <- "p-value"''')
    R('''res[1,2] <- x''')
    R('''res[2,2] <- m''')
    R('''res[3,2] <- n''')
    R('''res[4,2] <- k''')
    R('''res[5,2] <- p''')
    R('''write.table(as.data.frame(res), file = "%s", quote = F, sep = "\t", row.names = F)''' % outfile)

#########################################
#########################################
#########################################

@follows(buildGenusDiffList, buildGeneDiffList, buildDiffOverlap)
def diff_lists():
    pass

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@transform(buildGeneDiffList,
           suffix(".tsv"),
           add_inputs([os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
                      os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.norm.matrix")]),
           ".png")
def heatMapDiffGenes(infiles, outfile):
    '''
    heatmap differences between WT and HhaIL10R groups
    '''
    # we do this  for the different combinations to see how the
    # DNA and RNA differences compare

    R('''diff <- read.csv("%s", header = F, sep = "\t", stringsAsFactors = F)''' % infiles[0])
    R('''library(gplots)''')
    R('''library(gtools)''')

    for mat in infiles[1]:
        if "RNA" in mat:
            outname = os.path.dirname(outfile) + "/RNA_" + os.path.basename(outfile)
        elif "DNA" in mat:
            outname = os.path.dirname(outfile) + "/DNA_" + os.path.basename(outfile)
        R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % mat)

        R('''rownames(dat) <- dat$taxa''')
        R('''dat <- dat[, 1:ncol(dat)-1]''')
        R('''dat <- dat[diff[,1],]''')
        R('''dat <- na.omit(dat)''')
        R('''colnames(dat) <- unlist(strsplit(colnames(dat), ".diamond_count"))''')
        R('''dat <- dat[, mixedsort(colnames(dat))]''')
        R('''samples <- colnames(dat)''')
        R('''dat <- t(apply(dat, 1, scale))''')
        R('''colnames(dat) <- samples''')
        R('''cols <- colorRampPalette(c("blue", "white", "red"))''')
        R('''png("%s")''' % outname)
        R('''heatmap.2(as.matrix(dat), col = cols, scale = "row", trace = "none", Rowv = F, Colv = F, margins = c(15,15), 
                       distfun = function(x) dist(x, method = "manhattan"),
                       hclustfun = function(x) hclust(x, method = "ward.D2"))''')
        R["dev.off"]()

#########################################
#########################################
#########################################

@follows(mkdir("pca.dir"))
@jobs_limit(1, "R")
@transform([os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
            os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
            os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix"),
            os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix")],
           regex("(\S+)/(\S+).matrix"),
           r"pca.dir/\2.loadings.tsv")
def buildPCALoadings(infile, outfile):
    '''
    run PCA and heatmap the loadings
    '''
    if "RNA" in infile:
        suffix = "rna"
    else:
        suffix = "dna"
    if "gene" in infile:
        xlim, ylim = 40,40
    else:
        xlim, ylim = 12,7
    
    outname_plot = P.snip(outfile, ".loadings.tsv").replace("/", "/%s_" % suffix) + ".pca.pdf"
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[, 1:ncol(dat)-1]''')
    R('''pc <- prcomp(t(dat))''')
    R('''conds <- unlist(strsplit(colnames(dat), ".R[0-9]"))[seq(1, ncol(dat)*2, 2)]''')
    R('''conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]''')  
        
    # plot the principle components
    R('''library(ggplot2)''')
    R('''pcs <- data.frame(pc$x)''')
    R('''pcs$cond <- conds''')
    # get variance explained
    R('''imps <- c(summary(pc)$importance[2], summary(pc)$importance[5])''')
    R('''p <- ggplot(pcs, aes(x = PC1, y = PC2, colour = cond, size = 3)) + geom_point()''')
    R('''p2 <- p + xlab(imps[1]) + ylab(imps[2])''')
    R('''p3 <- p2 + scale_colour_manual(values = c("slateGrey", "green", "red", "blue"))''')
    R('''p3 + xlim(c(-%i, %i)) + ylim(c(-%i, %i))''' % (xlim, xlim, ylim, ylim))
    R('''ggsave("%s")''' % outname_plot)

    # get the loadings
    R('''loads <- data.frame(pc$rotation)''')
    R('''loads$taxa <- rownames(loads)''')

    # write out data
    R('''write.table(loads, file = "%s", sep = "\t", row.names = F, quote = F)''' % outfile.replace("/", "/%s_" % suffix))

    P.touch(outfile)

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@follows(buildPCALoadings)
@transform("pca.dir/*na_*loadings.tsv", suffix(".tsv"), ".pdf")
def plotPCALoadings(infile, outfile):
    '''
    plot the loadings associated with differentially
    expressed genes
    '''
    R('''library(ggplot2)''')
    R('''library(grid)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    
    R('''top5pc1 <- dat[order(-dat$PC1),][1:5,]''')
    R('''bottom5pc1 <- dat[order(dat$PC1),][1:5,]''')
    R('''top5pc2 <- dat[order(-dat$PC2),][1:5,]''')
    R('''bottom5pc2 <- dat[order(dat$PC2),][1:5,]''')
    R('''totext <- data.frame(rbind(top5pc1, bottom5pc1, top5pc2, bottom5pc2))''')
 
    R('''dat$x <- 0''')
    R('''dat$y <- 0''')
    R('''p <- ggplot(dat, aes(x = x, y = y, xend = PC1, yend = PC2, colour = taxa))''')
    R('''p2 <- p + geom_segment(arrow = arrow(length = unit(0.2, "cm")))''')
    R('''p2 + geom_text(data = totext, aes(x = PC1, y = PC2, label = totext$taxa, size = 6)) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.25))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("indicators.dir"))
@jobs_limit(1, "R")
@transform([os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix"),
            os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix")],
           regex("(\S+)/(\S+).matrix"),
           r"indicators.dir/\2.indicators.pdf")
def plotIndicatorGenera(infile, outfile):
    '''
    plot genera that load onto PCA
    '''
    if "RNA" in infile:
        suffix = "rna"
    else:
        suffix = "dna"

    R('''library(reshape)''')
    R('''library(ggplot2)''')
    R('''library(plyr)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[c("Akkermansia", "Escherichia", "Peptoniphilus"),]''')
    R('''dat2 <- melt(dat)''')
    R('''conds <- unlist(strsplit(as.character(dat2$variable), ".R[0-9]"))''')
    R('''conds <- conds[seq(1,length(conds),2)]''')
    R('''dat2$cond <- conds''')
    R('''dat2.sum <- ddply(dat2, c("cond", "taxa"), summarize, mean = mean(value), n = length(cond), sd = sd(value), se = sd/sqrt(n))''')
#    R('''dat2.sum <- dat2.sum[dat2.sum$cond == "stool.WT" | dat2.sum$cond == "stool.HhaIL10R",]''')
    R('''dodge = position_dodge(width=0.9)''')
    R('''plot1 <- ggplot(dat2.sum, aes(x = taxa, y = mean, fill = cond))''')
    R('''plot2 <- plot1 + geom_bar(stat = "identity", position = dodge)''') 
    R('''plot3 <- plot2 + geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.25, position = dodge)''') 
    R('''plot3  + scale_fill_manual(values = c("grey", "darkGreen","red", "darkBlue")) + theme_bw()''')
    
    outname = P.snip(outfile, ".pdf") + "." + suffix + ".pdf"
    R('''ggsave("%s")''' % outname)

#########################################
#########################################
#########################################

@jobs_limit(1, "R")
@transform(buildGenusDiffList,
           suffix(".tsv"),
           add_inputs([os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix"),
                      os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix")]),
           ".pdf")
def heatMapDiffGenera(infiles, outfile):
    '''
    heatmap differences between WT and HhaIL10R groups
    '''
    # we do this  for the different combinations to see how the
    # DNA and RNA differences compare

    R('''diff <- read.csv("%s", header = F, sep = "\t", stringsAsFactors = F)''' % infiles[0])
    R('''library(gplots)''')
    R('''library(gtools)''')

    for mat in infiles[1]:
        if "RNA" in mat:
            outname = os.path.dirname(outfile) + "/RNA_" + os.path.basename(outfile)
        elif "DNA" in mat:
            outname = os.path.dirname(outfile) + "/DNA_" + os.path.basename(outfile)
        R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % mat)
        R('''rownames(dat) <- dat$taxa''')
        R('''dat <- dat[, 1:ncol(dat)-1]''')
        R('''dat <- dat[diff[,1],]''')
        R('''dat <- na.omit(dat)''')
        R('''colnames(dat) <- unlist(strsplit(colnames(dat), ".diamond_count"))''')
        R('''dat <- dat[, mixedsort(colnames(dat))]''')
        R('''samples <- colnames(dat)''')
        R('''dat <- t(apply(dat, 1, scale))''')
        R('''colnames(dat) <- samples''')
        R('''cols <- colorRampPalette(c("blue", "white", "red"))''')
        R('''pdf("%s")''' % outname)
        R('''heatmap.2(as.matrix(dat), col = cols, scale = "row", trace = "none", Rowv = T, Colv = F, margins = c(15,15), 
                       distfun = function(x) dist(x, method = "manhattan"),
                       hclustfun = function(x) hclust(x, method = "ward.D2"))''')
        R["dev.off"]()

@follows(heatMapDiffGenera,heatMapDiffGenes)
def heatmap():
    pass

###################################################
###################################################
###################################################
# build ratio of RNA to DNA and see if there is
# a difference between those that are called 
# as differentially abundant in RNA analysis
# vs DNA analysis
###################################################
###################################################
###################################################

@follows(mkdir("rna_dna_ratio.dir"))
@merge(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.diff.tsv")) 
       + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.diff.tsv")),
       "rna_dna_ratio.dir/ratio_genes.tsv")
def buildRNADNARatio(infiles, outfile):
    '''
    build the ratio between RNA expression
    and DNA abundance - from normalised
    counts
    '''
    rna = [x for x in infiles if "RNA" in x][0]
    dna = [x for x in infiles if "DNA" in x][0]

    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''rna <- rna[rna$group1 == "HhaIL10R" & rna$group2 == "WT",]''')
    R('''dna <- dna[dna$group1 == "HhaIL10R" & dna$group2 == "WT",]''')
    
    R('''rownames(rna) <- rna$taxa''')
    R('''rownames(dna) <- dna$taxa''')

    R('''rna <- rna[,1:ncol(rna)-1]''')
    R('''dna <- dna[,1:ncol(dna)-1]''')

    # only look at those that are present in both
    R('''keep <- intersect(rownames(rna), rownames(dna))''')
    R('''rna <- rna[keep,]''')
    R('''dna <- dna[keep,]''')

    R('''rna.ratio <- rna$logFC''')
    R('''dna.ratio <- dna$logFC''')
    R('''rna.p <- rna$adj.P.Val''')
    R('''dna.p <- dna$adj.P.Val''')
    
    R('''ratio <- data.frame(gene = keep, dna = dna.ratio, rna = rna.ratio, pdna = dna.p, prna = rna.p, ratio = rna.ratio - dna.ratio)''')
    R('''write.table(ratio, file = "%s", sep = "\t", row.names = F, quote = F)''' % outfile)


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
    rna_diff = set([y[:-1] for y in open([x for x in infiles if "RNA" in x][0]).readlines()])
    dna_diff = set([y[:-1] for y in open([x for x in infiles if "DNA" in x][0]).readlines()])

    inf = IOTools.openFile(infiles[0])
    inf.readline()
    outf = IOTools.openFile(outfile, "w")
    outf.write("gene\tdna\trna\tpdna\tprna\tratio\tstatus\n")
    for line in inf.readlines():
        gene, dna, rna, pdna, prna, ratio = line[:-1].split("\t")
        gene = gene.strip('"')
        dna, rna = float(dna), float(rna)
        if gene in rna_diff and gene in dna_diff and dna > 0 and rna > 0:
            status = "up.both"
        elif gene in rna_diff and gene in dna_diff and dna < 0 and rna < 0:
            status = "down.both"
        elif gene in rna_diff and rna > 0:
            status = "up.RNA"
        elif gene in rna_diff and rna < 0:
            status = "down.RNA"
        elif gene in dna_diff and dna > 0:
            status = "up.DNA"
        elif gene in dna_diff and dna < 0:
            status = "down.DNA"
        else:
            status = "NS"
        outf.write("%(gene)s\t%(dna)s\t%(rna)s\t%(pdna)s\t%(prna)s\t%(ratio)s\t%(status)s\n" % locals())
    outf.close()

###################################################
###################################################
###################################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".pdf")
def plotSets(infile, outfile):
    '''
    plot the fold changes in RNA and DNA analyses
    and label by how they are regulated in DNA and
    RNA analyses
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)

    R('''cog2gene <- read.csv("goi.tsv", header = F, stringsAsFactors = F, sep = "\t", row.names = 1)''')

    # just get those signficant in either DNA or RNA or both
    R('''dat$status[dat$status == "NS"] = "z"''')
    R('''genes <- dat$gene''')

    # regression model
    R('''mod1 <- lm(dat$rna~dat$dna)''')
    R('''intercept <- mod1[[1]][1]''')
    R('''slope = mod1[[1]][2]''')

    # prediction intervals
    R('''pred.ints <- predict(mod1, interval = "prediction", level = 0.95)''')

    # add to data.frame
    R('''dat$lwr <- pred.ints[,2]''')
    R('''dat$upr <- pred.ints[,3]''')
    
    # add labels
    R('''dat$goi <- cog2gene[dat$gene,]''')
    R('''dat$pointsize <- ifelse(!(is.na(dat$goi)), 10, 1)''')

    # plot
    R('''plot1 <- ggplot(dat, aes(x = dna, y = rna, alpha = 1, colour = status)) + geom_point(shape = 18, aes(size = pointsize)) + scale_size_area() + xlim(c(-5,5))''') 
    R('''plot2 <- plot1 + scale_colour_manual(values = c("blue", "brown","darkGreen", "orange", "purple", "red", "grey"))''')
    R('''plot3 <- plot2 + geom_abline(yintercept = intercept, slope = slope)''')

    # prediction intervals
    R('''plot4 <- plot3 + geom_line(aes(x = dna, y = lwr), linetype = "dashed", colour = "black")''')
    R('''plot5 <- plot4 + geom_line(aes(x = dna, y = upr), linetype = "dashed", colour = "black")''')
    R('''plot5 + geom_text(aes(label = goi))''')

    # R('''plot6 <- plot5 + geom_vline(xintercept = 0) + geom_vline(xintercept = c(-1,1), linetype = "dashed")''')
    # R('''plot6 + geom_hline(yintercept = 0) + geom_hline(yintercept = c(-1,1), linetype = "dashed")''')
    
    R('''ggsave("%s")''' % outfile)

###################################################
###################################################
###################################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".outsidepi.tsv")
def buildGenesOutsidePredictionInterval(infile, outfile):
    '''
    annotate genes as being outside prediction
    interval
    '''
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)

    # just get those signficant in either DNA or RNA or both
    R('''genes <- dat$gene''')

    # regression model
    R('''mod1 <- lm(dat$rna~dat$dna)''')

    # prediction intervals
    R('''pred.ints <- predict(mod1, interval = "prediction", level = 0.95)''')

    # add to data.frame
    R('''dat$lwr <- pred.ints[,2]''')
    R('''dat$upr <- pred.ints[,3]''')
    
    # annotate with whether or not they are above 
    # prediction intervals
    R('''dat$pi_status[dat$rna > dat$upr & dat$status == "up.RNA"] <- "diff.up.rna"''')
    R('''dat$pi_status[dat$rna > dat$upr & dat$status == "down.DNA"] <- "diff.down.dna"''')
    R('''dat$pi_status[dat$rna > dat$upr & dat$status == "up.both"] <- "diff.up.rna"''')

    R('''dat$pi_status[dat$rna < dat$lwr & dat$status == "down.RNA"] <- "diff.down.rna"''')
    R('''dat$pi_status[dat$rna < dat$lwr & dat$status == "up.DNA"] <- "diff.up.dna"''')
    R('''dat$pi_status[dat$rna < dat$lwr & dat$status == "down.both"] <- "diff.down.rna"''')

    # write results
    R('''write.table(dat, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

###################################################
###################################################
###################################################

@follows(mkdir("annotations.dir"))
@merge([buildGenesOutsidePredictionInterval, 
        PARAMS.get("genes_annotations")],
       "annotations.dir/diff_genes_annotated.tsv")
def annotateFunctions(infiles, outfile):
    '''
    annotate cogs with functional descriptions
    from annotations file
    '''
    inf, annotations = infiles
    inf = IOTools.openFile(inf)
    annotations = IOTools.openFile(annotations)

    header = inf.readline()
    cog2data = {}
    E.info("reading classifications")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        cog = data[0]
        cog2data[cog] = data[1:]

    E.info("assigning descriptions to COGs")
    outf = IOTools.openFile(outfile, "w")
    outf.write(header[:-1] + "\tdescription\n")
    cogs = set()
    for line in annotations.readlines():
        cog, description = line[:-1].split("\t")
        if cog in cog2data:
            cogs.add(cog)
            outf.write("\t".join(cog2data[cog] + [cog] + [description]) + "\n")

    # collect those without descriptions
    for cog, data in cog2data.iteritems():
        if cog not in cogs:
            outf.write("\t".join(cog2data[cog] + [cog] + ["no description"]) + "\n")
    outf.close()

###################################################
###################################################
###################################################

@transform(buildGenesOutsidePredictionInterval, suffix(".tsv"), ".pdf")
def plotGenesOutsidePredictionInterval(infile, outfile):
    '''
    place those on plot that are outside
    prediction intervals
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''dat$c.status <- paste(dat$status, dat$pi_status, sep = "_")''')

    # regression model
    R('''mod1 <- lm(dat$rna~dat$dna)''')
    R('''intercept <- mod1[[1]][1]''')
    R('''slope = mod1[[1]][2]''')

    R('''plot1 <- ggplot(dat, aes(x = dna, y = rna, colour = pi_status)) + geom_point()''')
    R('''plot2 <- plot1 + geom_abline(yintercept = intercept, slope = slope)''')

    # prediction intervals
    R('''plot3 <- plot2 + geom_line(aes(x = dna, y = lwr), linetype = "dashed", colour = "black")''')
    R('''plot4 <- plot3 + geom_line(aes(x = dna, y = upr), linetype = "dashed", colour = "black")''')

    R('''plot5 <- plot4 + geom_vline(xintercept = 0) + geom_vline(xintercept = c(-1,1), linetype = "dashed")''')
    R('''plot5 + geom_hline(yintercept = 0) + geom_hline(yintercept = c(-1,1), linetype = "dashed")''')
    R('''ggsave("%s")''' % outfile)


###################################################
###################################################
###################################################
# testing for enrichment of high level functional
# categories
###################################################
###################################################
###################################################

@follows(mkdir("enrichment.dir"))
@transform(buildGenesOutsidePredictionInterval, 
           regex("(\S+)/(\S+).tsv"), 
           r"enrichment.dir/\2.foreground")
def buildForegroundSets(infile, outfile):
    '''
    build foreground sets of COGs for those
    that show larger fold change in RNA
    analysis
    '''
    inf = IOTools.openFile(infile)
    header = inf.readline()
    outf = IOTools.openFile(outfile, "w")
    outf.write("gene\trna_up\trna_down\n")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        cog, status = data[0], data[9]
        if status == "diff.up.rna":
            outf.write("%s\t1\t0\n" % cog)
        elif status == "diff.down.rna":
            outf.write("%s\t0\t1\n" % cog)
        else:
            outf.write("%s\t0\t0\n" % cog)
    outf.close()

###################################################
###################################################
###################################################

@split([buildForegroundSets,
        buildCommonGeneList,
        PARAMS.get("pathways_geneset")],
       "enrichment.dir/*.overall")
def runEnrichmentAnalysis(infiles, outfiles):
    '''
    run pathways analysis
    '''
    genes, background, gene2pathway = infiles

    # remove General function prediction only
    # and Function unknown
    temp = P.getTempFilename(".")
    statement = '''cat %(gene2pathway)s | grep -v "General function"
                                      | grep -v "Function unknown"
                                       > %(temp)s'''
    P.run()

    statement = '''python %(scriptsdir)s/runGO.py \                                                                                                                                                                                    
                   --background=%(background)s
                   --genes=%(genes)s \                                                                                                                                                                                                              
                   --filename-input=%(temp)s \                                                                                                                                                                                                     
                   --fdr \                                                                                                                                                                                                                   
                   -q BH \                                                                                                                                                                                                                   
                   --output-filename-pattern="enrichment.dir/%%(set)s.%%(go)s.%%(section)s" \                                                                                                                                                  
                   > enrichment.dir/pathways.log  \                                                                                                                                                                                            
                ; rm -rf %(temp)s
                ''' 
    P.run()

###################################################
###################################################
###################################################

@merge([buildForegroundSets, PARAMS.get("pathways_geneset")],
       "enrichment.dir/cogs_pathways.tsv")
def buildDiffCogsAndPathways(infiles, outfile):
    '''
    merge diff COGs and pathways
    '''
    R('''cogs <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''pathways <- read.csv("%s", header = F, stringsAsFactors = F, sep = "\t")''' % infiles[1])
    R('''dat <- merge(cogs, pathways, by.x = "gene", by.y = "V2", all.x = T, all.y = F)''')
    R('''write.table(dat, file = "%s", sep = "\t", row.names = F, quote = F)''' % outfile)

###################################################
###################################################
###################################################

@transform(buildDiffCogsAndPathways, suffix(".tsv"), ".pdf")
def plotPathways(infile, outfile):
    '''
    plot pathways associated with clusters
    '''
    R('''library(plyr)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''dat <- dat[,1:ncol(dat)-1]''')

    R('''dat$cluster <- ifelse(dat$rna_up == 1, "Up", NA)''')
    R('''dat$cluster <- ifelse(dat$rna_down == 1, "Down", dat$cluster)''')
    
    # remove NAs and unknown functions
    R('''dat <- na.omit(dat)''')

    # get counts per cluster
    R('''counts <- ddply(dat, c("cluster"), summarise, n = length(cluster))''')
    R('''rownames(counts) <- counts$cluster''')

    # add counts to main data
    R('''dat$count <- counts[dat$cluster,]$n''')

    # summarise per functional category
    R('''dat2 <- ddply(dat, c("V4", "cluster"), summarise, prop = (length(cluster)/mean(count))*100)''')

#    R('''dat2 <- dat2[grep("Function unknown", dat2$V4, invert = T),]''')
#    R('''dat2 <- dat2[grep("General function", dat2$V4, invert = T),]''')

    R('''plot1 <- ggplot(dat2, aes(x = V4, y = prop, fill = cluster)) + geom_bar(stat = "identity", position = "dodge") + facet_grid(~cluster) + coord_flip()''')
    R('''plot1 + scale_fill_manual(values = c("darkGrey", "black"))''')
    R('''ggsave("%s", width = 15)''' % outfile)

###################################################
###################################################
###################################################

@follows(mkdir("ratio_histograms.dir"))
@transform(annotateRNADNARatio, regex("(\S+)/(\S+).tsv"), r"ratio_histograms.dir/\2.pdf")
def plotHistogramPerSet(infile, outfile):
    '''
    plot a histogram of RNA/DNA fold change ratios
    per set
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)

    R('''dat$status[dat$status == "NS"] <- "z"''')
    R('''plot1 <- ggplot(dat, aes(x = status, y = ratio, fill = status))''')
    R('''plot2 <- plot1 + geom_boxplot(linewidth = 2) + geom_hline(yintercept = 0)''')
    R('''plot2 + scale_fill_manual(values = c("blue", "brown","darkGreen", "orange", "purple", "red", "grey"))''')
    R('''ggsave("%s", width = 12)''' % outfile)


    # R('''for (s in sets){
    #           outname <- paste("ratio_histograms.dir", s, sep = "/")
    #           outname <- paste(outname, "pdf", sep = ".")
    #           col <- as.character(colors.frame[s,]$color)
    #           dat2 <- dat[dat$status == s,]
    #           ggplot(dat2, aes(x = ratio)) + geom_histogram() + scale_fill_manual(values = c(col))
    #           ggsave(outname)}''')



###################################################
###################################################
###################################################

@jobs_limit(1, "R")
@transform(runEnrichmentAnalysis, suffix(".overall"), ".pdf")
def plotFunctions(infile, outfile):
    '''
    plot the enrichment results
    '''
    R('''library("ggplot2")''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''dat <- dat[order(dat$ratio, decreasing = F),]''')
    R('''dat$sig <- ifelse(dat$fdr < 0.05, "yes", "no")''')
    R('''plot1 <- ggplot(dat, aes(x = factor(description, levels = dat$description), y = ratio, fill = sig))''') 
    R('''plot2 <- plot1 + geom_bar(stat = "identity", position = "dodge") + coord_flip()''')
    R('''plot2 + scale_fill_manual(values = c("grey", "darkRed")) + ylim(c(0,3))''')
    R('''ggsave("%s")''' % outfile)

###################################################
###################################################
###################################################

@follows(mkdir("merged_diff.dir"))
@merge([os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
        os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
        buildGeneDiffList],
       "merged_diff.dir/merged_diff.norm.matrix")
def buildRNADNAMergedMatrix(infiles, outfile):
    '''
    merge the normalised matrices for DNA and RNA
    for looking at covariation across data sets
    '''
    R('''library(gtools)''')
    
    mat_dna, mat_rna, diff_dna, diff_rna = infiles

    # read in data
    # R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % mat)
    # R('''rownames(dat) <- dat$taxa''')
    # R('''dat <- dat[,1:ncol(dat)-1]''')
    

    # NB this was code for merging sets
    # read in dna matrix
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % mat_dna)
    R('''rownames(dna) <- dna$taxa''')
    R('''dna <- dna[,1:ncol(dna)-1]''')

    # read in rna matrix
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % mat_rna)
    R('''rownames(rna) <- rna$taxa''')
    R('''rna <- rna[,1:ncol(rna)-1]''')

    # get differentially abundant genes
    R('''dna.diff <- read.csv("%s", header = F, sep = "\t")''' % diff_dna)
    R('''rna.diff <- read.csv("%s", header = F, sep = "\t")''' % diff_rna)
    R('''diff <- union(dna.diff[,1], rna.diff[,1])''')

    # subset by diff genes
    #R('''dat <- dat[diff,]''')
    
    R('''dna <- dna[diff,]''')
    R('''rna <- rna[diff,]''')

    R('''dna <- dna[,mixedsort(colnames(dna))]''')
    R('''rna <- rna[,mixedsort(colnames(rna))]''')

    # combine sets
    R('''dat <- data.frame(cbind(dna, rna))''')

    # sort out column orders
    #R('''dat <- dat[, mixedsort(colnames(dat))]''')

    R('''dat$gene <- rownames(dat)''')

    # write out results
    R('''write.table(dat, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

###################################################
###################################################
###################################################

@jobs_limit(1, "R")
@transform(buildRNADNAMergedMatrix, 
           suffix(".matrix"), 
           add_inputs(annotateRNADNARatio),
           ".png")
def heatmapMergedMatrix(infiles, outfile):
    '''
    heatmap DNA and RNA together
    '''
    R('''library(gplots)''')

    # read in data
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''rownames(dat) <- dat$gene''')
    R('''dat <- dat[,1:ncol(dat)-1]''')

    R('''distfun <- function(x) dist(x, method = "manhattan")''')
    R('''hclustfun <- function(x) hclust(x, method = "ward.D")''')
    R('''colors <- colorRampPalette(c("blue", "white", "red"))(200)''')

    # read in annotations for sideColors
    R('''anno <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infiles[1])
    R('''anno <- anno[rownames(dat),]''')
    R('''anno <- anno[order(anno$status),]''')

    R('''sidecols <- ifelse(anno$status == "up.RNA", "purple", NA)''')
    R('''sidecols <- ifelse(anno$status == "down.RNA", "brown", sidecols)''')
    R('''sidecols <- ifelse(anno$status == "up.both", "red", sidecols)''')
    R('''sidecols <- ifelse(anno$status == "down.both", "darkGreen", sidecols)''')
    R('''sidecols <- ifelse(anno$status == "down.DNA", "blue", sidecols)''')
    R('''sidecols <- ifelse(anno$status == "up.DNA", "orange", sidecols)''')

    # scale data
    R('''dat.s <- data.frame(t(apply(dat,1,scale)))''')
    R('''dat.s <- dat.s[rownames(anno),]''')
    R('''colnames(dat.s) <- colnames(dat)''')

    # just get steady state and colitis
    # R('''print(dim(dat.s))''')
    R('''dat.s <- dat.s[,c(9:16, 25:32)]''')

    R('''png("%s")''' % outfile)
    R('''heatmap.2(as.matrix(dat.s), 
                   trace = "none",
                   Colv = F,
                   Rowv = T,
                   col = colors,
                   distfun = distfun,
                   hclustfun = hclustfun,
                   margins = c(15,15),
                   RowSideColors = sidecols)''')
    R["dev.off"]()

#########################################
#########################################
#########################################
# associate genes with taxonomic 
# assignments
#########################################
#########################################
#########################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".list")
def buildGeneListForTaxaAssociations(infile, outfile):
    '''
    build a list of COGs for use with annotating to
    taxa
    '''
    statement = '''cat %(infile)s | cut -f1 | tail -n+2 > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@follows(mkdir("associate_taxa.dir"), buildGeneListForTaxaAssociations)
@transform(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/*.genes.tsv.gz")),
            regex("(\S+)/(\S+).genes.tsv.gz"),
            add_inputs(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "diamond.dir/*.diamond.lca.gz"))),
            r"associate_taxa.dir/\2.ctaxa.tsv.gz")
def buildCountTaxaInCogsRNA(infiles, outfile):
    '''
    reads that map to differentially expressed genes are cross-referenced
    with their genera assignment - provides per taxa counts
    '''
    job_options="-l mem_free=20G"
    m = PARAMS.get("genes_map")
    alignment_genes = infiles[0]
    track = P.snip(os.path.basename(alignment_genes), ".diamond.genes.tsv.gz")
    alignment_taxa = [x for x in infiles[1] if os.path.basename(x).startswith(track)][0]

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

@merge(buildCountTaxaInCogsRNA, "associate_taxa.dir/associated_taxa_counts.tsv.gz")
def mergeCountTaxaInCogsRNA(infiles, outfile):
    '''
    merge tables of taxa associated with genes
    '''
    pattern = os.path.dirname(infiles[0]) + "/*.ctaxa.tsv.gz"
    prefixes = ",".join([P.snip(os.path.basename(x), ".diamond.ctaxa.tsv.gz") for x in glob.glob(pattern)])
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

@follows(mkdir("cog_diversity.dir"))
@transform(buildCountTaxaInCogsRNA, 
           regex("(\S+)/(\S+).tsv.gz"),
           r"cog_diversity.dir/\2.diversity.tsv")
def buildCogDiversity(infile, outfile):
    '''
    for each COG, estimate species richness in terms
    of expression
    '''
    R('''library("vegan")''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''dat <- dat[grep("unassigned", dat$taxa, invert = T),]''')
    R('''cogs <- unique(dat$cog)''')
    R('''chao.ests <- c()''')
    R('''evennesses <- c()''')
    R('''for (cog in cogs){
              dat2 <- dat[dat$cog == cog,]
              res <- data.frame(estimateR(dat2[,3]))
              chao1.est <- res["S.chao1",]
              chao.ests <- append(chao.ests, chao1.est)
              diversity.est <- diversity(dat2[,3])
              nspecies <- specnumber(dat2[,3])
              evenness <- diversity.est/log(nspecies)
              evennesses <- append(evennesses, evenness)
              }''')
    R('''result <- data.frame("gene" = cogs, "chao1" = chao.ests, "evenness" = evennesses)''')
    R('''write.table(result, file = "%s", quote = F, sep = "\t", row.names = F)''' % outfile)


#########################################
#########################################
#########################################

@merge(buildCogDiversity, "cog_diversity.dir/cog_transcriptional_diversity.tsv.gz")
def mergeCogDiversity(infiles, outfile):
    '''
    merge chao1 diversity estimates
    per COG
    '''
    pattern = os.path.dirname(infiles[0]) + "/*.diversity.tsv"
    prefixes = ",".join([P.snip(os.path.basename(x), ".diversity.tsv") for x in glob.glob(pattern)])
    statement = '''python %(scriptsdir)s/combine_tables.py
                   --glob=%(pattern)s
                   --missing=0
                   --prefixes=%(prefixes)s
                   --columns=1
                   --take=2,3
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@merge([mergeCogDiversity, buildGenesOutsidePredictionInterval],
       "cog_diversity.dir/cog_transcriptional_diversity_annotated.tsv")
def mergeCogDiversityWithAnnotations(infiles, outfile):
    '''
    merge cog diversity estimates with COG
    differential abundance status
    '''
    diversity = pandas.read_table(
        infiles[0],
        compression = "gzip"
        )
    d2 = diversity.filter(like = "chao1", axis = 1)
    d3 = diversity.filter(like = "evenness", axis = 1)
    diversity["average_richness"] = pandas.Series(d2.mean(axis = 1))
    diversity["average_evenness"] = pandas.Series(d3.mean(axis = 1))
    annotations = pandas.read_table(infiles[1])
    dat = pandas.merge(diversity, annotations)
    dat.to_csv(outfile, sep = "\t")

#########################################
#########################################
#########################################

@transform(mergeCogDiversityWithAnnotations, 
           suffix("_diversity_annotated.tsv"), 
           "_richness.pdf")
def plotCogRichnessDistributions(infile, outfile):
    '''
    plot the distribution of max(average) values
    per cluster
    '''
    R('''library(gplots)''')
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infile)

    R('''dat$pi_status <- ifelse(dat$status == "NS", "NS", dat$pi_status)''')
    R('''dat$pi_status[dat$pi_status == ""] <- "other_significant"''')

    #R('''dat$pi_status <- ifelse(dat$pi_status == "", dat$status, dat$pi_status)''')
    # # ks test
    # for x, y in itertools.product([1,2,3,4,5],[1,2,3,4,5]):
    #     outf = "associate_taxa.dir/cluster-%i_cluster-%i.sig" % (x,y)
    #     R('''k <- ks.test(averages$mean[averages$class == %i],averages$mean[averages$class == %i])''' % (x,y))
    #     R('''k <- data.frame("D" = k[[1]], "p-value" = k[[2]])''')
    #     R('''write.table(k, file = "%s", sep = "\t")''' % outf)

    R('''plot1 <- ggplot(dat, aes(x = average_richness, colour = factor(pi_status))) + stat_ecdf(size = 1.2)''')
    R('''plot1 + scale_colour_manual(values = c("cyan3", "darkorchid", "black", "darkgoldenrod2",  "grey", "darkBlue"))''')
#stat_ecdf(size = 0.3)''')
    # R('''plot1 + scale_colour_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@transform(mergeCogDiversityWithAnnotations, 
           suffix("_diversity_annotated.tsv"), 
           "_evenness.pdf")
def plotCogEvennessDistributions(infile, outfile):
    '''
    plot the distribution of max(average) values
    per cluster
    '''
    R('''library(gplots)''')
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infile)
    R('''dat$pi_status <- ifelse(dat$status == "NS", "NS", dat$pi_status)''')
    R('''dat$pi_status[dat$pi_status == ""] <- "other_significant"''')
    R('''dat <- na.omit(dat)''')
    R('''print(unique(dat$pi_status))''')

#    R('''dat$pi_status <- ifelse(dat$status == "NS", "NS", dat$pi_status)''')

    # R('''dat$pi_status <- ifelse(dat$pi_status == "", dat$status, dat$pi_status)''')
    # # ks test
    # for x, y in itertools.product([1,2,3,4,5],[1,2,3,4,5]):
    #     outf = "associate_taxa.dir/cluster-%i_cluster-%i.sig" % (x,y)
    #     R('''k <- ks.test(averages$mean[averages$class == %i],averages$mean[averages$class == %i])''' % (x,y))
    #     R('''k <- data.frame("D" = k[[1]], "p-value" = k[[2]])''')
    #     R('''write.table(k, file = "%s", sep = "\t")''' % outf)

    R('''plot1 <- ggplot(dat, aes(x = average_evenness, colour = factor(pi_status))) + stat_ecdf(size = 1.2)''')
    R('''plot1 + scale_colour_manual(values = c("cyan3", "darkorchid", "black", "darkgoldenrod2",  "grey", "darkBlue"))''')
    # R('''plot1 + scale_colour_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("associate_taxa.dir"), buildGeneListForTaxaAssociations)
@transform(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/*.genes.tsv.gz")),
            regex("(\S+)/(\S+).genes.tsv.gz"),
            add_inputs(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "diamond.dir/*.diamond.lca.gz"))),
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

@merge(buildProportionTaxaInCogsRNA, "associate_taxa.dir/associated_ptaxa.tsv.gz")
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

@transform(mergeProportionTaxaInCogsRNA, 
           suffix(".tsv.gz"), 
           ".average_ptaxa.matrix")
def buildTaxaCogCountsMatrix(infile, outfile):
    '''
    build cog x taxa matrix
    '''
    inf = IOTools.openFile(infile)
    header = inf.readline()
    result = {}

    # create container for results
    for line in inf.readlines():
        data = line[:-1].split("\t")
        cog, taxa = data[0], data[1]
        if taxa == "unassigned": continue
        result[cog] = {}

    # get average % taxa per cog
    inf = IOTools.openFile(infile)
    header = inf.readline()
    for line in inf.readlines():
        data = line[:-1].split("\t")
        
        if len(data) == 19:
            cog, taxa = data[0], data[1]
            values = map(float,data[3:])
        elif len(data) == 20:
            cog, taxa = data[0], data[1]
            values = map(float,data[4:])
        else:
            cog, taxa = data[0], data[1]
            values = map(float,data[2:])
        if taxa == "unassigned": continue
        ave = np.mean(values)
        try:
            result[cog][taxa] = ave
        except KeyError: continue
    df = pandas.DataFrame(result)
    df.to_csv(outfile, sep = "\t", na_rep = 0)

###################################################
###################################################
###################################################

@merge([buildTaxaCogCountsMatrix, buildGenesOutsidePredictionInterval], 
       "associate_taxa.dir/taxa_cogs_matrix_annotated.pdf")
def plotMaxTaxaContribution(infiles, outfile):
    '''
    plot the distribution of maximum genus
    contribution per gene set
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''annotations <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])

    R('''maximums <- apply(dat, 2, max)''')
    R('''dat2 <- data.frame("cog" = colnames(dat), "max" = maximums)''')
    R('''dat3 <- merge(dat2, annotations, by.x = "cog", by.y = "gene")''')
    R('''dat3$pi_status <- ifelse(dat3$status == "NS", "NS", dat3$pi_status)''')
    R('''dat3$pi_status[is.na(dat3$pi_status)] <- "other_significant"''')

    R('''plot1 <- ggplot(dat3, aes(x = as.numeric(as.character(max)), group = pi_status, colour = pi_status)) + stat_ecdf(size = 1.1)''')
    R('''plot1 + scale_colour_manual(values = c("cyan3", "darkorchid", "black", "darkgoldenrod2",  "grey", "darkBlue"))''')
    R('''ggsave("%s")''' % outfile)

###################################################
###################################################
###################################################

@merge([buildTaxaCogCountsMatrix, buildGenesOutsidePredictionInterval], 
       "associate_taxa.dir/taxa_cogs_matrix_annotated.sig")
def testSignificanceOfMaxTaxaContribution(infiles, outfile):
    '''
    plot the distribution of maximum genus
    contribution per gene set
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''annotations <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])

    R('''maximums <- apply(dat, 2, max)''')
    R('''dat2 <- data.frame("cog" = colnames(dat), "max" = maximums)''')
    R('''dat3 <- merge(dat2, annotations, by.x = "cog", by.y = "gene")''')
    R('''dat3$pi_status <- ifelse(dat3$status == "NS", "NS", dat3$pi_status)''')
    R('''diff.up.rna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.up.rna"]))''')
    R('''diff.down.rna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.down.rna"]))''')
    R('''diff.up.dna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.up.dna"]))''')
    R('''diff.down.dna <- as.numeric(as.character(dat3$max[dat3$pi_status == "diff.down.dna"]))''')
    R('''ns <- as.numeric(as.character(dat3$max[dat3$pi_status == "NS"]))''')
    
    # ks tests
    R('''ks1 <- ks.test(diff.up.rna, ns)''')
    R('''ks2 <- ks.test(diff.down.rna, ns)''')

    R('''ks3 <- ks.test(diff.up.dna, ns)''')
    R('''ks4 <- ks.test(diff.down.dna, ns)''')

    R('''res <- data.frame("RNAGreaterThanDNA.up.pvalue" = ks1$p.value,
                           "RNAGreaterThanDNA.up.D" = ks1$statistic,
                           "RNAGreaterThanDNA.down.pvalue" = ks2$p.value,
                           "RNAGreaterThanDNA.down.D" = ks2$statistic,
                           "DNAGreaterThanRNA.up.pvalue" = ks3$p.value,
                           "DNAGreaterThanRNA.up.D" = ks3$statistic,
                           "DNAGreaterThanRNA.down.pvalue" = ks4$p.value,
                           "DNAGreaterThanRNA.down.D" = ks4$statistic)''')
    R('''write.table(res, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)

#########################################
#########################################
#########################################

@transform(buildTaxaCogCountsMatrix, 
           suffix(".matrix"), 
           add_inputs(buildGenesOutsidePredictionInterval),
           ".pdf")
def heatmapTaxaCogProportionMatrix(infiles, outfile):
    '''
    plot the taxa associated with each cog on
    a heatmap
    '''
    R('''library(gplots)''')
    R('''library(gtools)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infiles[0])
    R('''print (ncol(dat))''')

    R('''annotations <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])
    R('''rownames(annotations) <- annotations$gene''')

    # get genes present in both - not sure why these are different
    # in the first place - need to check
    R('''genes <- intersect(rownames(annotations), colnames(dat))''')
    R('''dat <- dat[, genes]''')
    R('''dat <- dat[grep("unassigned", rownames(dat), invert = T),]''')

    R('''genera <- rownames(dat)''')
    R('''rownames(dat) <- genera''')
    R('''colnames(dat) <- genes''')
    R('''annotations <- annotations[genes,]''')

    R('''annotations <- annotations[order(annotations$pi_status),]''')
    
    # only for the COGs that have RNA fold > DNA fold up-regulated
    R('''annotations <- annotations[annotations$pi_status == "diff.up.rna",]''')

    R('''annotations <- na.omit(annotations)''')
    R('''dat <- dat[,rownames(annotations)]''')

    R('''annotation <- data.frame(cluster = as.character(annotations$pi_status))''')
    R('''rownames(annotation) <- rownames(annotations)''')
 
    R('''colors1 <- c("grey")''')
    R('''names(colors1) <- c("diff.up.rna")''')

    R('''anno_colors <- list(cluster = colors1)''')

    R('''cols <- colorRampPalette(c("white", "darkBlue"))(150)''')

    R('''dat <- dat[,colSums(dat > 50) >= 1]''')
    R('''dat <- dat[rowSums(dat > 10) >= 1,]''')

    R('''pdf("%s", height = 10, width = 15)''' % outfile)
    R('''library(pheatmap)''')
    R('''pheatmap(dat, 
                  clustering_distance_cols = "manhattan",
                  clustering_method = "ward",
                  annotation = annotation,
                  annotation_colors = anno_colors,
                  cluster_rows = T,
                  cluster_cols = F,
                  color = cols,
                  fontsize = 8)''')
    R["dev.off"]()

###################################################
###################################################
###################################################

@follows(mergeCountTaxaInCogsRNA, 
         mergeProportionTaxaInCogsRNA)
def associate_taxa_rna():
    pass

#########################################
#########################################
#########################################

@follows(mkdir("associate_taxa.dir"), buildGeneListForTaxaAssociations)
@transform(glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/*.genes.tsv.gz")),
            regex("(\S+)/(\S+).genes.tsv.gz"),
            add_inputs(glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "diamond.dir/*.diamond.lca.gz"))),
            r"associate_taxa.dir/\2.dna_ctaxa.tsv.gz")
def buildCountTaxaInCogsDNA(infiles, outfile):
    '''
    reads that map to differentially expressed genes are cross-referenced
    with their genera assignment - provides per taxa counts
    '''
    job_options="-l mem_free=20G"
    m = PARAMS.get("genes_map")
    alignment_genes = infiles[0]
    track = P.snip(os.path.basename(alignment_genes), ".diamond.genes.tsv.gz")
    alignment_taxa = [x for x in infiles[1] if os.path.basename(x).startswith(track)][0]

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
    prefixes = ",".join([P.snip(os.path.basename(x), ".diamond.dna_ctaxa.tsv.gz") for x in glob.glob(pattern)])
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

@follows(mkdir("scatterplot.dir"))
@split([runMetagenomeSeqPerCOGAndTaxa, 
        os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.diff.tsv"),
        os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.diff.tsv")],
       "scatterplot.dir/*scatters.pdf")
def plotPerCogTaxaDNAFoldRNAFold(infiles, outfiles):
    '''
    plot per COG ratio
    '''
    rna, dna, rna_cog, dna_cog = infiles

    R('''library(ggplot2)''')

    # read in cogs + taxa
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''dna <- dna[dna$group2 == "WT" & dna$group1 == "HhaIL10R",]''')
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''rna <- rna[rna$group2 == "WT" & rna$group1 == "HhaIL10R",]''')

    # read in cogs alone
    R('''dna.cog <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna_cog)
    R('''dna.cog <- dna.cog[dna.cog$group2 == "WT" & dna.cog$group1 == "HhaIL10R",]''')
    R('''rna.cog <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna_cog)
    R('''rna.cog <- rna.cog[rna.cog$group2 == "WT" & rna.cog$group1 == "HhaIL10R",]''')

    # merge data
    R('''dat <- merge(dna, rna, by.x = "taxa", by.y = "taxa", all.x = T, all.y = T, suffixes = c("dna", "rna"))''')

       # add gene column
    R('''dat$gene <- unlist(strsplit(dat$taxa, "-"))[seq(1, nrow(dat)*2, 2)]''')
    R('''dat$genus <- unlist(strsplit(dat$taxa, "-"))[seq(2, nrow(dat)*2, 2)]''')

    # more merging
    R('''dat <- merge(dat, dna.cog, by.x = "gene", by.y = "taxa", suffixes = c("taxa","dna.cog"), all.x = T, all.y = T)''')
    R('''dat <- merge(dat, rna.cog, by.x = "gene", by.y = "taxa", suffixes = c("dna.cog", "rna.cog"), all.x = T, all.y = T)''')

    # sub NA for 0
    R('''dat[is.na(dat)] <- 0''')

    # NOTE these are specified and hardcoded 
    # here
    R('''cogs <- c("COG0783", "COG2837", "COG0435","COG5520")''')

    # iterate over cogs and scatterplot
    # fold changes in DNA and RNA analysis.
    # if not present in one or other then fold change will
    # be 0
    R('''for (cog in cogs){
            dat2 <- dat[grep(cog, dat$gene),]
             
            # add the data for COG fold changes and abundance
            dat3 <- data.frame("genus" = append(dat2$genus, cog),
                               "dna.fold" = append(dat2$logFCdna, unique(dat2$logFCdna.cog)),
                               "rna.fold" = append(dat2$logFCrna, unique(dat2$logFCrna.cog)),
                               "abundance" = append(dat2$AveExprdna, sum(dat2$AveExprdna.cog)))
                               
            suffix <- paste(cog, "scatters.pdf", sep = ".")
            outname <- paste("scatterplot.dir", suffix, sep = "/")

            plot1 <- ggplot(dat3, aes(x = dna.fold, y = rna.fold, size = log10(abundance), label = genus))
            plot2 <- plot1 + geom_point(shape = 18) + geom_text(hjust = 0.5, vjust = 1, size = 5)
           
            plot3 <- plot2 + geom_abline(intercept = 0, slope = 1, colour = "blue") + geom_hline(yintercept = c(-1,1), linetype = "dashed")
            plot4 <- plot3 + geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
 #           plot4 
            ggsave(outname)
            }''')

###################################################
###################################################
###################################################

@transform(buildRNADNAMergedMatrix,
           regex("(\S+)/(\S+).tsv"),
           r"enrichment.dir/*heatmap.pdf")
def heatmapFunctions(infiles, outfile):
    '''
    heatmap per cogs per functional category
    '''
    print infiles

    # funcs, mat = infiles
    # functions = set()
    # for line in IOTools.openFile(funcs):
    #     data = line[:-1].split("\t")
    #     f = data[11]
    #     if f == "NA": continue
    #     functions.add(data[11])

    # R('''library(gplots)''')
    # R('''mat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % mat)
    # R('''rownames(mat) <- mat$gene''')
    # R('''mat <- mat[,1:ncol(mat)-1]''')
    # R('''funcs <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % funcs)
    # R('''sets <- c("down_dna", "down_rna", "up_dna", "up_rna", "up_both", "down_both", "up_all")''')
    # R('''colors <- colorRampPalette(c("blue", "white", "red"))(75)''')
    # for f in functions:
    #     R('''for (i in 2:8){
    #               setname <- paste(colnames(funcs)[i], "%s", sep = "_")
    #               setname <- gsub(" ", "_", setname)
    #               setname <- gsub("/", "_", setname)
    #               outname <- paste("enrichment.dir", setname, sep = "/")
    #               outname <- paste(outname, ".heatmap.pdf", sep = "")
    #               cogs <- funcs$gene[funcs[,i] == 1 & funcs$V4 == "%s"]
    #               cogs <- na.omit(cogs)
    #               if (length(cogs) <= 1){
    #                   next}                  
    #               mat2 <- mat[cogs,]
    #               mat2.s <- data.frame(t(apply(mat2,1,scale)))
    #               colnames(mat2.s) <- colnames(mat2)
    #               mat2.s <- mat2.s[, c(9:16, 25:32)]
    #               pdf(outname)
    #               heatmap.2(as.matrix(mat2.s), 
    #               trace = "none",
    #               Colv = F,
    #               Rowv = T,
    #               col = colors,
    #               #distfun = distfun,
    #               #hclustfun = hclustfun,
    #               margins = c(15,15))
    #               dev.off()}''' % (f, f))


###################################################
###################################################
###################################################

@follows(mkdir("pval_dist.dir"))
@split(annotateRNADNARatio, "pval_dist.dir/*.pvals.pdf")
def plotPvalDistributions(infile, outfiles):
    '''
    plot the pvalue distributions for rna and dna seq
    and where those unique COGs fall in thata distribution
    '''
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)

    # for DNA 
    R('''pdf("pval_dist.dir/dna.pvals.pdf")''')
    R('''hist(-log10(dat$pdna), breaks = 25)''')
    R('''abline(v = 1.3, lty = 2)''')
    R('''xup = -log10(dat$pdna[dat$status == "up.RNA"])''')
    R('''y0up = rep(0, length(xup))''')
    R('''y1up = rep(100, length(xup))''')
    R('''segments(x0 = xup, y0 = y0up, x1 = xup, y1 = y1up, col = "purple")''')
    R('''xdown = -log10(dat$pdna[dat$status == "down.RNA"])''')
    R('''y0down = rep(0, length(xdown))''')
    R('''y1down = rep(100, length(xdown))''')
    R('''segments(x0 = xdown, y0 = y0down, x1 = xdown, y1 = y1down, col = "brown")''')
    R["dev.off"]()

    # for RNA 
    R('''pdf("pval_dist.dir/rna.pvals.pdf")''')
    R('''hist(-log10(dat$prna), breaks = 25)''')
    R('''abline(v = 1.3, lty = 2)''')
    R('''xup = -log10(dat$prna[dat$status == "up.DNA"])''')
    R('''y0up = rep(0, length(xup))''')
    R('''y1up = rep(100, length(xup))''')
    R('''segments(x0 = xup, y0 = y0up, x1 = xup, y1 = y1up, col = "orange")''')
    R('''xdown = -log10(dat$prna[dat$status == "down.DNA"])''')
    R('''y0down = rep(0, length(xdown))''')
    R('''y1down = rep(100, length(xdown))''')
    R('''segments(x0 = xdown, y0 = y0down, x1 = xdown, y1 = y1down, col = "blue")''')
    R["dev.off"]()


###################################################
###################################################
###################################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".scatterplot.png")
def scatterplotFoldChanges(infile, outfile):
    '''
    scatterplot the fold changes bewteen RNA and DNA
    analysis - includes sized as -log10(p-values)
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''plot1 <- ggplot(dat, aes(x = dna, y = rna, size = -log10(prna)))''')
    R('''lm1 <- lm(dat$rna ~ dat$dna)''')
    R('''i <- lm1[[1]][[1]]''')
    R('''s <- lm1[[1]][[2]]''')
    R('''plot2 <- plot1 + geom_point()''')
    R('''plot3 <- plot2 + geom_abline(intercept = i, slope = s, linetype = "dashed")''')
    R('''plot4 <- plot3 + geom_abline(intercept = 0, slope = 1) + xlim(c(-5, 5)) + ylim(c(-5,5))''')
    R('''plot4''')
    R('''ggsave("%s")''' % outfile)

###################################################
###################################################
###################################################

@transform(annotateRNADNARatio, suffix(".tsv"), ".classes.tsv")
def mixtureModel(infile, outfile):
    '''
    run finite mixture models to cluster 
    genes based on their changes in DNA
    and RNA analysis. output classifications
    '''
    R('''library(mclust)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    # just get those signficant in either DNA or RNA or both
    R('''dat2 <- dat[dat$status != "NS",]''')
    R('''genes <- dat2$gene''')
    # get fold changes and p-values - these are adjusted p-values
    R('''dat2 <- dat2[, c("rna", "dna", "prna", "pdna")]''')

    # take -log10(p-values)
    R('''dat2$prna <- -log10(dat2$prna)''')
    R('''dat2$pdna <- -log10(dat2$pdna)''')
    R('''rownames(dat2) <- genes''')
    R('''clust <- Mclust(dat2)''')
    R('''classification <- data.frame(clust$classification)''')
    R('''colnames(classification) <- "class"''')
    R('''dat2$class <- ifelse(rownames(dat2) %in% rownames(classification), classification$class, NA)''')
    R('''dat2$gene <- rownames(dat2)''')
    R('''write.table(dat2, file = "%s", sep = "\t", row.names = F)''' % outfile)

###################################################
###################################################
###################################################

@transform(mixtureModel, suffix(".tsv"), ".png")
def plotClasses(infile, outfile):
    '''
    plot the classification of genes on a scatterplot
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''plot1 <- ggplot(dat, aes(x = dna, y = rna, size = prna, alpha = 0.5, colour = factor(class)))''')

    # regression line
    R('''lm1 <- lm(dat$rna ~ dat$dna)''')
    R('''i <- lm1[[1]][[1]]''')
    R('''s <- lm1[[1]][[2]]''')
    R('''plot2 <- plot1 + geom_point() + stat_density2d() + geom_abline(intercept = i, slope = s, linetype = "dashed")''')
    R('''plot3 <- plot2 + geom_vline(xintercept = c(-1,1), linetype = "dashed", colour = "darkGrey")''')
    R('''plot4 <- plot3 + geom_hline(yintercept = c(-1,1), linetype = "dashed", colour = "darkGrey")''')
    R('''plot5 <- plot4 + geom_abline(intercept = 0, slope = 1) + xlim(c(-5, 5)) + ylim(c(-5,5))''')
    R('''plot5''') #+ scale_colour_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''ggsave("%s")''' % outfile)


@follows(annotateFunctions)
def annotations():
    pass

###################################################
###################################################
###################################################
# Compare abundance estimations 
###################################################
###################################################
###################################################

@follows(mkdir("compare_abundance.dir"))
@split(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/*diamond*.norm.matrix")) 
       + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/*diamond*.norm.matrix")),
       "compare_abundance.dir/*.pdf")
def scatterplotAbundanceEstimates(infiles, outfiles):
    '''
    scatterplot abundance estimates for each sample
    '''
    levels = ["phylum", "class", "order", "family", "genus", "species"]
    for level in levels:
        rna, dna = [inf for inf in infiles if inf.find(level) != -1]
        R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
        R('''rownames(rna) <- rna$taxa''')
        R('''rna <- rna[,1:ncol(rna)-1]''')
        R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
        R('''rownames(dna) <- dna$taxa''')
        R('''dna <- dna[,1:ncol(dna)-1]''')

        # intersection of species present
        R('''keep <- intersect(rownames(rna), rownames(dna))''')

        # get data where there is rna and dna
        R('''rna <- rna[keep,]''')
        R('''dna <- dna[keep,]''')
        
        # take averages
        R('''rna.ave <- data.frame(apply(rna, 1, mean))''')
        R('''dna.ave <- data.frame(apply(dna, 1, mean))''')

        R('''print(cor(dna.ave,rna.ave))''')
        R('''png("compare_abundance.dir/average_abundance.%s.png")''' % level)
        R('''plot(dna.ave[,1], 
                  rna.ave[,1], 
                  pch = 16, 
                  col = "slateGrey",
                  xlab = "Mean DNA abundance",
                  ylab = "Mean RNA abundance",
                  main = paste("N = ", nrow(dna.ave), sep = ""))
             abline(lm(rna[,1]~dna[,1], na.rm = T))''')
        R["dev.off"]()


        R('''for (i in 1:ncol(dna)){
             name = paste("compare_abundance.dir", paste(colnames(dna[i]), ".%s.pdf", sep = ""), sep = "/")
             pdf(name)
             plot(dna[,i], rna[,i], 
                  pch = 16,
                  col = "slateGrey",
                  xlab = "DNA normalised counts",
                  ylab = "RNA normalised counts")
             abline(lm(rna[,i]~dna[,i], na.rm = T), lty = 2)
          # fold <- rna[,i]-dna[,i]
          # points(dna[,i][fold > 3], rna[,i][fold > 3], pch = 16, col = "blue")   
          # legend("bottomright", legend = rownames(dna[fold > 3,]), cex = 0.8)
          dev.off()
   
         }
     ''' % level)


###################################################
###################################################
###################################################
@follows(mkdir("compare_gene_abundance.dir"))
@split(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/*.norm.matrix")) 
       + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/*.norm.matrix")),
       "compare_gene_abundance.dir/*genes*.pdf")
def scatterplotGeneAbundanceEstimates(infiles, outfiles):
    '''
    scatterplot abundance estimates for each sample
    '''
    rna, dna = infiles
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''rownames(rna) <- rna$taxa''')
    R('''rna <- rna[,1:ncol(rna)-1]''')
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''rownames(dna) <- dna$taxa''')
    R('''dna <- dna[,1:ncol(dna)-1]''')

    # intersection of species genes
    R('''keep <- intersect(rownames(rna), rownames(dna))''')

    # get RNA where there is DNA
    # add labels where the RNA is higher
    # than the DNA by > 2 fold
    R('''rna <- rna[keep,]''')
    R('''dna <- dna[keep,]''')

    # take averages
    R('''rna.ave <- data.frame(apply(rna, 1, mean))''')
    R('''dna.ave <- data.frame(apply(dna, 1, mean))''')
    R('''print(cor(dna.ave,rna.ave))''')
    R('''png("compare_gene_abundance.dir/average_abundance.genes.png")''')
    R('''plot(dna.ave[,1], 
                  rna.ave[,1], 
                  pch = 16, 
                  col = "slateGrey",
                  xlab = "Mean DNA abundance",
                  ylab = "Mean RNA abundance",
                  main = paste("N = ", nrow(dna), sep = ""))
             abline(lm(rna[,1]~dna[,1], na.rm = T))''')
    R["dev.off"]()

    # histogram rna/dna ratio
    R('''pdf("compare_gene_abundance.dir/average_abundance.genes.ratio.pdf")''')
    R('''hist(log2(rna.ave[,1]/dna.ave[,1]), 
                  xlab = "log2(RNA/DNA)", breaks = 25)''')
    R('''abline(v = 0, lty = 2, lwd = 2)''')
    R["dev.off"]()


    R('''for (i in 1:ncol(dna)){
             name = paste("compare_gene_abundance.dir", paste(colnames(dna[i]), ".genes.pdf", sep = ""), sep = "/")
             pdf(name)
             plot(dna[,i], rna[,i], 
                  pch = 16,
                  col = "slateGrey",
                  xlab = "DNA normalised counts",
                  ylab = "RNA normalised counts")
             fold <- rna[,i]-dna[,i]
             abline(lm(rna[,i]~dna[,i], na.rm = T), lty = 2)
             # points(dna[,i][fold > 3], rna[,i][fold > 3], pch = 16, col = "blue")   
             # legend("bottomright", legend = rownames(dna[fold > 3,]), cex = 0.8)

          dev.off()
   
         }
     ''')

###################################################
###################################################
###################################################
@follows(mkdir("compare_detected.dir"))
@split(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/*diamond*.aggregated.counts.tsv.gz")) 
       + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/*diamond*.aggregated.counts.tsv.gz")),
       "compare_detected.dir/*.overlap.tsv")
def buildTaxaDetectionOverlap(infiles, outfiles):
    '''
    build species detection overlap
    '''
    levels = ["phylum", "class", "order", "family", "genus", "species"]
    for level in levels:
        rna, dna = [inf for inf in infiles if inf.find(level) != -1]
        R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
        R('''rownames(rna) <- rna$taxa''')
        R('''rna <- rna[,1:ncol(rna)]''')
        R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
        R('''rownames(dna) <- dna$taxa''')
        R('''dna <- dna[,1:ncol(dna)]''')

        R('''for (i in 1:ncol(dna)){
             name = paste("compare_detected.dir", paste(colnames(dna[i]), ".%s.overlap.tsv", sep = ""), sep = "/")
             taxa.rna = rownames(rna[rna[,i] !=0,])
             taxa.dna = rownames(dna[dna[,i] !=0,])
             nrna = length(taxa.rna)
             ndna = length(taxa.dna)
             noverlap = length(intersect(taxa.rna, taxa.dna))
             result = data.frame(nrna = nrna, ndna = ndna, noverlap = noverlap)
             write.table(result, file = name, sep = "\t", row.names = F)}
        ''' % level)

###################################################
###################################################
###################################################

@follows(mkdir("compare_detected.dir"))
@split(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "counts.dir/*diamond*.aggregated.counts.tsv.gz")) 
       + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "counts.dir/*diamond*.aggregated.counts.tsv.gz")),
       "compare_detected.dir/*.abundance.pdf")
def plotAbundanceLevelsOfTaxaOverlap(infiles, outfiles):
    '''
    build species detection overlap
    '''
    levels = ["genus"]
    for level in levels:
        rna, dna = [inf for inf in infiles if inf.find(level) != -1]
        outfile = "compare_detected.dir/genus.overlap.abundance.pdf"

        R('''library(ggplot2)''')

        # get rna reads per million
        R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
        R('''rownames(rna) <- rna$taxa''')
        R('''rna <- rna[,2:ncol(rna)]''')
        R('''rna <- sweep(rna, 2, colSums(rna)/1000000, "/")''')

        # get dna reads per million
        R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
        R('''rownames(dna) <- dna$taxa''')
        R('''dna <- dna[,2:ncol(dna)]''')
        R('''dna <- sweep(dna, 2, colSums(dna)/1000000, "/")''')

        # common and distinct sets
        R('''common <- intersect(rownames(dna), rownames(rna))''')
        R('''rna.only <- setdiff(rownames(rna), rownames(dna))''')
        R('''dna.only <- setdiff(rownames(dna), rownames(rna))''')

        # boxplot the abundance levels
        R('''rna.common <- apply(rna[common,], 1, mean)''')
        R('''dna.common <- apply(dna[common,], 1, mean)''')
        R('''rna.distinct <- apply(rna[rna.only,], 1, mean)''')
        R('''dna.distinct <- apply(dna[dna.only,], 1, mean)''')

        # test sig bewteen groups
        R('''ttest1 <- wilcox.test(rna.common, rna.distinct)''')
        R('''ttest2 <- wilcox.test(dna.common, dna.distinct)''')
        R('''ttest3 <- wilcox.test(rna.common, dna.distinct)''')
        R('''ttest4 <- wilcox.test(dna.common, rna.distinct)''')
        R('''ttest5 <- wilcox.test(dna.common, rna.common)''')

        R('''res <- data.frame("rna.common_vs_rna.distinct" = ttest1$p.value,
                               "dna.common_vs_dna.distinct" = ttest2$p.value,
                               "rna.common_vs_dna.distinct" = ttest3$p.value,
                               "dna.common_vs_rna.distinct" = ttest4$p.value,
                               "dna.common_vs_rna.common" = ttest5$p.value)''')
        outname_sig = P.snip(outfile, ".pdf") + ".sig"
        R('''write.table(res, file = "%s", row.names = F, sep = "\t", quote = F)''' % outname_sig)

        # create dataframe for plotting
        R('''dat <- data.frame(values = c(dna.distinct, dna.common, rna.common, rna.distinct),
                               status = c(rep("unique.dna", length(dna.distinct)),
                                        rep("common.dna", length(dna.common)),
                                        rep("common.rna", length(rna.common)),
                                        rep("unique.rna", length(rna.distinct))))''')
        R('''ggplot(dat, aes(x = factor(status, levels = status), y = values, stat = "identity")) + geom_boxplot() + scale_y_log10()''')
        R('''ggsave("%s")''' % outfile)

###################################################
###################################################
###################################################
@follows(mkdir("compare_detected_genes.dir"))
@split(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.tsv.gz")) 
       + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.tsv.gz")),
       "compare_detected_genes.dir/*.overlap.tsv")
def buildGeneDetectionOverlap(infiles, outfiles):
    '''
    build species detection overlap
    '''
    rna, dna = infiles
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''rownames(rna) <- rna$taxa''')
    R('''rna <- rna[,1:ncol(rna)-1]''')
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''rownames(dna) <- dna$taxa''')
    R('''dna <- dna[,1:ncol(dna)-1]''')
        
    R('''for (i in 1:ncol(dna)){
             name = paste("compare_detected_genes.dir", paste(colnames(dna[i]), ".overlap.tsv", sep = ""), sep = "/")
             taxa.rna = rownames(rna)[rna[,i] !=0]
             taxa.dna = rownames(dna)[dna[,i] !=0]
             nrna = length(taxa.rna)
             ndna = length(taxa.dna)
             noverlap = length(intersect(taxa.rna, taxa.dna))
             result = data.frame(nrna = nrna, ndna = ndna, noverlap = noverlap)
             write.table(result, file = name, sep = "\t", row.names = F)}
        ''')

###################################################
###################################################
###################################################
@follows(mkdir("compare_detected_genes.dir"))
@merge(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.tsv.gz")) 
       + glob.glob(os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.tsv.gz")),
       "compare_detected_genes.dir/genes.abundance.pdf")
def plotAbundanceLevelsOfGeneOverlap(infiles, outfile):
    '''
    plot abundances for unique and common genes
    '''
    rna, dna = infiles

    R('''library(ggplot2)''')

    # get rna reads per million
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''rownames(rna) <- rna$taxa''')
    R('''rna <- rna[,2:ncol(rna)]''')
    R('''rna <- sweep(rna, 2, colSums(rna)/1000000, "/")''')

    # get dna reads per million
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''rownames(dna) <- dna$taxa''')
    R('''dna <- dna[,2:ncol(dna)]''')
    R('''dna <- sweep(dna, 2, colSums(dna)/1000000, "/")''')

    # common and distinct sets
    R('''common <- intersect(rownames(dna), rownames(rna))''')
    R('''dna.only <- setdiff(rownames(dna), rownames(rna))''')

    # boxplot the abundance levels
    R('''rna.common <- apply(rna[common,], 1, mean)''')
    R('''dna.common <- apply(dna[common,], 1, mean)''')
    R('''dna.distinct <- apply(dna[dna.only,], 1, mean)''')

    # test sig bewteen groups
    R('''ttest1 <- wilcox.test(dna.common, dna.distinct)''')
    R('''ttest2 <- wilcox.test(rna.common, dna.distinct)''')
    R('''ttest3 <- wilcox.test(dna.common, rna.common)''')

    R('''res <- data.frame("dna.common_vs_dna.distinct" = ttest1$p.value,
                           "rna.common_vs_dna.distinct" = ttest2$p.value,
                           "dna.common_vs_rna.common" = ttest3$p.value)''')
    outname_sig = P.snip(outfile, ".pdf") + ".sig"
    R('''write.table(res, file = "%s", row.names = F, sep = "\t", quote = F)''' % outname_sig)

    # create dataframe for plotting
    R('''dat <- data.frame(values = c(dna.distinct, dna.common, rna.common),
                               status = c(rep("unique.dna", length(dna.distinct)),
                                        rep("common.dna", length(dna.common)),
                                        rep("common.rna", length(rna.common))))''')
                                        
    R('''ggplot(dat, aes(x = factor(status, levels = status), y = values, stat = "identity")) + geom_boxplot() + scale_y_log10()''')
    R('''ggsave("%s")''' % outfile)


###################################################
###################################################
###################################################

@follows(mkdir("compare_detected_genes.dir"))
@merge([os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.tsv.gz"),
        os.path.join(PARAMS.get("dna_communitiesdir"), "genes.dir/gene_counts.tsv.gz"),
        buildGeneDiffList],
       "compare_detected_genes.dir/diff_genes.abundance.pdf")
def plotAbundanceLevelsOfDiffGeneOverlap(infiles, outfile):
    '''
    plot abundances for unique and common genes
    '''
    rna, dna, diff_dna, diff_rna = infiles

    R('''library(ggplot2)''')

    # get rna reads per million
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % rna)
    R('''rna <- sweep(rna, 2, colSums(rna)/1000000, "/")''')

    # get dna reads per million
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % dna)
    R('''dna <- sweep(dna, 2, colSums(dna)/1000000, "/")''')

    # get diff sets
    R('''diff.dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % diff_dna)
    R('''diff.dna <- diff.dna[,1]''')
    R('''diff.rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % diff_rna)
    R('''diff.rna <- diff.rna[,1]''')

    # unique dna
    R('''unique.dna <- setdiff(diff.dna, diff.rna)''')
    
    # unique rna
    R('''unique.rna <- setdiff(diff.rna, diff.dna)''')

    # common
    R('''common <- intersect(diff.dna, diff.rna)''')

    # boxplot the abundance levels
    R('''rna.unique.dna <- apply(rna[unique.dna,], 1, mean)''')
    R('''dna.unique.dna <- apply(dna[unique.dna,], 1, mean)''')

    R('''rna.unique.rna <- apply(rna[unique.rna,], 1, mean)''')
    R('''dna.unique.rna <- apply(dna[unique.rna,], 1, mean)''')

    R('''rna.common <- apply(rna[common,], 1, mean)''')
    R('''dna.common <- apply(dna[common,], 1, mean)''')

    # test sig between groups
    R('''ttest1 <- wilcox.test(rna.unique.dna, dna.unique.dna)''')
    R('''ttest2 <- wilcox.test(rna.unique.rna, rna.unique.rna)''')
    R('''ttest3 <- wilcox.test(rna.common, dna.common)''')

    R('''res <- data.frame("rna.uniqueToDna_vs_dna.uniqueToDna" = ttest1$p.value,
                           "rna.uniqueToRna_vs_dna.uniqueToRna" = ttest2$p.value,
                           "rna.common_vs_dna.common" = ttest3$p.value)''')
    outname_sig = P.snip(outfile, ".pdf") + ".sig"
    R('''write.table(res, file = "%s", row.names = F, sep = "\t", quote = F)''' % outname_sig)

    # create dataframe for plotting
    R('''dat <- data.frame(values = c(rna.unique.dna, dna.unique.dna, 
                                      rna.common, dna.common,
                                      rna.unique.rna, dna.unique.rna),
                           status = c(rep("rna.unique.dna", length(rna.unique.dna)),
                                      rep("dna.unique.dna", length(dna.unique.dna)),
                                      rep("rna.common", length(rna.common)),
                                      rep("dna.common", length(dna.common)),
                                      rep("rna.unique.rna", length(dna.unique.rna)),
                                      rep("dna.unique.rna", length(dna.unique.rna))))''')
                                        
    R('''ggplot(dat, aes(x = factor(status, levels = status), y = values, stat = "identity")) + geom_boxplot() + scale_y_log10()''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("genera_profiles.dir"))
@split(runMetagenomeSeqPerCOGAndTaxa, "genera_profiles.dir/*.tsv")
def buildGenusProfiles(infiles, outfiles):
    '''
    get the normalised count matrices for
    specific genera
    '''
    genera = ["Lactobacillus", "Clostridium", "Bacteroides"]
    for infile in infiles:
        inf = infile.replace(".diff.tsv", ".norm.matrix")
        header = open(inf).readline()
        for genus in genera:
            outf = open("genera_profiles.dir/%s_%s.tsv" % (os.path.basename(inf), genus), "w")
            outf.write(header)
            inf1 = open(inf)
            for line in inf1.readlines():
                if genus in line:
                    outf.write(line)
            outf.close()

#########################################
#########################################
#########################################

@transform(buildGenusProfiles, suffix(".tsv"), ".pdf")
def heatmapGenusProfiles(infile, outfile):
    '''
    draw a heatmap of genus COG profiles
    '''
    R('''library(gtools)''')
    R('''library(gplots)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 17)''' % infile)

    R('''cols = colorRampPalette(c("darkGreen", "black", "purple"))(75)''')
    R('''dat <- dat[, mixedsort(colnames(dat))]''')
    R('''dat.s <- data.frame(t(apply(dat, 1, scale)))''')
    R('''colnames(dat.s) <- colnames(dat)''')
    R('''pdf("%s", height = 12, width = 10)''' % outfile)
    R('''heatmap.2(as.matrix(dat.s), Colv = F, scale = "none", trace = "none", margins = c(15,15), col = cols)''')
    R["dev.off"]()

#########################################
#########################################
#########################################

@merge(runMetagenomeSeqPerCOGAndTaxa, "taxa_cogs_diff.dir/cogs_taxa.tsv")
def getCommonCogsTaxa(infiles, outfile):
    '''
    output a list of common COG-Taxa pairs
    between RNA and DNA analysis
    '''
    inf1, inf2 = infiles
    set1 = set([x[:-1].split("\t")[-1].replace('"', '') for x in open(inf1)])
    set2 = set([x[:-1].split("\t")[-1].replace('"', '') for x in open(inf2)])

    outf = open(outfile, "w")
    for cog_taxa in set1.intersection(set2):
        if cog_taxa == "taxa": continue
        outf.write("%s\n" % cog_taxa)
    outf.close()

#########################################
#########################################
#########################################

@merge([runMetagenomeSeqPerCOGAndTaxa, getCommonCogsTaxa], "taxa_cogs_diff.dir/rna_dna_ratio.tsv")
def getRNADNARatioPerCogTaxa(infiles, outfile):
    '''
    for each COG-genus pair, get the fold change ratio
    between DNA and RNA
    '''
    rna = infiles[0]
    dna = infiles[1]
    common = infiles[2]
    R('''rna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % rna)
    R('''rna <- rna[rna$group1 == "HhaIL10R" & rna$group2 == "WT",]''')
    R('''dna <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % dna)
    R('''dna <- dna[dna$group1 == "HhaIL10R" & dna$group2 == "WT",]''')
    R('''common <- read.csv("%s", header = F, stringsAsFactors = F, sep = "\t")''' % common)

    R('''rownames(rna) <- rna$taxa''')
    R('''rownames(dna) <- dna$taxa''')

    R('''rna <- rna[common[,1],]''')
    R('''dna <- dna[common[,1],]''')

    R('''rna$fold <- 2^rna$logFC''')
    R('''dna$fold <- 2^dna$logFC''')

    R('''ratio <- rna$fold/dna$fold''')
    R('''dat <- data.frame("id" = rownames(rna), 
                            "cogtaxaRNAfold" = rna$logFC,
                            "cogtaxaDNAfold" = dna$logFC,
                            "RNAOverDNA" = ratio,
                            "abundDNA" = dna$AveExpr,
                            "abundRNA" = rna$AveExpr)''')
    R('''write.table(dat, file = "%s", row.names = F, sep = "\t")''' % outfile)
    
#########################################
#########################################
#########################################

@follows(mkdir("classes_ratio.dir"))
@merge([getRNADNARatioPerCogTaxa, mixtureModel], "classes_ratio.dir/classes_ratio.tsv")
def mergeClassesAndRatios(infiles, outfile):
    '''
    merge the gene classes and RNA/DNA fold ratios
    '''
    ratio, classes = infiles
    temp = P.getTempFile(".")

    # separate cogs from genus
    for line in open(ratio):
        if line.startswith('"id"'): continue
        data = line[:-1].split("\t")
        data[0] = data[0].replace("-", "\t").replace('"','')
        temp.write("\t".join(data) + "\n")
    temp.close()
    temp = temp.name

    # merge tables
    R('''ratios <- read.csv("%s", header = F, stringsAsFactors = F, sep = "\t")''' % temp)
    R('''classes <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % classes)
    R('''merged <- merge(classes, ratios, by.x = "gene",  by.y = "V1" , all.x = T, all.y = T)''')
    R('''colnames(merged) <- c("gene", "rna", "dna", "prna", "pdna", "class", "genus", 
                                "cogtaxaRNAfold", "cogtaxaDNAfold", "foldRatio", "abundDNA", "abundRNA")''')
    R('''write.table(merged, file = "%s", sep = "\t", quote = F, row.names = F)''' % outfile)
    os.unlink(temp)

#########################################
#########################################
#########################################

@transform(mergeClassesAndRatios, suffix(".tsv"), ".density.pdf")
def plotRatioByClassDensityPerGenus(infile, outfile):
    '''
    plot the histogram of ratios for each taxa-cog pair
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''p <- ggplot(dat, aes(x = log2(foldRatio), colour = "1", fill = factor(class), alpha = 0.5))''')
    R('''p1 <- p + geom_density() + scale_fill_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''p1 + geom_vline(xintercept = c(-1,1), linetype = "dashed") + scale_colour_manual(values = "grey") + scale_x_continuous(limits = c(-5,5))''') 
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@transform(mixtureModel, regex("(\S+)/(\S+).tsv"), r"classes_ratio.dir/\2.density.pdf")
def plotRatioByClassDensityPerCog(infile, outfile):
    '''
    plot the histogram of ratios for each COG
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    
    # untransform the data and get the ratio of RNA/DNA
    R('''dat$ratio <- 2^dat$rna/2^dat$dna''')
    R('''p <- ggplot(dat, aes(x = log2(ratio), colour = "1", fill = factor(class), alpha = 0.5))''')
    R('''p1 <- p + geom_density() + scale_fill_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''p1 + geom_vline(xintercept = c(-1,1), linetype = "dashed") + scale_colour_manual(values = "grey") + scale_x_continuous(limits = c(-5,5))''')
    R('''ggsave("%s")''' % outfile)


#########################################
#########################################
#########################################

@transform(getRNADNARatioPerCogTaxa, suffix(".tsv"), ".pdf")
def plotPerCOGRatio(infile, outfile):
    '''
    plot per COG ratio
    '''
    R('''library(gtools)''')
    R('''library(ggplot2)''')

    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$id''')
    R('''cogs <- unlist(strsplit(rownames(dat), "-", perl = TRUE))''')
    R('''cogs <- cogs[seq(1,length(cogs),2)]''')

    dirname = os.path.dirname(infile)

    R('''for (cog in cogs){
             outname <- paste("%s", paste(cog, "ratio.pdf", sep = "_"), sep = "/")
             dat2 <- dat[grep(cog, rownames(dat)),]
             dat2 <- dat2[mixedsort(rownames(dat2)),]
             pdf(outname, height = 9, width = 3)
             plot(log2(dat2[,2]), seq(1, nrow(dat2),1), pch = 17, xlim = c(-4,4))
             abline(v = 0)
             abline(v = c(-1,1), lty = 2)
             dev.off()}''' % dirname)
    
#########################################
#########################################
#########################################

@jobs_limit(1,"R")
@transform(runMetagenomeSeqPerCOGAndTaxa, 
           suffix(".diff.tsv"), 
           add_inputs(getCommonCogsTaxa),
           ".pdf")
def plotCandidates(infiles, outfile):
    '''
    barplot differences for different taxa for COGS
    just in WT vs. HhaIL10R
    '''
    infile, test_ids = infiles

    matrix = P.snip(infile, ".diff.tsv") + ".norm.matrix"
    if "dna" in infile:
        cog_prefix = "dna"
    else:
        cog_prefix = "rna"
    R('''library(ggplot2)''')
    R('''library(gtools)''')
    R('''library(plyr)''')
    R('''library(reshape)''')
    R('''candidates <- c("COG5520", "COG0364", "COG0362", "COG0649", "COG1007", "COG2009", "COG1294", "COG3278", "COG3256", "COG2033")''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % matrix)

    # just use common COG-taxa pairs
    R('''test_ids <- read.csv("%s", header = F, stringsAsFactors = F, sep = "\t")''' % test_ids)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[test_ids[,1],]''')

    # hack to get around the "-" present in a rowname - replace
    R('''rownames(dat)[grep("Phi29", rownames(dat))] <- gsub("Phi29-", "Phi29_", rownames(dat)[grep("Phi29", rownames(dat))])''')
    R('''dat <- dat[,1:ncol(dat)-1]''')
    R('''cogs <- unlist(strsplit(rownames(dat), "-", perl = TRUE))''')
    R('''cogs <- cogs[seq(1,length(cogs),2)]''')

    R('''dat$cog <- cogs''')
    R('''for (cog in cogs){
             outname <- paste(paste("taxa_cogs_diff.dir", paste("%s", cog, sep = "_"), sep = "/"), "pdf", sep = ".") 
             dat2 <- dat[dat$cog == cog,]
             dat2$test.id <- rownames(dat2)
             dat2 <- melt(dat2)
             conds <- unlist(strsplit(as.character(dat2$variable), ".R[0-9]"))
             conds <- conds[seq(1,length(conds),2)]
             dat2$cond <- conds
             dat2.sum <- ddply(dat2, c("cond", "test.id"), summarize, mean = mean(value), n = length(cond), sd = sd(value), se = sd/sqrt(n))
             dat2.sum <- dat2.sum[dat2.sum$cond == "stool.WT" | dat2.sum$cond == "stool.HhaIL10R",]
             dodge = position_dodge(width=0.9)
             plot1 <- ggplot(dat2.sum, aes(x = test.id, y = mean, fill = cond)) 
             plot2 <- plot1 + geom_bar(stat = "identity", position = dodge) 
             plot3 <- plot2 + geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.25, position = dodge) 
             plot3 + coord_flip() + scale_fill_manual(values = c("grey", "black")) + theme_bw()
             ggsave(outname)}
         ''' % cog_prefix)

    P.touch(outfile)



#########################################
#########################################
#########################################
# merge associated taxa with COG
# classifications/annotations
#########################################
#########################################
#########################################

@merge([annotateFunctions, mergeProportionTaxaInCogsRNA],
       "associate_taxa.dir/diff_genes_classified_taxa.tsv")
def mergeClassificationsAndTaxaAssociations(infiles, outfile):
    '''
    merge the annotated cogs with their per sample
    classifications
    '''
    annotations, taxa = infiles
    inf_anno = IOTools.openFile(annotations)
    inf_taxa = IOTools.openFile(taxa)
    inf_anno.readline()
    h_taxa = "\t".join(inf_taxa.readline().split("\t")[1:])

    cog2class = {}
    E.info("reading classifications")
    for line in inf_anno.readlines():
        data = line[:-1].split("\t")
        cog, classification, description = data[-2], data[-3], data[-1]
        cog2class[cog] = [classification, description]

    E.info("assigning classifications to genus assignments")
    outf = IOTools.openFile(outfile, "w")
    outf.write("cog\tclassification\tdescription\t%s" % h_taxa)
    for line in inf_taxa.readlines():
        data = line[:-1].split("\t")
        cog, dat = data[0], data[1:]
        if cog in cog2class:
            outf.write( "\t".join([cog] + cog2class[cog] + dat) + "\n")
    outf.close()

#########################################
#########################################
#########################################

@transform(annotateFunctions, suffix(".tsv"), ".pvals.pdf")
def boxplotClusterPvals(infile, outfile):
    '''
    boxplot RNA and DNA p-values for each cluster
    to make sure things all fit ok
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rna <- data.frame(fold = dat$rna, p = dat$prna, class = dat$class)''')
    R('''rna$gene <- dat$gene''')
    R('''rna$status = "rna"''')
    R('''dna <- data.frame(fold = dat$dna, p = dat$pdna, class = dat$class)''')
    R('''dna$gene <- dat$gene''')
    R('''dna$status = "dna"''')
    R('''dat2 <- data.frame(rbind(rna, dna))''')
    R('''plot1 <- ggplot(dat2, aes(x = factor(class), y = p, fill = status))''')
    R('''plot2 <- plot1 + geom_boxplot() + scale_fill_manual(values = c("grey", "white"))''')
    R('''plot3 <- plot2 + geom_hline(yintercept = 1.3, linetype = "dashed")''')
    R('''plot4 <- plot3 + geom_vline(xintercept = c(1.5,2.5,3.5,4.5), colour = "grey")''')
    R('''plot4''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("barplots.dir"))
@split(mergeClassificationsAndTaxaAssociations, "barplots.dir/*.pdf")
def barplotTaxaAssociations(infile, outfiles):
    '''
    barplot the relative contribution of genera
    to each COG
    '''
    R('''library(ggplot2)''')
    R('''library(reshape)''')
    R('''library(gtools)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''for (i in seq(5, ncol(dat),1)){
             dat[,i] <- as.numeric(dat[,i])}''')
    R('''dat$classification <- as.factor(dat$classification)''')
    R('''dat <- dat[, mixedsort(colnames(dat))]''')
    R('''dat2 <- melt(dat)''')

    R('''dat2 <- dat2[dat2$value >= 1,]''')
    R('''dat2 <- na.omit(dat2)''')
    R('''cogs <- unique(dat2$cog)''')
    R('''cogs <- cogs[!(is.na(cogs))]''')
    R('''for (cog in cogs){
                 cluster = unique(dat2$classification[dat2$cog == cog])
                 outname = paste(paste("barplots.dir", paste(cluster, cog, sep = "-"), sep = "/"), ".pdf", sep = "")
                 plot1 <- ggplot(dat2[dat2$cog == cog,], aes(x = variable, y = value, fill = taxa)) 
                 plot2 <- plot1 + geom_bar(stat = "identity") + coord_flip() 
                 ggsave(outname, height = 15, width = 15)}''')
         


#########################################
#########################################
#########################################

@transform(buildTaxaCogCountsMatrix, 
           suffix(".tsv"), 
           add_inputs(annotateFunctions),
           ".pdf")
def heatmapTaxaCogCountsMatrix(infiles, outfile):
    '''
    plot the taxa associated with each cog on
    a heatmap
    '''
    R('''library(gplots)''')
    R('''library(gtools)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infiles[0])
    R('''print (ncol(dat))''')

    R('''clusters <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])
#    R('''clusters <- clusters[clusters$class == 3,]''')
    R('''rownames(clusters) <- clusters$gene''')

    # get genes present in both - not sure why these are different
    # in the first place - need to check
    R('''genes <- intersect(rownames(clusters), colnames(dat))''')
    R('''dat <- dat[, genes]''')
    R('''dat <- dat[grep("unassigned", rownames(dat), invert = T),]''')

    R('''genera <- rownames(dat)''')
    R('''rownames(dat) <- genera''')
    R('''colnames(dat) <- genes''')
    R('''clusters <- clusters[genes,]''')

    R('''clusters <- clusters[order(clusters$class),]''')
    R('''dat <- dat[,rownames(clusters)]''')

    R('''print(ncol(dat))''')
    R('''annotation <- data.frame(cluster = as.character(clusters$class))''')
    R('''rownames(annotation) <- rownames(clusters)''')
 
    R('''colors1 <- c("orange", "purple", "red", "brown", "darkGreen")''')
    R('''names(colors1) <- c("1", "2", "3", "4", "5")''')

    R('''anno_colors <- list(cluster = colors1)''')

    R('''cols <- colorRampPalette(c("white", "darkBlue"))(150)''')

    R('''dat <- dat[rowSums(dat > 50) >= 1,]''')
    R('''dat <- dat[,colSums(dat > 50) >= 1]''')

    R('''print(ncol(dat))''')
    R('''pdf("%s", height = 10, width = 15)''' % outfile)
    R('''library(pheatmap)''')
    R('''pheatmap(dat, 
                  clustering_distance_cols = "manhattan",
                  clustering_method = "ward",
                  annotation = annotation,
                  annotation_colors = anno_colors,
                  cluster_rows = T,
                  cluster_cols = F,
                  color = cols,
                  fontsize = 8)''')
    R["dev.off"]()

#########################################
#########################################
#########################################

@transform(buildTaxaCogCountsMatrix, 
           suffix(".tsv"), 
           add_inputs(annotateFunctions),
           ".maximums.pdf")
def plotAveDistributions(infiles, outfile):
    '''
    plot the distribution of max(average) values
    per cluster
    '''
    R('''library(gplots)''')
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''Sys.setenv("DISPLAY"=":0.0")''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infiles[0])
    R('''dat <- dat[grep("unassigned", rownames(dat), invert = T), ]''')

    R('''clusters <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])
    R('''rownames(clusters) <- clusters$gene''')

    # get genes present in both - not sure why these are different
    # in the first place - need to check
    R('''genes <- intersect(rownames(clusters), colnames(dat))''')
    R('''dat <- dat[, genes]''')
    R('''genera <- rownames(dat)''')
    R('''rownames(dat) <- genera''')
    R('''colnames(dat) <- genes''')
    R('''clusters <- clusters[genes,]''')
    
    # get maximum (of the average) values
    R('''averages <- data.frame(apply(dat, 2, max))''')
    R('''colnames(averages) <- "mean"''')
    R('''averages$class <- clusters[rownames(averages),]$class''')

    # ks test
    for x, y in itertools.product([1,2,3,4,5],[1,2,3,4,5]):
        outf = "associate_taxa.dir/cluster-%i_cluster-%i.sig" % (x,y)
        R('''k <- ks.test(averages$mean[averages$class == %i],averages$mean[averages$class == %i])''' % (x,y))
        R('''k <- data.frame("D" = k[[1]], "p-value" = k[[2]])''')
        R('''write.table(k, file = "%s", sep = "\t")''' % outf)

    R('''plot1 = ggplot(averages, aes(x = mean, colour = factor(class))) + stat_ecdf()''')
    R('''plot1 + scale_colour_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("congruency.dir"))
@merge([os.path.join(PARAMS.get("dna_communitiesdir"), "csvdb"),
        mergeProportionTaxaInCogsRNA],
       "congruency.dir/congruency.tsv")
def buildCongruencyTable(infiles, outfile):
    '''
    build a table of whether or not the taxa for each gene
    are congruent i.e. abundance in the same direction
    or not - number up and number down
    '''
    E.info("reading taxa associations")
    cog2taxa = collections.defaultdict(set)

    inf = IOTools.openFile(infiles[1])
    header = inf.readline()
    for line in inf.readlines():
        data = line[:-1].split("\t")
        if len(data) > 18:
            cog, taxa = data[0], " ".join(data[1:3])
            values = map(float,data[3:])
        else:
            cog, taxa = data[0], data[1]
            values = map(float,data[2:])
    
        if len([x for x in values if x > 0]) >= 8:
            cog2taxa[cog].add(taxa)
                 
    E.info("assessing congruency")
    up = False
    dbh = sqlite3.connect(infiles[0])
    cc = dbh.cursor()
    result = collections.defaultdict(list)
    for data in cc.execute("""SELECT taxa, logFC 
                              FROM genus_diamond_aggregated_counts_diff
                              WHERE group1 == "HhaIL10R" 
                              AND group2 == "WT" 
                              """).fetchall():
        genus, logfc = data[0], data[1]
        if logfc > 0:
            status = 1
        elif logfc < 0:
            status = 2
        for cog, genera in cog2taxa.iteritems():
            if genus in genera:
                result[cog].append(logfc)

    outf = IOTools.openFile(outfile, "w")
    outf2 = open("testfcs.tsv", "w")
    outf.write("gene\tnup\tndown\tnetfc\tngenera\n")
    for cog, fcs in result.iteritems():
        outf2.write("\t".join([cog] + map(str,fcs)) + "\n")
        net = sum(fcs)
        nup = len([x for x in fcs if x > 0])
        ndown = len([x for x in fcs if x < 0])
        ngenera = len(fcs)
        if nup == 0:
            pup = "inf"
        elif ndown == 0:
            pup = "-inf"
        else:
            pup = float(nup)/float(ndown)

        outf.write("%(cog)s\t%(nup)s\t%(ndown)s\t%(net)s\t%(ngenera)s\n" % locals())
    outf.close()
    outf2.close()

#########################################
#########################################
#########################################

@merge([buildCongruencyTable, annotateFunctions],
       "congruency.dir/congruency_class.tsv")
def mergeCongruencyAndClasses(infiles, outfile):
    '''
    merge the congruency and cluster
    classifications
    '''
    df_cong = pandas.read_csv(infiles[0], sep = "\t")
    df_class = pandas.read_csv(infiles[1], sep = "\t")
    result = pandas.merge(df_cong, df_class, on = "gene")
    result.to_csv(outfile, sep = "\t", na_rep = 0, index = False)
    
#########################################
#########################################
#########################################

@transform(mergeCongruencyAndClasses, suffix(".tsv"), ".pdf")
def plotCongruencyByClass(infile, outfile):
    '''
    plot congruency by cluster class
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''plot1 <- ggplot(dat, aes(x = factor(class), y = log2(pup), colour = factor(class), alpha = 0.5))''')
    R('''plot2 <- plot1 + geom_boxplot()''')
    R('''plot2 + scale_colour_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@transform(mergeCongruencyAndClasses, suffix(".tsv"), ".ngenera.pdf")
def plotNumberOfGeneraByClass(infile, outfile):
    '''
    barplot the number of genera found for each
    class
    '''
    R('''library(ggplot2)''')
    R('''library(plyr)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''dat2 <- ddply(dat, c("class"),
                       summarise,
                       mean = mean(ngenera),
                       n = length(ngenera),
                       sd = sd(ngenera),
                       se = sd/sqrt(n))''')
    R('''plot1 <- ggplot(dat2, aes(x = factor(class), y = mean, fill = factor(class)))''')
    R('''plot2 <- plot1 + geom_bar(position = "dodge", stat = "identity")''')
    R('''plot3 <- plot2 + geom_errorbar(aes(ymax = mean + se, ymin = mean - se), position = "dodge", width = 0.25)''')
    R('''plot3 + scale_fill_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@transform(mergeCongruencyAndClasses, suffix(".tsv"), ".ngenera.stats.tsv")
def runTukeysNumberOfGeneraByClass(infile, outfile):
    '''
    run tukeys test on number of genera
    '''
    R('''library(ggplot2)''')
    R('''library(plyr)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''aov1 <- aov(dat$ngenera ~ factor(dat$class))''')
    R('''thsd <- TukeyHSD(aov1)''')
    R('''write.table(thsd[[1]], file = "%s", sep = "\t")''' % outfile)

#########################################
#########################################
#########################################

@transform(mergeCongruencyAndClasses, suffix(".tsv"), ".load")
def loadCongruencyAndClasses(infile, outfile):
    '''
    load the annotations file into database
    '''
    P.load(infile, outfile)

#########################################
#########################################
#########################################

# @transform(accumulateCounts, suffix(".accumulated.counts.tsv"), ".cog.pdf")
# def countDiffCogsPerGenus(infile, outfile):
#     '''
#     plot the number of COGs that each genus is associated
#     with
#     '''
#     R('''library(ggplot2)''')
#     R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infile)

#     # get the number of COGS hit by a genus
#     R('''d <- data.frame(rowSums(dat > 0))''')
#     R('''colnames(d) <- "count"''')
#     R('''d$taxa <- rownames(d)''')

#     # get number of COGs
#     R('''ncogs <- length(colnames(dat))''')
#     R('''d$count <- (d$count/ncogs)*100''')
    
#     # only plot those that have a max > 10%
#     R('''ggplot(d[d$count > 10,], aes(x = taxa, y = count)) + geom_bar(stat = "identity") + coord_flip()''')
#     R('''ggsave("%s", height = 15)''' % outfile)

#########################################
#########################################
#########################################

@transform(mergeProportionTaxaInCogsRNA, suffix(".tsv.gz"), ".load")
def loadAssociatedTaxa(infile, outfile):
    '''
    load taxa associated with diff genes
    '''
    P.load(infile, outfile)



#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))





