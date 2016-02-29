"""
=======================================
Compare original and new RNA-seq data
at Day14 vs. Day0
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
import CGATPipelines.PipelineMetagenomeCommunities as PipelineMetagenomeCommunities
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

@follows(mkdir("pca.dir"))
@jobs_limit(1, "R")
@transform([os.path.join(PARAMS.get("communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
            os.path.join(PARAMS.get("communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix")],
           regex("(\S+)/(\S+).matrix"),
           r"pca.dir/\2.loadings.tsv")
def buildPCALoadings(infile, outfile):
    '''
    run PCA and heatmap the loadings
    '''
    outname_plot = P.snip(outfile, ".loadings.tsv")  + ".pca.pdf"
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)

    # just get day14 and day0
    R('''remove <- c("day3", "day6", "day28")''')
    R('''for (day in remove){; dat <- dat[, grep(day, colnames(dat), invert=T)]}''')

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
    R('''p3 <- p2 + scale_colour_manual(values = c("slateGrey", "red"))''')

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

@transform([os.path.join(PARAMS.get("communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
            os.path.join(PARAMS.get("communitiesdir"), "counts.dir/genus.diamond.aggregated.counts.norm.matrix")],
           regex("(\S+)/(\S+).matrix"),
           r"pca.dir/\2.sig")
def testDistSignificance(infile, outfile):
    '''
    test whether the colitic samples
    cluster significantly
    '''
    PipelineMetagenomeCommunities.testDistSignificance(infile,
                                                       outfile)


#########################################
#########################################
#########################################

@follows(mkdir("correlation.dir"))
@merge(["original_gene_counts.diff.tsv",
        "replication_gene_counts.diff.tsv"],
       "correlation.dir/gene_abundance_scatter.png")
def scatterplotAbundanceEstimates(infiles, outfile):
    '''
    scatterplot abundance estimates for NOGs
    '''
    R('''dat.orig <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % infiles[0])
    R('''dat.orig <- dat.orig[dat.orig$group2 == "WT" & dat.orig$group1 == "HhaIL10R",]''')

    R('''dat.rep <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % infiles[1])
    R('''dat.rep <- dat.rep[dat.rep$group1 == "day0" & dat.rep$group2 == "day14",]''')
    
    R('''rownames(dat.orig) <- dat.orig$taxa''')
    R('''dat.orig <- dat.orig[dat.rep$taxa,]''')

    R('''png("%s")''' % outfile)
    R('''plot(dat.orig$AveExpr, dat.rep$AveExpr, pch=16, col="slateGrey")''')
    R('''abline(0,1)''')
    R('''text(x=3, y=15, labels=c(paste("r =", round(cor(dat.orig$AveExpr, dat.rep$AveExpr),2), sep=" ")))''')
    R["dev.off"]()

#########################################
#########################################
#########################################

@follows(mkdir("correlation.dir"))
@merge(["original_gene_counts.diff.tsv",
        "replication_gene_counts.diff.tsv"],
       "correlation.dir/gene_fold_changes_scatter.png")
def scatterplotFoldChanges(infiles, outfile):
    '''
    scatterplot abundance estimates for NOGs
    '''
    R('''dat.orig <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % infiles[0])
    R('''dat.orig <- dat.orig[dat.orig$group2 == "WT" & dat.orig$group1 == "HhaIL10R",]''')

    R('''dat.rep <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % infiles[1])
    R('''dat.rep <- dat.rep[dat.rep$group1 == "day0" & dat.rep$group2 == "day14",]''')
    
    R('''rownames(dat.orig) <- dat.orig$taxa''')
    R('''dat.orig <- dat.orig[dat.rep$taxa,]''')

    R('''png("%s")''' % outfile)
    R('''plot(dat.orig$logFC, -1*dat.rep$logFC, pch=16)''')
    R('''text(x=-4, y=5, labels=c(paste("r =", round(cor(dat.orig$logFC, -1*dat.rep$logFC),2), sep=" ")))''')
    R["dev.off"]()

#########################################
#########################################
#########################################

@follows(mkdir("correlation.dir"))
@merge(["original_gene_counts.diff.tsv",
        "replication_gene_counts.diff.tsv"],
       "correlation.dir/gene_diff_overlap.tsv")
def buildGeneDifferentialExpressionOverlap(infiles, outfile):
    '''
    scatterplot abundance estimates for NOGs
    '''
    R('''dat.orig <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % infiles[0])
    R('''dat.orig <- dat.orig[dat.orig$group2 == "WT" & dat.orig$group1 == "HhaIL10R",]''')

    R('''dat.rep <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % infiles[1])
    R('''dat.rep <- dat.rep[dat.rep$group1 == "day0" & dat.rep$group2 == "day14",]''')
    
    R('''rownames(dat.orig) <- dat.orig$taxa''')
    R('''dat.orig <- dat.orig[dat.rep$taxa,]''')

    R('''diff.orig <- dat.orig$taxa[dat.orig$adj.P.Val < 0.05]''')
    R('''diff.rep <- dat.rep$taxa[dat.rep$adj.P.Val < 0.05]''')
    
    R('''overlap <- intersect(diff.orig, diff.rep)''')
    R('''write.table(overlap, file="%s", sep="\t")''' % outfile)

    R('''norig <- length(diff.orig)''')
    R('''nrep <- length(diff.rep)''')
    R('''noverlap <- length(overlap)''')

    # significance testing
    R('''x <- length(intersect(dat.orig$taxa, dat.rep$taxa))''')
    R('''m <- nrep''')
    R('''n <- x - nrep''')
    R('''k <- norig''')
    R('''print(1-phyper(x,m,n,k))''')
    R('''write.table(data.frame(c(norig, nrep,noverlap)), file="correlation.dir/noverlap.tsv")''')

#########################################
#########################################
#########################################

@follows(mkdir("wolinella_weisella.dir"))
@transform("original_genus.diamond.aggregated.counts.norm.matrix",
           regex("(\S+).norm.matrix"),
           r"wolinella_weisella.dir/\1.wolinella.pdf")
def plotOriginalWolinella(infile, outfile):
    '''
    plot the abundance of Weisella and Wolinella;
    ones that replicated
    '''
    R('''library(reshape)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", 
                         header=T,
                         stringsAsFactors=F,
                         sep="\t")''' % infile)
    R('''dat <- melt(dat)''')
    R('''conds <- unlist(strsplit(as.character(dat$variable), ".R[0-9]"))''')
    R('''conds <- conds[seq(1,length(conds),2)]''')
    R('''dat$cond <- conds''')
    R('''dat <- dat[dat$taxa == "Wolinella",]''')
    R('''plot1 <- ggplot(dat, aes(x=factor(cond, levels=c("stool.WT","stool.aIL10R", "stool.Hh", "stool.HhaIL10R")), 
                         y=value, group=cond, colour=cond))''')
    R('''plot2 <- plot1 + geom_boxplot() + geom_jitter(size=3)''')
    R('''plot2 + scale_colour_manual(values=c("blue", "darkGreen", "red", "grey")) + ylim(c(0,3))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("wolinella_weisella.dir"))
@transform("replication_genus.diamond.aggregated.counts.norm.matrix",
           regex("(\S+).norm.matrix"),
           r"wolinella_weisella.dir/\1.wolinella.pdf")
def plotReplicationWolinella(infile, outfile):
    '''
    plot the abundance of Weisella and Wolinella;
    ones that replicated
    '''
    R('''library(reshape)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", 
                         header=T,
                         stringsAsFactors=F,
                         sep="\t")''' % infile)
    # just get day14 and day0
    R('''remove <- c("day3", "day6", "day28")''')
    R('''for (day in remove){; dat <- dat[, grep(day, colnames(dat), invert=T)]}''')

    R('''dat <- melt(dat)''')
    R('''conds <- unlist(strsplit(as.character(dat$variable), ".R[0-9]"))''')
    R('''conds <- conds[seq(1,length(conds),2)]''')
    R('''dat$cond <- conds''')
    R('''dat <- dat[dat$taxa == "Wolinella",]''')
    R('''plot1 <- ggplot(dat, aes(x=factor(cond, levels=c("stool.day0","stool.day14")), 
                         y=value, group=cond, colour=cond))''')
    R('''plot2 <- plot1 + geom_boxplot() + geom_jitter(size=3)''')
    R('''plot2 + scale_colour_manual(values=c("grey", "red")) + ylim(c(0,3))''')
    R('''ggsave("%s")''' % outfile)


#########################################
#########################################
#########################################

@follows(mkdir("wolinella_weisella.dir"))
@transform("original_genus.diamond.aggregated.counts.norm.matrix",
           regex("(\S+).norm.matrix"),
           r"wolinella_weisella.dir/\1.weissella.pdf")
def plotOriginalWeissella(infile, outfile):
    '''
    plot the abundance of Weisella and Wolinella;
    ones that replicated
    '''
    R('''library(reshape)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", 
                         header=T,
                         stringsAsFactors=F,
                         sep="\t")''' % infile)
    R('''dat <- melt(dat)''')
    R('''conds <- unlist(strsplit(as.character(dat$variable), ".R[0-9]"))''')
    R('''conds <- conds[seq(1,length(conds),2)]''')
    R('''dat$cond <- conds''')
    R('''dat <- dat[dat$taxa == "Weissella",]''')
    R('''plot1 <- ggplot(dat, aes(x=factor(cond, levels=c("stool.WT","stool.aIL10R", "stool.Hh", "stool.HhaIL10R")), 
                         y=value, group=cond, colour=cond))''')
    R('''plot2 <- plot1 + geom_boxplot() + geom_jitter(size=3)''')
    R('''plot2 + scale_colour_manual(values=c("blue", "darkGreen", "red", "grey")) + ylim(c(0,4))''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("wolinella_weisella.dir"))
@transform("replication_genus.diamond.aggregated.counts.norm.matrix",
           regex("(\S+).norm.matrix"),
           r"wolinella_weisella.dir/\1.weissella.pdf")
def plotReplicationWeissella(infile, outfile):
    '''
    plot the abundance of Weisella and Wolinella;
    ones that replicated
    '''
    R('''library(reshape)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", 
                         header=T,
                         stringsAsFactors=F,
                         sep="\t")''' % infile)
    # just get day14 and day0
    R('''remove <- c("day3", "day6", "day28")''')
    R('''for (day in remove){; dat <- dat[, grep(day, colnames(dat), invert=T)]}''')

    R('''dat <- melt(dat)''')
    R('''conds <- unlist(strsplit(as.character(dat$variable), ".R[0-9]"))''')
    R('''conds <- conds[seq(1,length(conds),2)]''')
    R('''dat$cond <- conds''')
    R('''dat <- dat[dat$taxa == "Weissella",]''')
    R('''plot1 <- ggplot(dat, aes(x=factor(cond, levels=c("stool.day0","stool.day14")), 
                         y=value, group=cond, colour=cond))''')
    R('''plot2 <- plot1 + geom_boxplot() + geom_jitter(size=3)''')
    R('''plot2 + scale_colour_manual(values=c("grey", "red")) + ylim(c(0,4))''')
    R('''ggsave("%s")''' % outfile)

@follows(plotOriginalWeissella,
         plotOriginalWolinella,
         plotReplicationWeissella,
         plotReplicationWolinella)
def plotWolinellaWeissella():
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
