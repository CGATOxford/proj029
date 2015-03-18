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

@follows(mkdir("metrics.dir"))
@merge([os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
        os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.diff.tsv")],
       "metrics.dir/expression_cv.png")
def plotAveExprVsCv(infiles, outfile):
    '''
    plot the average expression of COGS vs
    coefficient of variation
    '''
    R('''library(ggplot2)''')
    R('''diff <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])
    R('''mat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''rownames(mat) <- mat$taxa''')
    R('''mat <- mat[, 1:ncol(mat)-1]''')
    R('''mat <- mat[diff$taxa,]''')
    R('''cvs <- apply(mat, 1, function(x) sd(x)/mean(x))''')
    R('''diff$cv <- cvs''')
    R('''ggplot(diff, aes(x = AveExpr, y = cv, size = -log10(P.Value))) + geom_point()''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("metrics.dir"))
@merge([os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
        os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.diff.tsv")],
       "metrics.dir/MA.png")
def plotMA(infiles, outfile):
    '''
    plot the average expression of COGS vs
    coefficient of variation
    '''
    R('''library(ggplot2)''')
    R('''diff <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])
    R('''diff <- diff[diff$group1 == "HhaIL10R" & diff$group2 == "WT",]''')
    R('''mat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''rownames(mat) <- mat$taxa''')
    R('''mat <- mat[, 1:ncol(mat)-1]''')
    R('''mat <- mat[diff$taxa,]''')
    R('''ggplot(diff, aes(x = AveExpr, y = logFC, size = -log10(P.Value))) + geom_point() + geom_hline(yintercept = c(-1,1), linetype = "dashed")''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################

@follows(mkdir("metrics.dir"))
@merge([os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.norm.matrix"),
       os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.diff.tsv")],
       "metrics.dir/heatmap.png")
def heatmapCandidates(infiles, outfile):
    '''
    heatmap those COGS that are:
    expressed > 15
    '''
    R('''library(gplots)''')
    R('''library(gtools)''')
    R('''diff <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[1])
    R('''diff <- diff[diff$group1 == "HhaIL10R" & diff$group2 == "WT",]''')
    R('''cogs <- diff$taxa[diff$AveExpr > 15]''')

    R('''mat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''rownames(mat) <- mat$taxa''')
    R('''mat <- mat[, 1:ncol(mat)-1]''')
    R('''mat <- mat[, mixedsort(colnames(mat))]''')
    R('''mat <- mat[cogs,]''')

    R('''mat <- mat[,9:16]''')
    R('''mat <- mat[order(rowMeans(mat), decreasing = T),]''')

    R('''cols <- colorRampPalette(c("darkGreen", "black", "purple"))(75)''')
    R('''mat.s <- data.frame(t(apply(mat, 1, scale)))''')
    R('''png("%s")''' % outfile)
    R('''heatmap.2(as.matrix(mat), trace = "none", col = cols, Rowv = F, Colv = F, margins = c(15,15))''')
    R["dev.off"]()

#########################################
#########################################
#########################################

@follows(mkdir("diff.dir"))
@transform(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/gene_counts.diff.tsv"),
           regex("(\S+)/(\S+).diff.tsv"),
           r"diff.dir/\2.candidate.list.tsv")
def buildCandidateGeneList(infile, outfile):
    '''
    build candidate list
    '''
    outf = open(outfile, "w")
    inf = open(infile)
    inf.readline()
    for line in inf:
        data = line[:-1].split("\t")
        if data[6] == '"HhaIL10R"' and data[7] == '"WT"':
            if float(data[1]) < 15: continue
            outf.write("%s\n" % data[8].replace('"',""))
    outf.close()

#########################################
#########################################
#########################################
# associate genes with taxa and see
# taxa specific gene expression patterns
#########################################
#########################################
#########################################

@follows(mkdir("associate_taxa.dir"), buildCandidateGeneList)
@transform(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "genes.dir/*.genes.tsv.gz")),
            regex("(\S+)/(\S+).genes.tsv.gz"),
            add_inputs(glob.glob(os.path.join(PARAMS.get("rna_communitiesdir"), "diamond.dir/*.diamond.lca.gz"))),
            r"associate_taxa.dir/\2.ctaxa.tsv.gz")
def associateCandidateGenesWithTaxaCounts(infiles, outfile):
    '''
    reads that map to differentially expressed genes are cross-referenced
    with their genera assignment - provides per taxa counts
    '''
    job_options="-l mem_free=25G"
    m = PARAMS.get("genes_map")
    alignment_genes = infiles[0]
    track = P.snip(os.path.basename(alignment_genes), ".diamond.genes.tsv.gz")
    alignment_taxa = [x for x in infiles[1] if os.path.basename(x).startswith(track)][0]

    statement = '''python %(projscripts)s/diff2genera.py 
                   -m %(m)s
                   -d diff.dir/gene_counts.candidate.list.tsv
                   --alignment-taxa=%(alignment_taxa)s
                   --alignment-genes=%(alignment_genes)s
                   --counts
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()


#########################################
#########################################
#########################################

@merge(associateCandidateGenesWithTaxaCounts, "associate_taxa.dir/associated_taxa_counts.tsv.gz")
def mergeAssociatedTaxaCounts(infiles, outfile):
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
                    > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@follows(mkdir("taxa_cogs_diff.dir"))
@transform(mergeAssociatedTaxaCounts,
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

@jobs_limit(1,"R")
@transform(runMetagenomeSeqPerCOGAndTaxa, 
           suffix(".diff.tsv"), 
           ".pdf")
def plotCandidates(infile, outfile):
    '''
    barplot differences for different taxa for candidate
    COGS
    just in WT vs. HhaIL10R
    '''

    matrix = P.snip(infile, ".diff.tsv") + ".norm.matrix"

    R('''library(ggplot2)''')
    R('''library(gtools)''')
    R('''library(plyr)''')
    R('''library(reshape)''')

    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % matrix)
    R('''rownames(dat) <- dat$taxa''')

    # hack to get around the "-" present in a rowname - replace
    R('''rownames(dat)[grep("SPO1", rownames(dat))] <- gsub("SPO1-", "SPO1_", rownames(dat)[grep("SPO1", rownames(dat))])''')
    R('''rownames(dat)[grep("Lambda", rownames(dat))] <- gsub("Lambda-", "Lambda_", rownames(dat)[grep("Lambda", rownames(dat))])''')
    R('''dat <- dat[,1:ncol(dat)-1]''')

    # filter for expression level
    R('''dat <- dat[rowMeans(dat) > 5,]''')

    R('''cogs <- unlist(strsplit(rownames(dat), "-", perl = TRUE))''')
    R('''cogs <- cogs[seq(1,length(cogs),2)]''')


    cog_prefix = "rna"

    R('''dat$cog <- cogs''')

    R('''cogs <- unique(cogs)''')
    R('''for (cog in cogs){
             outname <- paste(paste("taxa_cogs_diff.dir", paste("%s", cog, sep = "_"), sep = "/"), "pdf", sep = ".") 
             dat2 <- dat[dat$cog == cog,]
             dat2$test.id <- rownames(dat2)
             dat2 <- melt(dat2)
             conds <- unlist(strsplit(as.character(dat2$variable), ".R[0-9]"))
             conds <- conds[seq(1,length(conds),2)]
             dat2$cond <- conds
             dat2.sum <- ddply(dat2, c("cond", "test.id"), summarize, mean = mean(value), n = length(cond), sd = sd(value), se = sd/sqrt(n))
             #dat2.sum <- dat2.sum[dat2.sum$cond == "stool.WT" | dat2.sum$cond == "stool.HhaIL10R",]
             dodge = position_dodge(width=0.9)
             plot1 <- ggplot(dat2.sum, aes(x = test.id, y = mean, fill = cond)) 
             plot2 <- plot1 + geom_bar(stat = "identity", position = dodge) 
             plot3 <- plot2 + geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.25, position = dodge) 
             plot3 + coord_flip() + scale_fill_manual(values = c("purple", "darkGreen", "red", "blue")) + theme_bw()
             ggsave(outname)}
         ''' % cog_prefix)

    P.touch(outfile)

#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))




