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
    R('''p3 <- p2 + scale_colour_manual(values = c("slateGrey", "green", "red", "blue"))''')
#    R('''p3 + xlim(c(-%i, %i)) + ylim(c(-%i, %i))''' % (xlim, xlim, ylim, ylim))
    R('''ggsave("%s")''' % outname_plot)

    get the loadings
    R('''loads <- data.frame(pc$rotation)''')
    R('''loads$taxa <- rownames(loads)''')

    # write out data
    R('''write.table(loads, file = "%s", sep = "\t", row.names = F, quote = F)''' % outfile.replace("/", "/%s_" % suffix))

    P.touch(outfile)

#########################################
#########################################
#########################################

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
