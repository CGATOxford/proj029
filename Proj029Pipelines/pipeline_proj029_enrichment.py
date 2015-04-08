"""
=======================================
Perform functional enrichment testing
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
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R
import pandas
import sys
import CGATPipelines.Pipeline as P

#######################
# parameters
#######################
P.getParameters(
    ["pipeline.ini"])
PARAMS = P.PARAMS

#########################################
#########################################
#########################################

@follows(mkdir("pathways.dir"))
@transform("ratio_genes.annotated.outsidepi.tsv", regex("(\S+).tsv"), r"pathways.dir/\1.foreground.tsv.gz")
def buildForegroundSet(infile, outfile):
    '''
    build foreground set of COGs
    '''
    status = PARAMS.get("group_status")
    statement = '''cat %(infile)s | grep %(status)s | cut -f1 | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################

@split([buildForegroundSet,
        "common_genes.tsv",
        PARAMS.get("pathways_geneset")],
       "pathways.dir/*.overall")
def runPathwaysAnalysis(infiles, outfiles):
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
                   --output-filename-pattern="pathways.dir/%%(set)s.%%(go)s.%%(section)s" \                                                                                                                                                  
                   > pathways.dir/pathways.log  \                                                                                                                                                                                            
                ; rm -rf %(temp)s
                ''' 
    P.run()

#########################################
#########################################
#########################################

@merge([buildForegroundSet, PARAMS.get("pathways_geneset")],
       "pathways.dir/cogs_pathways.tsv")
def buildDiffCogsAndPathways(infiles, outfile):
    '''
    merge diff COGs and pathways
    '''
    R('''cogs <- read.csv("%s", header = F, stringsAsFactors = F, sep = "\t")''' % infiles[0])
    R('''colnames(cogs) <- "gene"''')
    R('''pathways <- read.csv("%s", header = F, stringsAsFactors = F, sep = "\t")''' % infiles[1])
    R('''dat <- merge(cogs, pathways, by.x = "gene", by.y = "V2", all.x = T, all.y = F)''')
    R('''write.table(dat, file = "%s", sep = "\t", row.names = F, quote = F)''' % outfile)

#########################################
#########################################
#########################################

@transform(buildDiffCogsAndPathways, suffix(".tsv"), ".pdf")
def plotPathways(infile, outfile):
    '''
    plot pathways associated with clusters
    '''
    R('''library(plyr)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''dat <- dat[,1:ncol(dat)-1]''')

    # put NAs into "Function unknown"
    # category
    R('''dat$V4[is.na(dat$V4)] <- "Function unknown"''')

    # summarise per functional category
    R('''dat2 <- ddply(dat, c("V4"), summarise, n = length(V4))''')
    R('''dat2 <- dat2[order(dat2$n),]''')
    R('''plot1 <- ggplot(dat2, aes(x = factor(V4, levels = V4), y = n)) + geom_bar(stat = "identity", position = "dodge") + coord_flip()''')
    R('''ggsave("%s", width = 15)''' % outfile)

#########################################
#########################################
#########################################


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
