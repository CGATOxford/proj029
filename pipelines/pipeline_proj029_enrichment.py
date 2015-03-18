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
import CGAT.Pipeline as P

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
@transform("classes_ratio.tsv", regex("(\S+).tsv"), r"pathways.dir/\1.foreground.tsv.gz")
def buildForegroundSets(infile, outfile):
    '''
    build foreground set of COGs
    '''
    inf = open(infile)
    inf.readline()
    cogs = set()

    outf = IOTools.openFile(outfile, "w")
    outf.write("gene\tC1\tC2\tC3\tC4\tC5\tup\tdown\n")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        cog, cluster = data[0], data[5]
        if cog in cogs: continue
        cogs.add(cog)
        if cluster == "1":
            outf.write("%s\t1\t0\t0\t0\t0\t0\t1\n" % cog)
        elif cluster == "2":
            outf.write("%s\t0\t1\t0\t0\t0\t1\t0\n" % cog)
        elif cluster == "3":
            outf.write("%s\t0\t0\t1\t0\t0\t1\t0\n" % cog)
        elif cluster == "4":
            outf.write("%s\t0\t0\t0\t1\t0\t0\t1\n" % cog)
        elif cluster == "5":
            outf.write("%s\t0\t0\t0\t0\t1\t1\t0\n" % cog)
    outf.close()

#########################################
#########################################
#########################################

@split([buildForegroundSets,
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

@merge([buildForegroundSets, PARAMS.get("pathways_geneset")],
       "pathways.dir/cogs_pathways.tsv")
def buildDiffCogsAndPathways(infiles, outfile):
    '''
    merge diff COGs and pathways
    '''
    R('''cogs <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infiles[0])
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

    R('''dat$cluster <- ifelse(dat$C1 == 1, "C1", NA)''')
    R('''dat$cluster <- ifelse(dat$C2 == 1, "C2", dat$cluster)''')
    R('''dat$cluster <- ifelse(dat$C3 == 1, "C3", dat$cluster)''')
    R('''dat$cluster <- ifelse(dat$C4 == 1, "C4", dat$cluster)''')
    R('''dat$cluster <- ifelse(dat$C5 == 1, "C5", dat$cluster)''')
    
    # remove NAs and unknown functions
    R('''dat <- na.omit(dat)''')
    R('''dat <- dat[grep("Function unknown", dat$V4, invert = T),]''')
    R('''dat <- dat[grep("General function", dat$V4, invert = T),]''')

    # get counts per cluster
    R('''counts <- ddply(dat, c("cluster"), summarise, n = length(cluster))''')
    R('''rownames(counts) <- counts$cluster''')

    # add counts to main data
    R('''dat$count <- counts[dat$cluster,]$n''')

    # summarise per functional category
    R('''dat2 <- ddply(dat, c("V4", "cluster"), summarise, prop = (length(cluster)/mean(count))*100)''')
    R('''plot1 <- ggplot(dat2, aes(x = V4, y = prop, fill = cluster)) + geom_bar(stat = "identity", position = "dodge") + facet_grid(~cluster) + coord_flip()''')
    R('''plot1 + scale_fill_manual(values = c("orange", "purple", "red", "brown", "darkGreen"))''')
    R('''ggsave("%s", width = 15)''' % outfile)

#########################################
#########################################
#########################################


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
