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

import os
import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R
import pandas
import sys
import CGATPipelines.Pipeline as P
import collections

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


@merge(["common_genes.tsv", PARAMS.get("pathways_geneset")],
       "pathways.dir/background_cogs_pathways.tsv")
def buildBackgroundCogsAndPathways(infiles, outfile):
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


@transform(runPathwaysAnalysis, suffix(".overall"), ".bar.pdf")
def plotPathways(infile, outfile):
    '''
    plot pathways associated with clusters
    '''
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % infile)
    R('''dat <- dat[order(dat$ratio),]''')
    R('''plot1 <- ggplot(dat, aes(x = factor(description, levels=description), y=ratio, stat="identity"))''')
    R('''plot1 + geom_bar(stat="identity") + coord_flip()''')
    R('''ggsave("%s")''' % outfile)

#########################################
#########################################
#########################################


@follows(mkdir("heatmaps.dir"))
@split(buildDiffCogsAndPathways, "heatmaps.dir/*.tsv")
def splitPathways(infile, outfiles):
    '''
    map cogs to pathways in separate files
    '''
    inf = IOTools.openFile(infile)
    inf.readline()
    pathway2nogs = collections.defaultdict(set)
    for line in inf.readlines():
        data = line[:-1].split("\t")
        nog, pathway = data[0], data[3]
        pathway = pathway.replace(" ", "_").replace("/", "_")
        pathway2nogs[pathway].add(nog)
    for pathway, nogs in pathway2nogs.iteritems():
        outname = os.path.join("heatmaps.dir", pathway + ".tsv")
        outf = IOTools.openFile(outname, "w")
        outf.write("NOG\tpathway\n")
        for nog in nogs:
            outf.write("%s\t%s\n"% (nog, pathway))
    outf.close()

#########################################
#########################################
#########################################


@follows(mkdir("annotations.dir"))
@transform(splitPathways,
           regex("(\S+)/(\S+).tsv"),
           add_inputs(PARAMS.get("annotations_eggnog")),
           r"annotations.dir/\2.annotated.tsv")
def annotateNogs(infiles, outfile):
    '''
    annotate the NOGs with their descriptions
    '''
    pathways, annotations = infiles
    anno = {}
    # read annotations
    for line in open(annotations).readlines():
        data = line[:-1].split("\t")
        nog, description = data
        anno[nog] = description

    # write out annotations
    p = IOTools.openFile(pathways)
    p.readline()
    outf = IOTools.openFile(outfile, "w")
    outf.write("NOG\tpathway\tdescription\n")
    for line in p.readlines():
        data = line[:-1].split("\t")
        nog, pathway = data
        try:
            outf.write("%s\t%s\t%s\n" % (nog, pathway, anno[nog]))
        except KeyError:
            outf.write("%s\t%s\t%s\n" % (nog, pathway, "NA"))
    outf.close()

#########################################
#########################################
#########################################


@jobs_limit(1,"R")
@transform(splitPathways,
           suffix(".tsv"),
           add_inputs(PARAMS.get("matrix_file")),
           ".heatmap.pdf")
def heatmapNogs(infiles, outfile):
    '''
    heatmap nogs per functional category
    '''
    # check files is compatible with heatmaps
    # i.e. >=2 NOGs
    nogs, matrix = infiles
    inf = IOTools.openFile(nogs)

    # header
    inf.readline()

    # check lines
    if len(inf.readlines()) == 1:
        P.touch(outfile)
    else:
        R('''library(gtools)''')
        R('''library(gplots)''')

        # read in matrix
        R('''mat <- read.csv("%s", header=T,stringsAsFactors=F, sep="\t")''' % matrix)
        R('''rownames(mat) <- mat$taxa''')
        R('''mat <- mat[,1:ncol(mat)-1]''')

        # read in NOGs
        R('''nogs <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % nogs)
        
        # subset matrix
        R('''mat <- mat[nogs$NOG,]''')
        
        # scale
        R('''mat.s <- data.frame(t(apply(mat, 1, scale)))''')
        R('''colnames(mat.s) <- colnames(mat)''')
        R('''mat.s <- mat.s[, mixedsort(colnames(mat.s))]''')
        R('''mat.s <- mat.s[order(rownames(mat.s)),]''')

        # heatmap
        R('''pdf("%s")''' % outfile)
        R('''cols <- colorRampPalette(c("blue", "white", "red"))(75)''')
        R('''heatmap.2(as.matrix(mat.s),
                       trace="none",
                       Colv=F,
                       Rowv=F,
                       col=cols,
                       margins=c(15,15))''')
        R["dev.off"]()


#########################################
#########################################
#########################################


@jobs_limit(1,"R")
@transform(splitPathways,
           suffix(".tsv"),
           add_inputs(PARAMS.get("matrix_taxa")),
           ".taxa.heatmap.pdf")
def heatmapTaxaAssociatedWithNogs(infiles, outfile):
    '''
    heatmap nogs per functional category
    '''
    # check files is compatible with heatmaps
    # i.e. >=2 NOGs
    nogs, matrix = infiles
    inf = IOTools.openFile(nogs)

    # header
    inf.readline()

    # check lines
    if len(inf.readlines()) == 1:
        P.touch(outfile)
    else:
        R('''library(gtools)''')
        R('''library(gplots)''')
        R('''library(pheatmap)''')

        # read in matrix
        R('''mat <- t(read.csv("%s", header=T,stringsAsFactors=F, sep="\t", row.names=1))''' % matrix)

        # read in NOGs
        R('''nogs <- read.csv("%s", header=T, stringsAsFactors=F, sep="\t")''' % nogs)
        
        # subset matrix
        R('''mat <- mat[nogs$NOG,]''')
        R('''mat <- mat[, colSums(mat > 5) >= 1]''')
        R('''mat <- mat[order(rownames(mat)),]''')

        # heatmap
        R('''pdf("%s")''' % outfile)
        R('''cols <- colorRampPalette(c("white", "blue"))(75)''')
        R('''pheatmap(as.matrix(mat),
                       color=cols,
                       cluster_cols=F,
                       cluster_rows=F)''')

        R["dev.off"]()

@follows(heatmapNogs,
         heatmapTaxaAssociatedWithNogs)
def heatmaps():
    pass

@follows(runPathwaysAnalysis,annotateNogs,heatmaps)
def full():
    pass


#########################################
#########################################
#########################################


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
