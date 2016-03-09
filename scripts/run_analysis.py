'''
run_analysis.py
=================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will run through the analysis steps as described in the documentation. It will not 
run the computationally intensive steps but will link to files that have already been created.

Usage
-----



Example::

   python run_analysis.py 

Type::

   python run_analysis.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import glob

# import CGAT modules
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineMetagenomeCommunities as PipelineMetagenomeCommunities
import Proj029Pipelines.PipelineMetaomics as PipelineMetaomics


#######################
# function definitions
#######################


def linkData(datadir):
    '''
    link data from data directory as specified
    on the commandline
    '''
    # link goi.tsv to cwd
    
    if not os.path.exists("goi.tsv"):
        os.symlink(os.path.abspath(os.path.join(datadir, "goi.tsv")),
                   os.path.abspath("goi.tsv"))

    # create directories - assume we are in the top-level
    # working directory
    to_make = [
        "RNA", 
        "DNA", 
        "compare_datasets", 
        "compare_datasets/scatterplot_genus_cog_fold.dir"]
    for directory in to_make:
        try:
            os.mkdir(directory)
        except OSError:
            continue

    # link files that are present in data directory

    for nt in ["RNA", "DNA"]:
        for data in glob.glob(os.path.abspath(os.path.join(datadir, nt))+"/*"):
            try:
                os.symlink(os.path.abspath(data),
                           os.path.abspath(os.path.join(nt, 
                                                        os.path.basename(data))))
            except OSError:
                continue


##############################################################
##############################################################
##############################################################


def loadCounts(scriptsdir):
    '''
    load Counts tables for genus and 
    for genes(NOGs)
    '''
    directories = ["RNA", "DNA"]
    filenames = ["genus.diamond.aggregated.counts.tsv.gz",
                 "gene_counts.tsv.gz"]
    for directory in directories:
        for filename in filenames:
            filename = os.path.abspath(os.path.join(directory, filename))
            outfile = P.snip(filename, ".tsv.gz") + ".load"
            table = P.toTable(outfile)
            database = os.path.abspath(os.path.join(directory, "csvdb"))
            statement = """zcat %(filename)s |
                           python %(scriptsdir)s/csv2db.py \
                                 --database-backend=sqlite \
                                 --retry \
                                 --database=%(database)s \
                                 --table=%(table)s \
                                 > %(outfile)s
                        """ % locals()
            os.system(statement)


##############################################################
##############################################################
##############################################################


def runMetagenomeSeq(cgatdir, 
                     rdir, 
                     filtered=False,
                     genus_level=False):
    '''
    run metagenomeSeq on counts tables
    '''
    directories = ["RNA", "DNA"]
    filenames = ["genus.diamond.aggregated.counts.tsv.gz",
                 "gene_counts.tsv.gz",
                 "genus.diamond.aggregated.counts.restricted.tsv.gz",
                 "gene_counts.restricted.tsv.gz",
                 "associated_taxa_counts.tsv.gz"]

    for directory in directories:
        for filename in filenames:
            filename = os.path.abspath(os.path.join(directory, filename))
            outfilename = P.snip(filename, ".gz")
            
            # check if file exists
            if not os.path.exists(filename):
                continue
            
            # only apply to genus level
            if genus_level:
                if "associated_taxa" not in filename:
                    continue
            else:
                if "associated_taxa" in filename:
                    continue

            # switch to skip re-running of 
            # previously run analysis
            if not filtered:
                if "restricted" in filename:
                    continue
            else:
                if "restricted" not in filename:
                    continue

            if "restricted" in filename:
                outprefix = P.snip(filename, ".restricted.tsv.gz")
            else:
                outprefix = P.snip(filename, ".tsv.gz")

            os.system("""gunzip -c \
                         %(filename)s > %(outfilename)s""" % locals())
            
            statement = ("""%(rdir)s/Rscript \
                         %(cgatdir)s/R/run_metagenomeseq.R --k 4 --a 0.1 -c %(outfilename)s -p %(outprefix)s""" % locals())
            os.system(statement)


##############################################################
##############################################################
##############################################################


def loadDiffTables(scriptsdir):
    '''
    load differential abundance tables for genus and 
    for genes(NOGs)
    '''
    directories = ["RNA", "DNA"]
    filenames = ["genus.diamond.aggregated.counts.diff.tsv",
                 "gene_counts.diff.tsv"]
    for directory in directories:
        for filename in filenames:
            filename = os.path.abspath(os.path.join(directory, filename))
            outfile = P.snip(filename, ".tsv") + ".load"
            table = P.toTable(os.path.basename(outfile))
            database = os.path.abspath(os.path.join(directory, "csvdb"))
            statement = """cd %(directory)s;
                           cat %(filename)s |
                           python %(scriptsdir)s/csv2db.py \
                                 --database-backend=sqlite \
                                 --retry \
                                 --database=%(database)s \
                                 --table=%(table)s \
                                 > %(outfile)s;
                           cd ../
                        """ % locals()
            os.system(statement)


##############################################################
##############################################################
##############################################################


def buildRestrictedSet(proj029scriptsdir, 
                       infile, 
                       restrictto):
    '''
    build filtered files for proper metagenomeseq
    analysis
    '''
    outfile = P.snip(infile, ".tsv.gz") + ".restricted.tsv.gz"
    logfile = infile + ".log"
    os.system("""zcat %(infile)s | python %(proj029scriptsdir)s/counts2restrictedcounts.py \
                                     --restrict-to=%(restrictto)s \
                                     --log=restrict.log \
                                     | gzip > %(outfile)s""" % locals())


##############################################################
##############################################################
##############################################################


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--datadir", dest="datadir", type="string",
                      help="provide data directory")

    parser.add_option("-r", "--rdir", dest="rdir", type="string",
                      help="provide R directory [default: %default]")

    parser.add_option("--cgatdir", dest="cgatdir", type="string",
                      help="provide cgat directory")

    parser.add_option("--proj029dir", dest="proj029dir", type="string",
                      help="provide proj029 directory")

    parser.set_defaults(rdir="/usr/bin")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    assert options.rdir, "must specify R location"
    assert options.datadir, "must specify data location"
    assert options.cgatdir, "must specify cgat location"
    assert options.proj029dir, "must specify proj029 location"

    # data and script directories
    datadir = options.datadir
    cgatdir = options.cgatdir
    proj029dir = options.proj029dir
    rdir = options.rdir
    cgatscriptsdir = "%s/scripts" % cgatdir
    proj029scriptsdir = "%s/scripts" % proj029dir

    # linking data from data directory
    E.info("creating directories and linking data")
    linkData(datadir)
    
    # load the counts tables
    E.info("loading counts tables into database")
    loadCounts(cgatscriptsdir)

    # initial normalisation of count data
    E.info("running initial metagenomeSeq normalisation")
    runMetagenomeSeq(cgatdir, rdir)

    # load differential abundance tables
    E.info("loading differential abundance tables")
    loadDiffTables(cgatscriptsdir)

    # build detection overlaps
    E.info("building genus-level overlaps between data sets")
    PipelineMetaomics.buildDetectionOverlap("RNA/genus.diamond.aggregated.counts.tsv.gz",
                                            "DNA/genus.diamond.aggregated.counts.tsv.gz",
                                            "compare_datasets/genus_overlap.tsv")

    E.info("building NOG-level overlaps between data sets")
    PipelineMetaomics.buildDetectionOverlap("RNA/gene_counts.tsv.gz",
                                            "DNA/gene_counts.tsv.gz",
                                            "compare_datasets/gene_overlap.tsv")

    # abundances of overlapping and distinct sets
    E.info("plotting abundances of overlapping and distinct genus sets")
    PipelineMetaomics.plotAbundanceLevelsOfOverlap("RNA/genus.diamond.aggregated.counts.tsv.gz",
                                                  "DNA/genus.diamond.aggregated.counts.tsv.gz",
                                                   "compare_datasets/genus_abundance_distributions.png")

    E.info("plotting abundances of overlapping and distinct gene sets")
    PipelineMetaomics.plotAbundanceLevelsOfOverlap("RNA/gene_counts.tsv.gz",
                                                   "DNA/gene_counts.tsv.gz",
                                                   "compare_datasets/gene_abundance_distributions.png",
                                                   of="genes")

    
    E.info("correlating abundance estimates bewteen genus abundances")
    PipelineMetaomics.scatterplotAbundanceEstimates("DNA/genus.diamond.aggregated.counts.norm.matrix",
                                                    "RNA/genus.diamond.aggregated.counts.norm.matrix",
                                                    "compare_datasets/genus_abundance_correlation.png")


    E.info("correlating abundance estimates between NOG abundances")
    PipelineMetaomics.scatterplotAbundanceEstimates("DNA/gene_counts.norm.matrix",
                                                    "RNA/gene_counts.norm.matrix",
                                                    "compare_datasets/gene_abundance_correlation.png")

    # build common sets of genera and NOGs
    E.info("building common genus set")
    PipelineMetaomics.buildCommonList("RNA/csvdb",
                                      "DNA/csvdb",
                                      "compare_datasets/common_genera.tsv")

    PipelineMetaomics.buildCommonList("RNA/csvdb",
                                      "DNA/csvdb",
                                      "compare_datasets/common_genes.tsv")
    

    directories = ["RNA", "DNA"]
    # build filtered sets for metagenomeSeq analysis
    E.info("build filtered sets (>0.1RPM in both data sets)")
    for directory in directories:
        buildRestrictedSet(proj029scriptsdir,
                           os.path.abspath(
                os.path.join(directory, "genus.diamond.aggregated.counts.tsv.gz")),
                           "compare_datasets/common_genera.tsv")
        buildRestrictedSet(proj029scriptsdir,
                           os.path.abspath(
                os.path.join(directory, "gene_counts.tsv.gz")),
                           "compare_datasets/common_genes.tsv")

    # re-rerun metagenomeSeq
    E.info("re-running metagenomeSeq on filtered data")
    runMetagenomeSeq(cgatdir, rdir, filtered = True)

    # load overwritten differential abundance tables
    E.info("loading differential abundance tables")
    loadDiffTables(cgatscriptsdir)

    # run principle components analysis
    E.info("running principle components analysis")
    filenames = ["genus.diamond.aggregated.counts.norm.matrix",
                 "gene_counts.norm.matrix"]
    for directory in directories:
        # change into correct directory
        os.chdir(os.path.abspath(directory))
        for filename in filenames:
            outfile = P.snip(filename, ".norm.matrix") + ".loadings.tsv"
            PipelineMetaomics.runPCA(filename, 
                                     outfile)
        # return to top level
        os.chdir(os.path.abspath("../"))
    
    E.info("plotting PCA loadings")
    for directory in directories:
        os.chdir(os.path.abspath(directory))
        filename = "genus.diamond.aggregated.counts.loadings.tsv"
        outfile = P.snip(filename, ".tsv") + ".png"
        PipelineMetaomics.plotPCALoadings(filename,
                                          outfile)
        os.chdir(os.path.abspath("../"))

    # predicting colitis-responsive NOGs
    E.info("building colitis-responsive NOG sets")
    PipelineMetaomics.buildRNADNARatio("DNA/gene_counts.diff.tsv",
                                       "RNA/gene_counts.diff.tsv",
                                       "compare_datasets/rna_dna_ratio.tsv")

    PipelineMetaomics.buildGeneDiffList("RNA/csvdb",
                                        "compare_datasets/common_genes.tsv",
                                        "compare_datasets/rna_diff_genes.tsv")

    PipelineMetaomics.buildGeneDiffList("DNA/csvdb",
                                        "compare_datasets/common_genes.tsv",
                                        "compare_datasets/dna_diff_genes.tsv")

    PipelineMetaomics.annotateRNADNARatio("compare_datasets/rna_dna_ratio.tsv",
                                          "compare_datasets/dna_diff_genes.tsv",
                                          "compare_datasets/rna_diff_genes.tsv",
                                          "compare_datasets/rna_dna_ratio.annotated.tsv")

    PipelineMetaomics.plotSets("compare_datasets/rna_dna_ratio.annotated.tsv",
                               "compare_datasets/rna_dna_ratio.annotated.png")


    PipelineMetaomics.buildGenesOutsidePredictionInterval("compare_datasets/rna_dna_ratio.annotated.tsv", 
                                                          "compare_datasets/rna_dna_ratio.annotated.outsidepi.tsv")

    # assigning NOGs to genera. Plot the maximum contribution
    # of genera to NOGs
    E.info("building proportion contribution matrix")
    PipelineMetaomics.buildGenusCogCountsMatrix("RNA/associated_ptaxa.tsv.gz", 
                                                "RNA/associated_ptaxa_average.matrix")
    
    E.info("plotting maximum genus contribution to NOG expression")
    PipelineMetaomics.plotMaxTaxaContribution("RNA/associated_ptaxa_average.matrix",
                                              "compare_datasets/rna_dna_ratio.annotated.outsidepi.tsv",
                                              "RNA/associated_ptaxa_max_contribution.png")

    
    E.info("plotting major genera attributed to colitis-responsive NOGS (up-regulates)")
    PipelineMetaomics.heatmapTaxaCogProportionMatrix("RNA/associated_ptaxa_average.matrix",
                                                     "compare_datasets/rna_dna_ratio.annotated.outsidepi.tsv",
                                                     "RNA/associated_ptaxa_heatmap.pdf")
    
    # run metagenomeSeq on per genus/NOG counts
    E.info("running metagenomeSeq on genus-NOG counts")
    runMetagenomeSeq(cgatdir, rdir, filtered=False, genus_level=True)

    E.info("plotting DNA vs. RNA fold changes for genes of interest")
    os.chdir(os.path.abspath("compare_datasets"))
    PipelineMetaomics.scatterplotPerCogTaxaDNAFoldRNAFold("../RNA/associated_taxa_counts.diff.tsv",
                                                          "../DNA/associated_taxa_counts.diff.tsv",
                                                          "../RNA/gene_counts.diff.tsv",
                                                          "../RNA/gene_counts.diff.tsv")
    os.chdir(os.path.abspath(".."))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
