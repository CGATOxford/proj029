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

   python run_analysis.py --rdir=/user/bin/R

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
                                 --backend=sqlite \
                                 --retry \
                                 --database=%(database)s \
                                 --table=%(table)s \
                                 > %(outfile)s
                        """ % locals()
            os.system(statement)


##############################################################
##############################################################
##############################################################


def runMetagenomeSeq(cgatdir, rdir):
    '''
    run metagenomeSeq on counts tables
    '''
    directories = ["RNA", "DNA"]
    filenames = ["genus.diamond.aggregated.counts.tsv.gz",
                 "gene_counts.tsv.gz"]
    for directory in directories:
        for filename in filenames:
            filename = os.path.abspath(os.path.join(directory, filename))
            outfilename = P.snip(filename, ".gz")
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
            table = P.toTable(outfile)
            database = os.path.abspath(os.path.join(directory, "csvdb"))
            statement = """cat %(filename)s |
                           python %(scriptsdir)s/csv2db.py \
                                 --backend=sqlite \
                                 --retry \
                                 --database=%(database)s \
                                 --table=%(table)s \
                                 > %(outfile)s
                        """ % locals()
            os.system(statement)


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

    parser.set_defaults(rdir="/usr/bin")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # data and script directories
    datadir = options.datadir
    cgatdir = options.cgatdir
    rdir = options.rdir
    cgatscriptsdir = "%s/scripts" % cgatdir

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


    E.info("plotting abundances of overlapping and distinct genus sets")
    PipelineMetaomics.plotAbundanceLevelsOfOverlap("RNA/gene_counts.tsv.gz",
                                                   "DNA/gene_counts.tsv.gz",
                                                   "compare_datasets/gene_abundance_distributions.png")

    
    E.info("correlating abundance estimates bewteen genus abundances")
    PipelineMetaomics.scatterplotAbundanceEstimates("DNA/genus.diamond.aggregated.counts.norm.matrix",
                                                    "RNA/genus.diamond.aggregated.counts.norm.matrix",
                                                    "compare_datasets/genus_abundance_correlation.png")


    E.info("correlating abundance estimates bewteen NOG abundances")
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
    




    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
