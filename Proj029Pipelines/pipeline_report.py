"""
================
Analysis update
================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python


This is a report that monitors the progress of the data analysis. Many different analysis
steps are run - involving multiple different pipelines. As such there is not a single pipeline
report that covers all of the areas of analysis. Therefore this report is an overview of each
of the processing steps and data analysis procedures that have been performed - pulling in 
example plots and data from the results directories of other pipelines.

"""

import os, sys, re
import glob
from ruffus import *
import CGAT.Experiment as E
import CGATPipelines.PipelineMapping as PipelineMapping

# pipeline configuration
import CGATPipelines.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'annotations_dir': "",
        'paired_end': False})

PARAMS = P.PARAMS

#########################################################################
#########################################################################
#########################################################################
# Get the raw read counts
#########################################################################
DATADIR_RAW = PARAMS["input_raw"]
SEQUENCESUFFIXES = ("stool*.fastq.1.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR_RAW, suffix_name)
                      for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r".*/(\S+).(fastq.1.gz|sra|csfasta.gz|csfasta.F3.gz|export.txt.gz)")

#########################################################################
#########################################################################
#########################################################################

@follows(mkdir("nreads.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"nreads.dir/\1.nreads")
def countReads(infile, outfile):
    '''count number of reads in input files.'''
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run()

#########################################################################
#########################################################################
#########################################################################

@merge(countReads, "nreads.dir/reads_summary.load")
def loadReadCounts(infiles, outfile):
    ''' 
    load read counts summary
    '''
    tempname = "reads_summary.tsv"
    temp = open(tempname, "w")
    for infile in infiles:
        track = P.snip(os.path.basename(infile), ".nreads")
        nreads = open(infile).readline().split("\t")[1]
        temp.write("%s\t%s\n" % (track, nreads))
    temp.close()
    P.load(tempname, outfile, "--header=track,nreads")


#########################################################################
#########################################################################
#########################################################################

@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


