"""
=======================================================
Filter fastq files based on reads unmapped in bam file
=======================================================

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

import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import CGATPipelines.PipelineMapping as PipelineMapping
#import CGATPipelines.PipelineMetagenomeAssembly as PipelineMetagenomeAssembly
import CGAT.FastaIterator as FastaIterator
import CGAT.Metaphlan as Metaphlan
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import pysam
import CGAT.Fastq as Fastq

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
###################################################################
###################################################################
@transform(glob.glob("*.fastq.gz"), 
           regex("(\S+).fastq.gz"),
           add_inputs(glob.glob("*.bam")),
           r"filtered.\1.fastq.gz")
def filterFastq(infiles, outfile):
    '''
    filter fastq file based on read tyep 
    speicifed in bam file
    '''
    t = PARAMS.get("reads_type")
    fastq = infiles[0]
    track = P.snip(fastq, ".fastq.gz")
    bams = infiles[1]
    if len(bams) == 1:
        bam = bams[0]
    else:
        bam = [bam for bam in bams if bam.find(track) != -1][0]

    job_options="-l mem_free=30G"
    
    if PARAMS.get("reads_invert") == 1:
        invert = "--invert"
    else:
        invert = ""
    statement = '''zcat %(fastq)s | python /ifs/projects/proj029/src/scripts/fastq2filteredfastq.py
                                    -b %(bam)s
                                    %(invert)s
                                    --reads=%(t)s
                                    --log=%(outfile)s.log
                                    | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
