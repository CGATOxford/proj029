'''
fastq2filteredfastq.py
=======================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

filter a fastq file based on reads in a bam file.

Usage
-----

Example::

   zcat in.fastq.gz | python fastq2filteredfastq.py --bamfile=x.bam

Type::

   python fastq2filteredfastq.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import pysam
import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq

import CGAT.Experiment as E


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-b", "--bamfile", dest="bamfile", type="string",
                      help="input bamfile to filter reads from")
    parser.add_option("-r", "--reads", dest="reads", type="choice",
                      choices=("mapped", "unmapped"),
                      help="type of read to keep")
    parser.add_option("-s", "--scriptsdir", dest="scriptsdir", type="string",
                      help="CGAT scripts directory")
    parser.add_option("-i", "--invert", dest="invert", action="store_true",
                      help="invert selection - if for example unmapped reads \
                            aren't output")


    parser.set_defaults(bamfile = None,
                        reads = "mapped",
                        scriptsdir = "/ifs/devel/nicki/cgat_git/cgat/scripts",
                        invert = False)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    c = E.Counter()
    c.input_alignments = 0
    c.input_reads = 0
    c.output_reads = 0

    # output text file for reads TO KEEP
    bam = pysam.Samfile(options.bamfile, "rb")
    temp = P.getTempFile(".")
    E.info("iterating over bam file")

    for alignment in bam.fetch(until_eof = True):
        c.input_alignments += 1
        if options.reads == "unmapped":
            if alignment.is_unmapped:
                #c.input_alignments += 1
                temp.write(alignment.qname + "\n")
        elif options.reads == "mapped":
            if not alignment.is_unmapped:
                #c.input_alignments += 1
                temp.write(alignment.qname + "\n")
    temp.close()            
    
    tempname = temp.name

    E.info("filtering fastq file")
    # filter fastq file
    ids = set(IOTools.readList(IOTools.openFile(tempname).readlines()))
    c.input_alignments = len(ids)
    for fastq in Fastq.iterate(options.stdin):
        c.input_reads += 1
        if fastq.identifier.endswith("/1") or fastq.identifier.endswith("/2"):
            identifier = fastq.identifier[:-2]
        elif len(fastq.identifier.split(" ")) == 2:
            identifier = fastq.identifier.split(" ")[0]
        else:
            identifier = fastq.identifier
        if not options.invert:
            if identifier in ids: 
                c.output_reads += 1
                options.stdout.write("%s\n" % fastq)
        else:
            if identifier in ids: continue
            c.output_reads += 1
            options.stdout.write("%s\n" % fastq)
            
    E.info(c)

    os.unlink(tempname)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
