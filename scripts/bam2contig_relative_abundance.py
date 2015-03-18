'''
bam2contig_relative_abundance.py
=================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------
Take a bam file and output a tab delimited text file of all contigs and the
% of reads that map to each.

Usage
-----

Example::

   python bam2contig_relative_abundance.py 

Type::

   python bam2contig_relative_abundance.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import pysam

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-b", "--bamfile", dest="bamfile", type="string",
                      help="supply bam file"  )
    parser.add_option("-n", "--nreads", dest="nreads", type="int",
                      help="supply number of mapped reads"  )
    parser.add_option("-l", "--scale-length", dest="scale_length", action="store_true",
                      help="scale by contig length - i.e. per Kb"  )


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    samfile = pysam.Samfile(options.bamfile)
    options.stdout.write("contig\tpct_mapped\ttotal_mapped\tfpkm\n")
    for contig, length in zip(samfile.references, samfile.lengths):
        c = 0
        for alignment in samfile.fetch(contig):
            c += 1
        pct_mapped = (float(c)/options.nreads)*100
        fpkm = c/((length/1000)/(float(options.nreads)/1000000))
        total_mapped = c
        options.stdout.write("%s\t%f\t%i\t%f\n" % (contig, pct_mapped, total_mapped, fpkm))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
