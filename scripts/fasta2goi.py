'''
fasta2goi.py
=============

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

get fasta sequences for genes of interest specified 
as a list

Usage
-----

Example::

   python fasta2goi.py

Type::

   python fasta2goi.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.FastaIterator as FastaIterator


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-g", "--genes", dest="genes", type="string",
                      help="supply list of gene ids to extract"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    gois = options.genes.split(",")

    for fasta in FastaIterator.iterate(options.stdin):
        header = fasta.title.split(":")
        gene_id = header[7].split(" ")[0]
        if len(header) == 9:
            description = "NA"
        else:
            description = header[9].strip('"')
        if gene_id in gois:
            options.stdout.write(">gene_id: %s description: %s\n%s\n" % (gene_id, description, fasta.sequence))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


