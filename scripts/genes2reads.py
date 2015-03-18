'''
genes2reads.py
===============

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

From a list of genes and alignment file return all read ids mapping to genes.

Usage

-----

Example::

   python genes2reads.py 

Type::

   python genes2reads.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Diamond as Diamond

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-g", "--gene-list", dest="gene_list", type="string",
                      help="provide gene list"  )
    parser.add_option("-a", "--alignment-file", dest="alignment_file", type="string",
                      help="provide alignment file"  )


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    E.info("loading gene list")
    genes = set([x[:-1] for x in IOTools.openFile(options.gene_list).readlines()])
    E.info("loaded gene list")

    E.info("reading alignments and returning results")
    for alignment in Diamond.best_alignment_iterator(Diamond.query_iterator(Diamond.alignment_iterator(IOTools.openFile(options.alignment_file)))):
        if alignment.ref in genes:
            options.stdout.write("%s\n" % alignment.qid)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
