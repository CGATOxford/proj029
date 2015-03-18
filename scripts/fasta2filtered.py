'''
fasta2filtered.py
=======================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

filter fasta file based on fileof seq ids

Usage
-----

Example::


Type::


for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import CGAT.Pipeline as P
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator

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

    parser.add_option("-f", "--filename", dest="filename",
                      help="filename with ids")


    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    ids = set([x[:-1] for x in open(options.filename).readlines()])
    for fasta in FastaIterator.iterate(options.stdin):
        if fasta.title in ids:
            options.stdout.write(">%s\n%s\n" % (fasta.title, fasta.sequence))
            

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
