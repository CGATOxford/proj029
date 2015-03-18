'''
metaphlan2species.py - template for CGAT scripts
====================================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

create list of species from output of metaphlan2table.py --relab

Usage
-----

Example::

   python metaphlan2species.py 

Type::

   python metaphlan2species.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    header = options.stdin.readline()
    for line in options.stdin.readlines():
        data = line[:-1].split("\t")
        taxa, species, relab = data[0], data[1], data[2]
        if taxa == "species":
            options.stdout.write("%s\n" % species)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
