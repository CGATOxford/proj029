'''

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

combine marker genes and their clades with fasta information - cpg content etc

Usage
-----

Example::

   python cgat_script_template.py 

Type::

   python cgat_script_template.py --help

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

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-m", "--marker-file", dest="marker_file", type="string",
                      help="supply marker file with marker to clade annotations")

    parser.add_option("-f", "--fasta-table", dest="fasta_table", type="string",
                      help="supply table with annotation information")

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    marker2clade = {}
    for line in IOTools.openFile(options.marker_file).readlines():
        data = line[:-1].split("\t")
        level, clade = data[1].split("__")
        marker2clade[data[0]] = [level, clade]

    inf = IOTools.openFile(options.fasta_table)
    header = inf.readline()[:-1] + "\tlevel\tclade\n"
    options.stdout.write(header)
    for line in inf.readlines():
        data = line[:-1].split("\t")
        if data[0] == "total": continue
        options.stdout.write(line[:-1] + "\t" + "\t".join(marker2clade[data[0]]) + "\n")

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

