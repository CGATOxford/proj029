'''
counts2restrictedcounts.py
===========================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take a counts table and filter it based on features in a file

Take a 

Usage
-----

.. Example use case

Example::

   python counts2restrictedcounts.py

Type::

   python counts2restrictedcounts.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-r", "--restrict-to", dest="restrict_to", type="string",
                      help="supply file of features to restrict counts table to")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    features = set([x[:-1] for x in open(options.restrict_to).readlines()])
    header = options.stdin.readline()
    options.stdout.write(header)
    for line in options.stdin.readlines():
        data = line[:-1].split("\t")
        if data[0] in features:
            options.stdout.write(line)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
