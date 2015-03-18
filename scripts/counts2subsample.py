'''
counts2subsample.py
====================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

subsample a count table of form:

    ref   count
    COG1  23
    COG2  1000

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
import numpy as np
import collections
import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-s", "--sample", dest="sample", type="int",
                      help="subsample to this number"  )

    parser.set_defaults(sample = 1000000)

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    E.info("getting probabilities")
    options.stdin.readline()
    counts = []
    ids = []
    for line in options.stdin:
        data = line[:-1].split("\t")
        counts.append(int(data[1]))
        ids.append(data[0])

    probs = [float(x)/sum(counts) for x in counts]
    E.info("loaded probabilities")

    E.info("sampling %i random samples" % options.sample)
    samples = np.random.choice(ids, size = options.sample, p = probs, replace = True)
    E.info("finished sampling")

    result = collections.defaultdict(int)
    for s in samples:
        result[s] += 1

    E.info("writing results")
    options.stdout.write("ref\tcount\n")
    for s, count in result.iteritems():
        options.stdout.write("%s\t%s\n" % (s, count))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
