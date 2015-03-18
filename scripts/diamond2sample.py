'''
diamond2sample.py
==================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Pseudo-random subsampling of diamond alignments

Usage
-----

.. Example use case

Example::

   python diamond2sample.py

Type::

   python diamond2sample.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import random
import CGAT.Diamond as Diamond
import CGAT.IOTools as IOTools

def countAlignments(infile):
    '''
    count total number of alignments in file
    '''
    total = 0
    for alignment in Diamond.best_alignment_iterator( \
        Diamond.query_iterator( \
            Diamond.alignment_iterator( \
                IOTools.openFile(infile)))):
        total += 1
    return total

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-a", "--alignment-file", dest="alignment_file", type="string",
                      help="supply alignment file to sample from")
    parser.add_option("-s", "--sample", dest="sample", type="int",
                      help="nuber of alignments to output")

    parser.set_defaults(sample = 10000)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    
    E.info("counting alignments")
    total = countAlignments(options.alignment_file)
    E.info("counted %i alignments" % total)

    # which rows to take
    E.info("goint to sample %i alignments" % options.sample)
    take = random.sample(range(1, total + 1), options.sample)

    c = 1
    E.info("sampling")
    for alignment in Diamond.best_alignment_iterator( \
        Diamond.query_iterator( \
            Diamond.alignment_iterator( \
                IOTools.openFile(options.alignment_file)))):
        if c in take:
            options.stdout.write("\t".join([alignment.qid, 
                                            alignment.gi, 
                                            alignment.ref, 
                                            alignment.ngaps,
                                            alignment.length, 
                                            alignment.evalue, 
                                            alignment.nmismatches, 
                                            alignment.identity, 
                                            alignment.score]) + 
                                 "\n")
        c += 1

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
