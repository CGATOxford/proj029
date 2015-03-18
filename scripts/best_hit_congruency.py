'''
best_hit_congruency.py
=======================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Take a file of reads and test for congruency of hits in terms
of COG annotations.

Usage
-----

.. Example use case

Example::

   python best_hit_congruency.py

Type::

   python best_hit_congruency.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import random
import CGAT.Diamond as Diamond
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

    parser.add_option("-r", "--reads", dest="reads", type="string",
                      help="supply list of reads")

    parser.add_option("-m", "--map", dest="map", type="string",
                      help="gene2cog mapping file")

    parser.add_option("-a", "--alignment-file", dest="alignment_file", type="string",
                      help="alignment file")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    
    E.info("reading mapping file %s" % options.map)
    gene2cog={}
    for line in IOTools.openFile(options.map):
        gene, cog = line[:-1].split("\t")
        gene2cog[gene] = cog
    E.info("loaded mapping file")
        
    E.info("loading reads")
    reads = set([x[:-1] for x in IOTools.openFile(options.reads).readlines()])

    E.info("loaded reads")

    E.info("iterating over reads and testing congruency")
    options.stdout.write("read\tnhits\tncongruent\tpcongruent\n")
    for alignments in Diamond.query_iterator( \
        Diamond.alignment_iterator( \
            IOTools.openFile(options.alignment_file))):
        
        # skip those not in reads list
        if alignments[0].qid not in reads: continue

        # total number of hits
        n_hits = float(len(alignments))
        if n_hits == 1.0:
            continue

        # get the best alignment
        scores = []
        for alignment in alignments:
            scores.append(alignment.score)
        best = max(scores)
        best_alignments = [x for x in alignments if x.score == best]
        
        # if there are more than one best hit
        if len(best_alignments) > 1:
            best_alignments = random.sample(best_alignments, 1)
        best_alignment = best_alignments[0]

        # get the COG associated with the alignment 
        best_ref = best_alignment.ref
        best_cog = gene2cog[best_alignment.ref]
        n_congruent = float(len([x.ref for x in alignments if gene2cog[x.ref] == best_cog and x.ref != best_ref]))
        p_congruent = (n_congruent/(n_hits-1))*100

        options.stdout.write("%s\t%i\t%i\t%f\n" % (alignments[0].qid, n_hits, n_congruent, p_congruent))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
