'''
GSEAPrepranked.py
==================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Perform gene set enrichment analysis on a pre-reanked list.

Usage
-----

.. Example use case

Example::

   python GSEAPreranked.py

Type::

   python GSEAPreranked.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import random
import collections
import tempfile
import CGAT.Experiment as E
from rpy2.robjects import r as R

#############################################
#############################################
#############################################

def readGeneset(geneset):
    '''
    return set of identifiers in gene set
    '''
    geneset = set([x[:-1] for x in open(geneset).readlines()])
    return geneset

#############################################
#############################################
#############################################

def readRankedList(rankedList):
    '''
    read in ranked list
    '''
    rl = open(rankedList)
    ranked_list = collections.OrderedDict()
    for line in rl.readlines():
        identifier, fold = line[:-1].split("\t")
        ranked_list[identifier] = float(fold)
    return ranked_list

#############################################
#############################################
#############################################

def calculateEs(ranked_list, geneset):
    '''
    calculate the enrichment score
    '''
    ng = len(geneset)
    nr = 0.0
    nl = 0.0
    es = []
    for identifier, fold in ranked_list.iteritems():
        fold = float(fold)
        nl += 1
        if identifier in geneset:
            nr += abs(fold)**1

    phit = 0
    pmiss = 0
    position_in_list = []
    i = 1
    for identifier, fold in ranked_list.iteritems():
        i += 1
        fold = float(fold)
        if identifier in geneset:
            phit += abs(fold)**1/nr
            position_in_list.append(i)
        else:
            pmiss += 1/(nl - ng)
        es.append(phit-pmiss)
    return max(es), es, position_in_list

#############################################
#############################################
#############################################

def shuffleRankedList(ranked_list):
    '''
    shuffle ranks
    '''
    x = ranked_list.keys()
    random.shuffle(x)
    return collections.OrderedDict(zip(x, ranked_list.values()))

#############################################
#############################################
#############################################

def calculateSignificance(observed, permuted):
    '''
    permutation p-value
    '''
    return float(len([x for x in permuted if x >= observed]))/len(permuted)

#############################################
#############################################
#############################################

def buildLeadingEdgeSet(ranked_list, running_es, geneset):
    '''
    return the set of features that contribute
    to the enrichment i.e. before the max(es)
    '''
    leading_edge = []
    ind_at_max = [i for i in range(len(running_es)) \
                  if running_es[i] == max(running_es)][0]
    for i in range(len(running_es)):
        feature = ranked_list.keys()[i]
        if i <= ind_at_max and feature in geneset:
            leading_edge.append(ranked_list.keys()[i])
    return leading_edge

#############################################
#############################################
#############################################

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-l", "--ranked-list", dest="ranked_list", type="string",
                      help="supply ranked list with names and values")
    parser.add_option("-g", "--geneset", dest="geneset", type="string",
                      help="geneset to test for enrichment")
    parser.add_option("-n", "--nperm", dest="nperm", type="int",
                      help="number of permutations to compute p-value")
    parser.add_option("-o", "--outfile-prefix", dest="outfile_prefix", type="string",
                      help="prefix for output plots")

    parser.set_defaults(nperm=10000)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    ranked_list, geneset, nperm = readRankedList(options.ranked_list), \
                                  readGeneset(options.geneset), options.nperm
    observed, running_es, position_in_list = calculateEs(ranked_list, geneset)

    # wrrite out leading edge subset
    outf = open("%s_leadingedge.tsv" % options.outfile_prefix, "w")
    for l in buildLeadingEdgeSet(ranked_list, running_es, geneset):
        outf.write("%s\n" % l)
    outf.close()

    permuted = []
    for i in range(nperm):
        shuffled = shuffleRankedList(ranked_list)
        es, running_es_permuted, position_permuted = calculateEs(shuffled, geneset)
        permuted.append(es)

    # calculate significance
    p = calculateSignificance(observed, permuted)
    sys.stdout.write("p = %f\nES = %f" % (p, observed))

    # plot the running ES scores
    temp = tempfile.NamedTemporaryFile(dir=".", delete=False)
    for es in running_es:
        temp.write(str(es) + "\n")
    temp.close()

    R('''library(ggplot2)''')
    R('''es <- read.csv("%s", header=F, sep="\t", stringsAsFactors=F)''' % temp.name)
    R('''enrichment <- data.frame(cbind(seq(1, nrow(es), 1), es[,1]))''')
    R('''colnames(enrichment) <- c("rank", "ES")''')
    R('''rank.in.list <- data.frame(x=c(%s), y=-0.25*max(enrichment$ES))''' % ",".join(map(str,position_in_list)))

    R('''plot1 <- ggplot(enrichment, aes(x=rank, y=ES)) + geom_line(size=2)''')
    R('''plot2 <- plot1 + geom_segment(data=rank.in.list, 
                                       aes(x=x, 
                                           y=0, 
                                           xend=x, 
                                           yend=y, 
                                           alpha=0.5), 
                                           col="red")''')
    R('''plot2 + geom_vline(xintercept=0) + geom_hline(yintercept=0)''')
    outname = options.outfile_prefix + "_running_es.pdf"
    R('''ggsave("%s")''' % outname)
    os.unlink(temp.name)

    # plot the significance - permuted values histogram 
    temp = tempfile.NamedTemporaryFile(dir=".", delete=False)
    for es in permuted:
        temp.write(str(es) + "\n")
    temp.close()

    R('''dat <- read.csv("%s", header=F, sep="\t", stringsAsFactors=F)''' % temp.name)
    R('''plot1 <- ggplot(dat, aes(x=V1)) + geom_histogram(colour="grey")''')
    R('''obs <- data.frame(x=%f, y=0)''' % observed)
    R('''plot1 + geom_point(data=obs, aes(x=x, y=y), size=5, col="red")''')
    outname = options.outfile_prefix + "_histogram.pdf"
    R('''ggsave("%s")''' % outname)
    os.unlink(temp.name)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
