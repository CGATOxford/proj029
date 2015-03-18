'''
genematrix2sample.py
=====================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take a matrix with gene counts and subsample to a specified number.


Usage
-----

Example::

   python genematrix2sample.py 

Type::

   python genematrix2sample.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import random
import collections
import pandas
import itertools

import CGAT.Experiment as E


# function definitions

###################################################
###################################################
###################################################

def nCols(infile):
    '''
    get number of count columns
    '''
    return len(infile.readline().split("\t")) - 1

###################################################
###################################################
###################################################

def readData(infile):
    '''
    read pandas data frame object
    '''
    data = pandas.read_csv(infile, sep = "\t", index_col = 0)
    return data

###################################################
###################################################
###################################################

def buildPopulations(data):
    '''
    build "bag" of gene names based on counts 
    in the data. Return list of lists - one
    for each column in the data
    '''
    # read in data as a data frame
    E.info("building populations to sample from")
    populations = []
    for i in range(len(data.columns)):
        population = []
        for ind in data.index:
            gene = ind
            population += itertools.repeat(gene, data.loc[ind][i])
        populations.append(population)
    E.info("built populations")
    return populations

###################################################
###################################################
###################################################

def sample(population, sample):
    '''
    sample from the population
    '''
    E.info("subsampling from population")
    subsample = random.sample(population, sample)
    E.info("finished subsampling")
    return subsample

###################################################
###################################################
###################################################

def buildCounts(subsample):
    '''
    build a dictionary of counts for each
    COG/gene
    '''
    counts = collections.defaultdict(int)
    E.info("counting sampled data")
    for gene in subsample:
        counts[gene] += 1
    E.info("finished counting")
    return counts

###################################################
###################################################
###################################################

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-s", "--sample", dest="sample", type="int",
                      help="number to sample to")  

    # set default options
    parser.set_defaults(sample = 10000000)

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    infile = options.stdin

    header = infile.readline().split("\t")
    result = {}
    populations = buildPopulations(readData(infile))
    for i in range(len(populations)):
        counts = buildCounts(sample(populations[i], options.sample))
        result[header[i]] = counts
    
    data = pandas.DataFrame(result)
    data.to_csv(options.stdout, na_rep = 0, sep = "\t")

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

