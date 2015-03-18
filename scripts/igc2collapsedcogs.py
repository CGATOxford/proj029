'''
igc2collapsedcogs.py
=====================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Collapse genes from the IGC into a single representative COG based on
expression levels in a matrix

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
import numpy as np
import CGAT.FastaIterator as FastaIterator

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-a", "--annotations", dest="annotations", type="string",
                      help="supply IGC annotation file"  )
    parser.add_option("-e", "--expression-matrix", dest="expression_matrix", type="string",
                      help="supply expression matrix file"  )
    parser.add_option("-f", "--fasta", dest="fasta", type="string",
                      help="supply expression matrix file"  )
    parser.add_option("-m", "--map", dest="map", action="store_true",
                      help="output the final gene to COG map as a text file"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    c = E.Counter()
    E.info("loading annotations")
    gene2cog = {}
    for line in IOTools.openFile(options.annotations).readlines():
        c.input_annotations += 1
        data = line[:-1].split("\t")
        gene, cog = data[1], data[8]

        # ignore genes that have an undefined COG function
        if cog == "unknown": continue

        # we remove genes that have > 1 COG identifier
        # These represent ~1.5% of the database and 
        # we expect to be able to align to another 
        # gene with any one of the given identifiers
        if len(cog.split(";")) > 1: 
            continue

        gene2cog[gene] = cog
        c.filtered_annotations += 1
    E.info("annotations loaded")

    E.info("calculating max COG expression")
    mat = open(options.expression_matrix)
    h = mat.readline()
    result = {}
    for line in mat.readlines():
        c.input_genes += 1
        data = line[:-1].split("\t")
        values, gene = map(float,data[:-1]), data[-1].strip('"')

        try:
            cog = gene2cog[gene]
        except KeyError:
            continue

        # get average abundance for gene
        ave = np.mean(values)
        if cog not in result.keys():
            result[cog] = (gene, ave)
        else:
            if ave > result[cog][1]:
                result[cog] = (gene, ave)
            else:
                result[cog] = result[cog] 

    E.info("max COG found")
    genes = set()
    for gene_ave in result.values():
        genes.add(gene_ave[0])

    E.info("resolving gene names and writing fasta file")
    for fasta in FastaIterator.iterate(IOTools.openFile(options.fasta)):
        c.input_fasta += 1
        gene = fasta.title.split(" ")[0]
        if gene in genes:
            options.stdout.write(">%s\n%s\n" % (fasta.title, fasta.sequence))
        else:
            continue

    if options.map:
        outf = IOTools.openFile("gene2cog.map.gz", "w")
        for gene in genes:
            outf.write("%s\t%s\n" % (gene, gene2cog[gene]))
        outf.close()
            
    c.output_genes = len(genes)

    E.info("%s" % c)
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
