'''
metaref2gtf.py
===============

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take as input the metaref gene fasta file and output a faux gtf file:
The gtf file is faux in the sense that each gene is assigned to its 
own artificial chromosome with coordinates ranging from 1-length(gene).

Usage
-----

Example::

   python metaref2gtf.py

Type::

   python metaref2gtf.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools

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

    parser.add_option("-p", "--protein-file", dest="protein_file", type="string",
                      help="supply protein file - used to annotate gene names, clades etc")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    assert options.protein_file, "must supply protein file"

    # get gene id to gene name, species name map
    id2attributes = {}

    E.info("getting attributes for each gene in %s" % options.protein_file)
    for fasta in FastaIterator.iterate(IOTools.openFile(options.protein_file)):
        identifier = fasta.title
        gene_id = str(identifier.split("|")[1].split(" ")[0]).strip()
        gene_name = identifier.split("|")[0]
        description = " ".join(identifier.split("|")[1].split(" ")[1:]).split("[")[0]
        species = identifier.split("[")[1].split("]")[0]
        id2attributes[gene_id] = (gene_id, gene_name, description, species)
    E.info("finished reading attributes")

    E.info("iterating over input and creating gtf entries")
    
    for fasta in FastaIterator.iterate(options.stdin):
        gene = str(fasta.title)
        try:
            chrom = "chr%s" % gene
            start, end = 1, len(fasta.sequence)
            gene_id = "gene_id %s" % id2attributes[gene][0]
            gene_name = "gene_name %s" % id2attributes[gene][1]
            description = "description %s" % id2attributes[gene][2]
            species = "species %s" % id2attributes[gene][3]
            attributes = "; ".join([gene_id, gene_name, description, species])
        except KeyError:
            E.warn("no matching annotation for gene %s" % gene)
            continue

        # strand is alway +
        options.stdout.write("\t".join([chrom, 
                                        "metaref", 
                                        "gene", 
                                        str(start), 
                                        str(end), 
                                        ".", 
                                        "+", 
                                        "0", 
                                        attributes]) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
