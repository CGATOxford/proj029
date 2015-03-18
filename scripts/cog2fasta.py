'''
cog2fasta.py
=============

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a COG identifier as input and outputs the fasta sequence for all
genes in the Integrated Gene Catalogue (IGC) that are associated with that particular
COG.

Usage
-----

Example::

   python cog2fasta.py 

Type::

   python cog2fasta.py

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
import CGAT.FastaIterator as FastaIterator

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-m", "--gene2cog", dest="gene2cog", type="string",
                      help="provide gene to COG mapping file"  )
    parser.add_option("-c", "--cog", dest="cog", type="string",
                      help="provide COG identifier (genes will be returned that map to this COG)")
    parser.add_option("-f", "--fasta", dest="fasta", type="string",
                      help="provide full fasta file from which to retrieve sequences")
    parser.add_option("--output-list", dest="output_list", action="store_true",
                      help="output gene names associated with COG?")

    parser.set_defaults(gene2cog = None,
                        cog = None,
                        fasta = None,
                        output_list=False)

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    assert options.gene2cog and options.cog and options.fasta, "must specify all command line options"
        
    E.info("retrieving gene names associated with %s" % options.cog)
    genes = set()

    c = E.Counter()
    for line in IOTools.openFile(options.gene2cog).readlines():
        data = line[:-1].split("\t")
        gene, cog = data[0], data[1]
        if cog == options.cog:
            c.ngenes += 1
            genes.add(gene)

    E.info("retrieving gene sequences")
    if options.output_list:
        outf = IOTools.openFile("gene_names.tsv.gz", "w")
    else:
        outf = None
    for fasta in FastaIterator.iterate(IOTools.openFile(options.fasta)):
        identifier = fasta.title.split(" ")[0]
        if identifier in genes:
            options.stdout.write(">%s\n%s\n" % (identifier, fasta.sequence))
            if outf:
                outf.write("%s\n" % identifier)
    if outf:
        outf.close()

    E.info("%s" % c)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
