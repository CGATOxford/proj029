'''
species2multifasta.py - template for CGAT scripts
====================================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Taxonomic profiling is a general first step in a metagenomics analysis. However, the majority
of reads do not contribute to the profile. In addition, while we could align all of our reads
to the entire set of reference genomes that are currently available this is likely to be
time consuming and often redundant. Therefore, to determine the number of reads that align to
known genomes, it is possible to use the taxonomic profiling analysis to subset the number of
genomes we will align to. This script takes a list of species produced from a taxonomic profile
and outputs a single multifasta file containing genomes of the listed species.

species names must be separated by an underscore e.g heliocbacter_hepaticus

Usage
-----

Example::

   python species2multifasta.py --species=species_list.tsv --genomes-dir=genomes.dir

Type::

   python species2multifasta.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import CGAT.FastaIterator as FastaIterator
import glob
import collections
import CGAT.IOTools as IOTools

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-s", "--species", dest="species", type="string",
                      help="list of species to build genomes file"  )
    parser.add_option("-g", "--genomes-dir", dest="genomes_dir", type="string",
                      help="directory containing genomes. This directory follows the NCBI structure \
                            each species has its own subdirectory"  )
    parser.add_option("--avoid-plasmids", dest="avoid_plasmids", action="store_true",
                      help="use the name of the genome sequence to avoid plasmids in the analysis")
    parser.add_option("-l", "--level", dest="level", type="choice", choices=("species", "genus"),
                      help="return all species genomes in a genus, or just the species")
    parser.add_option("--masked", action = "store_true",
                      help="should masked sequences be returned - uses ensembl nomenclature")

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # strip newline character for species
    if options.level == "species":
        species_list = set([x[:-1] for x in open(options.species).readlines()])
    elif options.level == "genus":
        species_list = set([x[:-1].split("_")[0] for x in open(options.species).readlines()])

    if options.masked:
        pattern = "*dna_rm.genome*"
    else:
        pattern = "*dna.genome*"

    for species in species_list:
        E.info("checking for %s in genomes directory" % species)
        for species_path, dirs, files in os.walk(options.genomes_dir):
            if os.path.basename(species_path) != "dna": continue
            if species_path.find(species) != -1:
                genome_file = glob.glob(os.path.join(species_path, pattern))
                E.info("adding genome %s to multifasta" % ",".join(genome_file))
                if len(genome_file) != 1:
                    #E.warn("more than one genome found in directory %s: avoiding plasmids" % species_path)
                    for g in genome_file:
                        identifier = IOTools.openFile(g).readline()[0]
                        if identifier.find("plasmid") != -1: continue
                        genome_file = g
                else:
                    genome_file = genome_file[0]
                if options.avoid_plasmids:
                    if IOTools.openFile(genome_file).readline().find("plasmid") != -1: continue
                for line in IOTools.openFile(genome_file).readlines():
                    options.stdout.write(line)
        E.info("finished writing genomes for %s" % species)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
