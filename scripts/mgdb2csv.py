'''
mgdb2csv.py
============

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

mgdg stands for metagenome database. This script therefore outputs
data from a central database created by pipeline_metagenomeassembly.py
of vatious kinds.

Usage
-----

Example::

   python mgdb2csv.py -d <database> --data=coverage,gc

Type::

   python mgdb2csv.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import gzip
import sqlite3
import CGAT.Experiment as E

def connect(database):
    '''
    connect to database
    '''
    dbh = sqlite3.connect(database)
    cc = dbh.cursor()
    return cc


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """



    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-d", "--database", dest="database", type="string",
                      help="supply database created by pipeline_metagenomeassembly.py"  )
    parser.add_option("--data", dest="data", type="choice",
                      choices=("coverage", "essential_genes"), help="what data is required"  )
    parser.add_option("--taxa-level", dest="taxa_level", type="choice",
                      choices=("genus", "family", "order", "class", "phylum"), help="what data is required"  )
    parser.add_option("-a", "--assembler", dest="assembler", type="choice",
                      choices=("ray", "metavelvet", "idba", "sga", "soapdenovo"), help="what data is required"  )
    parser.add_option("--track-assembly", dest="track_assembly", type="string",
                      help="track name for contigs"  )
    parser.add_option("--track-one", dest="track_one", type="string",
                      help="first coverage track"  )
    parser.add_option("--track-two", dest="track_two", type="string",
                      help="second coverage track"  )

    parser.set_defaults(database = None
                        , data = "coverage"
                        , taxa_level = "phylum"
                        , assembler = "idba"
                        , track_assembly = None
                        , track_one = None
                        , track_two = None)

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if not options.track_assembly:
        raise ValueError("must specify assembly track")
    if not options.track_one:
        raise ValueError("must specify track of coverage for condition 1")
    if not options.track_two:
        raise ValueError("must specify track of coverage for condition 2")

    # general variables set here
    track_one, track_two = options.track_one.replace("-", "_"), options.track_two.replace("-", "_")
    taxa_level = options.taxa_level
    assembler = options.assembler
    track_contigs = options.track_assembly.replace("-", "_")

    # connect to database
    cc = connect(options.database)
    if options.data == "coverage":

        E.info("script will output contig, length, pGC, coverage_1, coverage_2 and taxa")
        
        sys.stdout.write("contig\tlength\tgc\t%s_cov\t%s_cov\ttaxa\n" % (track_one, track_two))
        statement = """select length.scaffold_name \
                       , length.length \
                       , gc.pGC as gc \
                       , cov1.cov_mean as %(track_one)s_cov \
                       , cov2.cov_mean as %(track_two)s_cov \
                       , taxa.%(taxa_level)s \
                       from %(assembler)s_%(track_contigs)s_filtered_contigs_lengths as length \
                       , %(assembler)s_%(track_contigs)s_filtered_contigs_gc as gc \
                       , %(assembler)s_%(track_one)s_filtered_contigs_coverage_stats as cov1 \
                       , %(assembler)s_%(track_two)s_filtered_contigs_coverage_stats as cov2 \
                       , %(track_contigs)s_filtered_contigs_taxa as taxa \
                       where length.scaffold_name = gc.id \
                       and cov1.contig = length.scaffold_name \
                       and cov2.contig = length.scaffold_name \
                       and taxa.id = length.scaffold_name""" % locals()

        print statement
        for data in cc.execute(statement).fetchall():
            sys.stdout.write("\t".join(map(str, data)) + "\n")

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )



