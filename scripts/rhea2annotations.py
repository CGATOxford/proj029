'''
rhea2annotions
====================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

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
import collections
import CGAT.IOTools as IOTools

import CGAT.Experiment as E


def getRheaIdFromKeggIds(kegg_file, kegg_id):
    '''
    get all Rhea ids associated with user-specified
    kegg ids
    '''
    for line in kegg_file:
        data = line[:-1].split("\t")
        try:
            rhea_id, kid = data[0], data[3]
            if kid == kegg_id:
                r = rhea_id
                break
        except IndexError:
            continue
    return r

def getGenomesFromMicrome(microme_file, rhea_id):
    '''
    at the moment are just getting annotations
    at the genus level
    '''
    genomes = set()
    for line in IOTools.openFile(microme_file).readlines():
        data = line[:-1].split("\t")
        genome, rid = data[3].split(" ")[0], data[7]
        if rid == rhea_id:
            genomes.add(genome)
    return genomes

def resolveAnnotations(genomeswithreaction, testgenomes):
    '''
    resolve the genomes in the testgenomes file
    with those genomes that have the reaction
    '''
    toannotate = set([line[:-1] for line in open(testgenomes)])
    with_reaction = genomeswithreaction.intersection(toannotate)
    without_reaction = toannotate.difference(with_reaction)
    for genome in list(with_reaction) + list(without_reaction):
        if genome in with_reaction:
            yield "%s\t%s\n" % (genome, "yes")
        else:
            yield "%s\t%s\n" % (genome, "no")
            
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("--microme-file", dest="microme_file", type="string",
                      help="provide microme file"  )

    parser.add_option("-k", "--kegg-id", dest="kegg_id", type="string",
                      help="provide kegg reaction ID's from which to get genomes"  )

    parser.add_option("-g", "--genomes", dest="genomes", type="string",
                      help="list of genus-level genomes to annotate"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    assert options.kegg_id, "must specify kegg id"
    assert options.genomes, "must specify genomes to annotate"
    assert options.microme_file, "must specify microme file from which to annoate genomes"

    E.info("obtaining genomes with reaction %s" % options.kegg_id)
    genomes_with_reaction = getGenomesFromMicrome(
        options.microme_file, getRheaIdFromKeggIds(options.stdin, options.kegg_id)
        )

    E.info("resolving annotations for %s" % options.genomes)
    
    for annotation in resolveAnnotations(genomes_with_reaction, options.genomes):
        options.stdout.write(annotation)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
