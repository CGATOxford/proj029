'''
lcakegg2counts.py
==================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Count the number of alignments assigned to KEGG orthology groups. The alignments
can be counted at different levels of the KO hierarchy e.g. category A-D => high level
pathway - KO gene.

KO = Kegg orthology.

Usage
-----

Example::

   cat in.lcakegg | python lcakegg2counts.py --kegg-table=kegg.table 

Type::

   python lcakegg2counts.py --help

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

########################################
########################################
########################################
class KeggAssignment(object):
    '''
    class mapping alignment (read) with
    kegg assignment
    '''
    def __init__(self):
        '''
        initialise class
        '''
        self.identifier = None
        self.ko = []

    def getKo(self, identifier, ko):
        '''
        return data
        '''
        self.identifier = identifier
        self.ko = ko
        return self

#######################################
def lcakegg_iterator(infile):
    '''
    return iterator for kegg
    output from lcamapper.sh
    '''
    for line in infile.readlines():
        data = line[:-1].split(";")
        identifier = data[0]
        ko = re.findall("K\S+", data[2])
        
        # remove trailing ":"
        ko = [x[:-1] for x in ko]
        yield KeggAssignment().getKo(identifier, ko)

#######################################
#######################################
#######################################
class KeggTableEntry(object):
    '''
    class for assigning pathways to KO
    '''
    def __init__(self):

        '''
        initialise class. cata etc refers
        to categories in the kegg KO
        hierarchy - category A, B, C
        '''
        self.ko = None
        self.cata = None
        self.catb = None
        self.catc = None

    def getMapping(self, ko, cata, catb, catc):
        '''
        map KO to pathways
        '''
        self.ko, self.cata, self.catb, self.catc = ko, cata, catb, catc 
        return self

#######################################
def kegg_table_iterator(infile):
    '''
    iterator for kegg table. Kegg table is
    output from the keggtre2table.py script
    '''
    for line in infile.readlines():
        data = line[:-1].split("\t")
        ko, cata, catb, catc = data[-1], data[0], data[1], data[2]
        yield KeggTableEntry().getMapping(ko, cata, catb, catc)

#######################################
#######################################
#######################################
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-k", "--kegg-table", dest="kegg_table", type="string",
                      help="supply table of pathways to KO mappings. Usually by \
                            keggtre2table.py ")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("count", "proportion"), help="method for computing counts")

    parser.add_option("-l", "--level", dest="level", type="choice",
                      choices=("A", "B", "C", "D"), help="what level of the KEGG KO hierarchy \
                                                          to perform the analysis?")


    # set default options
    parser.set_defaults(method="proportion"
                        , level="A")

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if not options.kegg_table:
        raise ValueError("must specify a KEGG table")

    E.info("reading KEGG table")
    # container mapping pathway/gene to KO
    kegg = {}

    # iterate over kegg table and add kos into list
    # at the desired level
    c = E.Counter()
    for ko in kegg_table_iterator(IOTools.openFile(options.kegg_table)):        
        c.ko_input += 1
        kegg[ko.ko] = (ko.cata, ko.catb, ko.catc)
    E.info("read in KEGG table")
    
    E.info("iterating over alignments")
    result = collections.defaultdict(int)
    for alignment in lcakegg_iterator(options.stdin):
        c.input += 1
        if not alignment.ko:
            c.unmapped += 1
        else:
            c.mapped += 1
            if options.level == "A":
                catindex = 0
            elif options.level == "B":
                catindex = 1
            elif options.level == "C":
                catindex = 2
            elif options.level == "D":
                catindex = None
         
            if not catindex:
                if len(alignment.ko) > 1:
                    # doesn't matter if it is mapped
                    # to a pathway
                    for k in alignment.ko:
                            result[k] += 1
                else:
                    result[alignment.ko[0]] += 1

            else:
                # if the read maps to multiple kegg entries
                # add count for each one
                if len(alignment.ko) > 1:
                    try:
                        for k in alignment.ko:
                            result[kegg[k][catindex]] += 1
                    # if the KO is not mapped to a pathway
                    # count as no hit
                    except KeyError:
                        c.no_hit += 1
                else:
                    try:
                        result[kegg[alignment.ko[0]][catindex]] += 1
                    except KeyError:
                        c.no_hit += 1
                    
    total_mapped = c.mapped
    for pathway, count in result.iteritems():
        if options.method == "proportion":
            options.stdout.write("\t".join([pathway, str(float(count)/total_mapped)]) + "\n")
        else:
            options.stdout.write("\t".join([pathway, str(count)]) + "\n")

    options.stdout.write("no_pathway_hit\t%s" % str(float(c.no_hit)/ c.mapped))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main(sys.argv) )
    

