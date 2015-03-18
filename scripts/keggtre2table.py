'''
keggtre2table.py
=================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take the kegg.tre file that comes as part of the mtools package (Dan Huson)
and output a table of pathways 2 KEGG orthology mappings.

Usage
-----

Example::

   python keggtre2table.py  --kegg-tre=kegg.tre --map=kegg.map

Type::

   python keggtre2table.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import collections

import CGAT.Experiment as E

#####################################
#####################################
#####################################
def readMap(mapfile):
     '''
     return a dictionary mapping of KO code
     to pathway/gene ko
     '''
     code2pathway = {}
     for line in mapfile.readlines():
          data = line[:-1].split("\t")
          if len(data) == 1: continue
          code2pathway[data[0]] = data[1]
     return code2pathway

#####################################
#####################################
#####################################
def tre2table(trefile, code2pathway):
     '''
     write the mapping of KO ids
     to pathways
     '''

     # iterate over the tree file using the 
     # trainling KO id as the key for 
     # subsequent iteration over lower level 
     # trees. The order as defined by KEGG
     # categories is A, B, C, D - hence atree,
     for alldata in open(trefile).readlines():

          # seperate levels of the tree are 
          # separated by "(" e.g. "(((" = category A
          alldata = alldata.split("(((")
          for atree in alldata:
               atree = atree.split("((")

               # for this and subsequent trees branches in the
               # tree the trailing KO is the identifier for 
               # lower level branches
               a = atree[-1].split(")")[-1][:-1]

               # The consequnces of removing the -3 A
               # pathways has not been extensively investigated
               # although they are likely to be no hits or
               # uncharacterised into KEGG pathways
               if a == '' or a == "-3;": continue
               atree = [re.sub("\([0-9]{7,},", ",", x) for x in atree]

               bdict = collections.defaultdict(list)
               for btree in atree:
                    btree = btree.split("(")
                    # remove trailing spaces at start of
                    # btree list
                    if '' in btree: 
                         btree = btree[1:]
                    b = btree[-1].split(")")[-1][:-1]

                    # need to remove the last KO entry
                    btree = [re.sub("\)20000.*,", ",", x) for x in btree]
                    for ctree in btree:
                         ctree = ctree.split(")")
                         c = ctree[-1][:-1]
                         for dtree in ctree:
                              dtree = dtree.split(",")
                              if '' in dtree: continue
                              for d in dtree:
                                   try:
                                        sys.stdout.write("\t".join([code2pathway[a], 
                                                                    code2pathway[b], 
                                                                    code2pathway[c], 
                                                                    code2pathway[d]]) + "\n")
                                   except KeyError:

                                        # where multiple C KOs map to the 
                                        # same set of genes in D they are 
                                        # assumed to be differenct pathways
                                        # containing the same set of genes
                                        cs = c.split(",")
                                        for c in cs:
                                             sys.stdout.write("\t".join([code2pathway[a], 
                                                                         code2pathway[b], 
                                                                         code2pathway[c], 
                                                                         code2pathway[d]]) + "\n")


#####################################
#####################################
#####################################
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )
    parser.add_option("-k", "--kegg-tre", dest="kegg_tre", type="string",
                      help="specify kegg.tre file (from mtools)")
    parser.add_option("-m", "--map", dest="map", type="string",
                      help="specify map file between KO code and pathway/gene name")

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    E.info("reading map file")
    code2pathway = readMap(open(options.map))
    E.info("writing map")
    tre2table(options.kegg_tre, code2pathway)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )












