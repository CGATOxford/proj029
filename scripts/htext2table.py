'''
htext2table.py
==============

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Parse Kegg orthology data from an htext file into readable and
useful format. The resulting format is a mapping of each kegg
pathway with the genes within that pathway.

Usage
-----

Example::

   cat in_htext.keg | python mtext2table.py --category=C 

   This example results in the outputting of a pathway
   to gene map of KEGG category C pathways.

Type::

   python mtext2table.py --help

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

def getData(inf):
    '''
    parse KO's from infile
    '''
    # first 8 lines are not useful
    for i in range(8):
        inf.readline()
    
    # read data into memory 
    # This is not memory efficient
    # but I am struggling!
    all_data = []
    split_data = []
    for line in inf.readlines():
        line = line.replace(" ", "_")
        data = line[:-1]
        if data == "#":
            all_data.append(split_data)
            split_data = []
        else:
            split_data.append(data)

    return all_data

#################################
#################################
#################################
def c2ko(all_data):
    '''
    return mapping of C category to 
    KOs, genes and descriptions.
    '''
    result = {}
    all_indices = []
    kos = set()
    for data in all_data:
        C_indices = []
        c = 0
        for annotation in data:
            c += 1
            if annotation.startswith("C_"):
                if "PATH" not in annotation: continue
                C_indices.append(c - 1)
        stop = 0
        for index in C_indices:
            stop += 1
            C = data[index]
            C = [x for x in C.split("_") if x != '']
            C = "_".join(C[2:-1])
            genes = []
            kos = []
            descs = []
            try:
                D = data[index + 1:C_indices[stop]]
            except IndexError:
                continue
            for description in D:
                if description.startswith("B"): continue
                d = [x for x in description.split("_") if x != '']
                ko = d[1]
                kos.append(ko)
                try:
                    representative_gene = [x[:-1] for x in d if x.endswith(";")][0]
                except IndexError:
                    representative_gene = ko
                desc = description.split("_[EC")[0]
                desc = " ".join([x for x in desc.split("_") if x != '']).split(";")
                if len(desc) == 1:
                    desc = desc[0]
                else:
                    desc = desc[1]
                genes.append(representative_gene)
                descs.append(desc)
            result[C] = [kos, genes, descs]
    return result

#################################
#################################
#################################
def b2ko(all_data):
    '''
    return mapping of B category to 
    KOs, genes and descriptions.
    '''
    result = {}
    all_indices = []
    for data in all_data:
        B_indices = []
        c = 0
        for annotation in data:
            c += 1
            if annotation == "B": continue 
            if annotation.startswith("B_"):
                B_indices.append(c - 1)
        stop = 0
        for index in B_indices:
            stop += 1
            B = data[index]
            B = re.sub(r"\W", "", B).split("__")[1][1:-1]
            genes = []
            kos = []
            descs = []
            try:
                D = data[index + 1:B_indices[stop]]
            except IndexError:
                continue
            for description in D:
                if description == "B" or description.startswith("C_"): continue
                d = [x for x in description.split("_") if x != '']
                ko = d[1]
                kos.append(ko)
                try:
                    representative_gene = [x[:-1] for x in d if x.endswith(";")][0]
                except IndexError:
                    representative_gene = ko
                desc = description.split("_[EC")[0]
                desc = " ".join([x for x in desc.split("_") if x != '']).split(";")
                if len(desc) == 1:
                    desc = desc[0]
                else:
                    desc = desc[1]
                genes.append(representative_gene)
                descs.append(desc)
            result[B] = [kos, genes, descs]
    return result

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

    parser.add_option("--category", dest="category", type="choice",
                      choices = ("B", "C"), help="supply help"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    data = getData(options.stdin)
    if options.category == "B":
        options.stdout.write("Category B pathway\tKO\tGenes\tDescriptions\n")
        for pathway, descriptions in b2ko(data).iteritems():
            options.stdout.write("\t".join([pathway, "; ".join(descriptions[0]), "; ".join(descriptions[1]), "; ".join(descriptions[2])]) + "\n")

    elif options.category == "C":
        options.stdout.write("Category C pathway\tKO\tGenes\tDescriptions\n")
        for pathway, descriptions in c2ko(data).iteritems():
            options.stdout.write("\t".join([pathway, "; ".join(descriptions[0]), "; ".join(descriptions[1]), "; ".join(descriptions[2])]) + "\n")
    else:
        raise ValueError("must specify the category of pathway")


    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


    

                
