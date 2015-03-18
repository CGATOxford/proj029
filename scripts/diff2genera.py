'''
diff2genera.py
===============

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take differentially expressed genes from a metatranscriptomic study and associate genes
with genera 

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

import CGAT.Diamond as Diamond
import CGAT.LCA as LCA
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import pandas


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-m", "--gene2cog", dest="gene2cog", type="string",
                      help="supply gene to COG mapping file"  )

    parser.add_option("-d", "--diff", dest="diff", type="string",
                      help="supply differentially expressed genes list")

    parser.add_option("--alignment-genes", dest="alignment_genes", type="string",
                      help="supply alignment to genes")

    parser.add_option("--alignment-taxa", dest="alignment_taxa", type="string",
                      help="supply alignment to genera")

    parser.add_option("-l", "--level", dest="level", type="choice",
                      choices = ("species", "genus", "family", "order", "phylum", "class"),
                      help="what level to classify to?")

    parser.add_option("--counts", dest="counts", action="store_true",
                      help="if set, counts are output. Default is proportion")


    parser.set_defaults(counts = False,
                        level = "genus")

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # load gene set of differentially expressed genes
    E.info("reading differentially expressed COGS")
    diff = set([x[:-1] for x in IOTools.openFile(options.diff).readlines()])
    E.info("loaded COG set")
    

    # read in gene2cog map
    E.info("reading gene to COG map")
    cog2genes = collections.defaultdict(set)
    for line in IOTools.openFile(options.gene2cog).readlines():
        data = line[:-1].split("\t")
        gene, cog = data
        if cog not in diff: continue
        cog2genes[cog].add(gene)
    E.info("loaded COG to gene mapping")

    # get total cog counts to get how many are not assigned to 
    # a genus
    cog2full_counts = collections.defaultdict(int)

    # read in reads that are mapped to each cog - this may be quite memory 
    # heavy
    read2cog = {}
    E.info("reading reads associated with COGS")
    agenes = IOTools.openFile(options.alignment_genes)
    for alignment in Diamond.best_alignment_iterator(Diamond.query_iterator(Diamond.alignment_iterator(agenes))):
        cog = [x for x in cog2genes.keys() if alignment.ref in cog2genes[x]]
        if len(cog) == 1:
            cog = cog[0]
        else:
            continue
        if cog not in diff: continue
        read2cog[alignment.qid] = cog
        cog2full_counts[cog] += 1
    E.info("loaded read to COG assignment")

    # for each cog, count the number of reads that come from
    # each genera
    cog2counts = {}
    for cog in read2cog.values():
        cog2counts[cog] = collections.defaultdict(int)

    E.info("counting reads associated with taxa")
    for lca in LCA.iterate(IOTools.openFile(options.alignment_taxa)):
        if options.level == "species":
            lev = lca.species
        elif options.level == "genus":
            lev = lca.genus
        elif options.level == "family":
            lev = lca.family
        elif options.level == "order":
            lev = lca.order
        elif options.level == "class":
            lev = lca._class
        elif options.level == "phylum":
            lev = lca.phylum
        try:
            if lev == "NA": continue
            genera, cog = lev, read2cog[lca.identifier]
            cog2counts[cog][genera] += 1
        except KeyError: 
            continue

    E.info("writing results")
    options.stdout.write("cog\ttaxa\tpreads\n")
    for cog, taxa_count in cog2counts.iteritems():
        # get total counts for each cog
        total = sum(cog2counts[cog].values())
        # output unassigned
        unassigned = cog2full_counts[cog] - total
        total = cog2full_counts[cog]

        if options.counts:
            options.stdout.write("%s\t%s\t%f\n" % (cog, "unassigned", float(unassigned)))
        else:
            options.stdout.write("%s\t%s\t%f\n" % (cog, "unassigned", (float(unassigned)/total)*100))
        for taxa, count in taxa_count.iteritems():
            if options.counts:
                options.stdout.write("%s\t%s\t%f\n" % (cog, taxa, float(count)))
            else:
                options.stdout.write("%s\t%s\t%f\n" % (cog, taxa, (float(count)/total)*100))
                                     
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
