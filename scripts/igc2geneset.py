###################################################
# convert IGC annotations to geneset for
# enrichment analysis
###################################################

import sys, re, os
import gzip

inf = sys.stdin

for line in inf.readlines():
    data = line[:-1].split("\t")
    ontology = "cog_nog"
    evidence = "NA"
    gene_id, cog_id, description = data[1], data[8], data[12]
    if description == "Function unknown": continue
    descriptions = description.split(";")
    cogs = cog_id.split(";")
    if len(descriptions) > 1:
        descriptions = description.split(";")
        if len(cogs) > 1:
            for c in cogs:
                for d in descriptions:
                    sys.stdout.write("\t".join([ontology, gene_id, c, d, evidence]) + "\n")
        else:
            for d in descriptions:
                sys.stdout.write("\t".join([ontology, gene_id, cog_id, d, evidence]) + "\n")
    else:
        if len(cogs) > 1:
            for c in cogs:
                sys.stdout.write("\t".join([ontology, gene_id, c, description, evidence]) + "\n")
        else:
            sys.stdout.write("\t".join([ontology, gene_id, cog_id, description, evidence]) + "\n")
