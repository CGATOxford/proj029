################################################
################################################
################################################
# build functional annoations for COGs
################################################
################################################
################################################

import sys, re, os

functions = {}
fun = open("fun2003-2014.tab")
fun.readline()
for line in fun.readlines():
    data = line[:-1].split("\t")
    functions[data[0]] = data[1]


genes = open("all.funccat.txt")
genes.readline()
for line in genes.readlines():
    data = line[:-1].split("\t")
    cog, fun = data[0], data[1]
    if len(fun) > 1:
        for f in list(fun):
            try:
                sys.stdout.write("COG\t%s\t%s\t%s\tNA\n" % (cog, f, functions[f]))
            except KeyError:
                continue
    else:
        sys.stdout.write("COG\t%s\t%s\t%s\tNA\n" % (cog, fun, functions[fun]))

