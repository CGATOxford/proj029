import sys, os
import gzip

inf = sys.argv[1]

for line in gzip.open(inf):
    data = line[:-1].split("\t")
    sys.stdout.write("%s\t%s\t%s\t%s\t%s\n" % (data[0], data[2],data[3], data[3], data[4]))
