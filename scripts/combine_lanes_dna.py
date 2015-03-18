##############################################
# combine reads over run over multiple lanes
##############################################

import os, sys, re
import collections
import glob

# script for linking files
scriptsdir = "/ifs/projects/proj029/src"

# first link to data in working directory
os.system("python %s/map_samples_dna.py" % scriptsdir)

# to remove at the end
to_remove = glob.glob("*.fastq*")


# iterate over files and combine those that have
# the same index i.e. condition
file_list = collections.defaultdict(list)
for inf in to_remove:
    name_split = inf.split("-")
    index = name_split[1]
    file_list[index].append(inf)

reps = ["R1", "R2", "R3", "R4"]
for condition, files in file_list.iteritems():
    for rep in reps:
        outprefix = "stool-" + condition + "-%s" % rep
        p1 = [inf for inf in files if inf.endswith(".1.gz") and inf.find(rep) != -1]
        p1.sort()
        outname1 = outprefix + ".fastq.1.gz"
        statement = "zcat %s | gzip > %s&" % (" ".join(p1), outname1)
        os.system(statement)
        p2 = [inf for inf in files if inf.endswith(".2.gz") and inf.find(rep) != -1]
        p2.sort()
        outname2 = outprefix + ".fastq.2.gz"
        statement = "zcat %s | gzip > %s&" % (" ".join(p2), outname2)
        os.system(statement)

to_remove = " ".join(to_remove)
#os.system("rm -rf %s" % to_remove)
