1##############################################
# combine reads over run over multiple lanes
##############################################

import os, sys, re
import collections
import glob

# script for linking files
scriptsdir = "/ifs/projects/proj029/src"

# first link to data in working directory
os.system("python %s/map_samples_rna.py" % scriptsdir)

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
        p1 = " ".join(p1)
        outname1 = outprefix + ".fastq.1.gz"
        if os.path.exists(outname1): continue
        statement = "zcat %(p1)s | gzip > %(outname1)s" % locals()
        os.system(statement)
        p2 = [inf for inf in files if inf.endswith(".2.gz") and inf.find(rep) != -1]
        p2.sort()
        p2 = " ".join(p2)
        outname2 = outprefix + ".fastq.2.gz"
        
        if os.path.exists(outname2): continue
        statement = "zcat %(p2)s | gzip > %(outname2)s" % locals()
        os.system(statement)

to_remove = " ".join(to_remove)
#os.system("rm -rf %s" % to_remove)
