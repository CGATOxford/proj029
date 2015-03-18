#####################################################
# This script provides the mapping between samples
# that have come from the WTCHG and conditions
#####################################################

import glob
import sys, re, os

# mapping between ids and conditons
code2condition = {"A":("WT", "R1"),
                  "B":("WT", "R2"),
                  "C":("WT", "R3"),
                  "D":("aIL10R", "R1"),
                  "E":("aIL10R","R2"),
                  "F":("aIL10R","R3"),
                  "G":("aIL10R","R4"),
                  "H":("Hh","R1"),
                  "I":("Hh","R2"),
                  "J":("Hh","R3"),
                  "K":("Hh","R4"),
                  "L":("HhaIL10R","R1"),
                  "M":("HhaIL10R","R2"),
                  "N":("HhaIL10R","R3"),
                  "2":("WT","R4"),
                  "11":("HhaIL10R", "R4")}


wtchgnumber2code = {"201":"A",
                    "202":"B",
                    "203":"C",
                    "204":"D",
                    "205":"E",
                    "206":"F",
                    "207":"G",
                    "208":"H",
                    "209":"I",
                    "210":"J",
                    "211":"K",
                    "212":"L",
                    "213":"M",
                    "214":"N",
                    "215":"2",
                    "216":"11"}

# there are two separate projects for the data so two
# separate data directories are specified
runs = ["lane1", "lane2", "lane3", "lane4"]
datadir1="/ifs/projects/proj029/data/RNA/P130502"
datadir2="/ifs/projects/proj029/data/RNA/P140205"

# iterate over directories and assign new names
for directory in [datadir1, datadir2]:
    for run in runs:
        for infile in glob.glob(os.path.join(os.path.join(directory, run), "*fastq*gz")):
            newname = os.path.basename(infile).replace("WTCHG_", "")
            if newname.endswith("_1.fastq.gz"):
                newname = newname.replace("_1.fastq.gz", ".fastq.1.gz")
            elif newname.endswith("_2.fastq.gz"):
                newname = newname.replace("_2.fastq.gz", ".fastq.2.gz")
            data = newname.split("_")
            r, identifier = data[0], data[1].split(".")[0]
            try:
                condition = code2condition[wtchgnumber2code[identifier]][0]
                rep = code2condition[wtchgnumber2code[identifier]][1]
            except KeyError: 
                continue
            newname = "-".join([r, condition, rep]) + ".%s" % ".".join(data[1].split(".")[1:])
            if os.path.exists(os.path.abspath(newname)): continue
            os.symlink(os.path.abspath(infile), os.path.abspath(newname))
            
            

