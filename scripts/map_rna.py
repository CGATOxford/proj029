

import os, sys
import glob
import collections

datadir="/ifs/projects/proj029/data/rna"

files = glob.glob(os.path.join(datadir, "*.fastq.gz"))
suffices = {"1.fastq.gz":"fastq.1.gz", "2.fastq.gz":"fastq.2.gz"}

reps = {"WTCHG_79893_224_": "control-R1.", "WTCHG_79893_225_":"hh_il10-R1."}

for fastq in files:
    old_name = os.path.basename(fastq[:len(fastq)-10])
    suffix = fastq[len(fastq)-10:]
    new_name = "-".join(["gut", reps[old_name]]) + suffices[suffix]
    os.symlink(fastq, new_name)

    
    
