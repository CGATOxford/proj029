################################################
# link to genomes from ifs/mirror/ncbi/bacteria
################################################

import os, sys
import glob

dir = "/ifs/mirror/ncbi/bacteria"
inf = open(sys.argv[1])
inf.readline()

for line in inf.readlines():
    species = line[:-1]
    print species
    print glob.glob(os.path.join(dir, species))
