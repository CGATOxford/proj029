# build design file

import sys

mapping = {"stool-Hh-R1": "Hh",
           "stool-Hh-R2": "Hh",
           "stool-Hh-R3": "Hh",
           "stool-Hh-R4": "Hh",
           "stool-HhaIL10R-R1": "HhaIL10R",
           "stool-HhaIL10R-R2": "HhaIL10R",
           "stool-HhaIL10R-R3": "HhaIL10R",
           "stool-HhaIL10R-R2": "HhaIL10R",
           "stool-WT-R1": "WT",
           "stool-WT-R2": "WT",
           "stool-WT-R3": "WT",
           "stool-WT-R4": "WT",
           "stool-aIL10R-R1": "aIL10R",
           "stool-aIL10R-R2": "aIL10R",
           "stool-aIL10R-R3": "aIL10R",
           "stool-aIL10R-R4": "aIL10R"}

sys.stdout.write("track\tinclude\tgroup\tpair\n")
for x, y in mapping.iteritems():
    sys.stdout.write("\t".join([x, "1", y, "0"]) + "\n")
