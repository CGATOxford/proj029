###############################################
# mapping the cgaes and mothers to samples
###############################################

import glob
import sys, re, os

# conditions to cages mapping
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


code2cage = {"2": "C1",
             "A": "C1",
             "F": "C1",
             "G": "C1",
             "H": "C2",
             "I": "C2",
             "11":"C2",
             "L": "C2",
             "B": "C3",
             "C": "C3",
             "D": "C3",
             "E": "C3",
             "J": "C4",
             "K": "C4",
             "M": "C4",
             "N": "C4"}

code2mother = {"2": "M1",
               "A": "M1",
               "F": "M1",
               "G": "M1",
               "H": "M1",
               "I": "M1",
               "11":"M1",
               "L": "M1",
               "B": "M1",
               "C": "M2",
               "D": "M2",
               "E": "M2",
               "J": "M2",
               "K": "M2",
               "M": "M2",
               "N": "M2"}

sys.stdout.write("track\tcage\tdam\n")
for code, condition in code2condition.iteritems():
    condition = "stool_" + condition[0] + "_%s" % condition[1] 
    cage = code2cage[code]
    mother = code2mother[code]
    sys.stdout.write("\t".join( [condition, cage, mother]) + "\n")





