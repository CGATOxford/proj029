#########################################
# build design file for mouse expression
# data
#########################################

import glob
import CGAT.Pipeline as P

bams = glob.glob("*.bam")

outf = open("design.tsv", "w")
outf.write("track\tinclude\tgroup\tpair\n")

for bam in bams:
    track = P.snip(bam, ".bam")
    if "HhaIL10R" in track:
        group = "HhaIL10R"
    else:
        group = "WT"
    outf.write("%(track)s\t1\t%(group)s\t0\n" % locals())
outf.close()

