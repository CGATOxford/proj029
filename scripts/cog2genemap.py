#################################################
# build tsv file with interesting COGs to 
# gene map
#################################################

import sys

cog2gene = {"COG0783": "COG0783: Dps/Ferritin",
            "COG2837": "COG2837: Peroxidase",
            "COG0435": "COG0435: GST"}

for x, y in cog2gene.iteritems():
    sys.stdout.write("%s\t%s\n" % (x,y))
