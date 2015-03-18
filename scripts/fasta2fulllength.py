################################################
# filter infile fasta sequences based on whether
# they contain conanical start and stop codons
################################################

import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import sys

start = "ATG"
stops = ["TAA", "TAG", "TGA"]

for fasta in FastaIterator.iterate(sys.stdin):
    if not fasta.sequence.startswith(start):
        continue
    if fasta.sequence[-3:] not in stops:
        continue
    sys.stdout.write(">%s\n%s\n" % (fasta.title, fasta.sequence))
