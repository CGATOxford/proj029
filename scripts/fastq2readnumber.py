import CGAT.Fastq as Fastq
import sys, re, os
import CGAT.IOTools as IOTools

fastqfile = sys.argv[1]
number = sys.argv[2]

for fastq in Fastq.iterate(IOTools.openFile(fastqfile)):
    new_id = fastq.identifier + "/%s" % str(number)
    sys.stdout.write("@%s\n%s\n+\n%s\n" % (new_id, fastq.seq, fastq.quals))


