############################################
############################################
############################################
# combine data across runs
############################################
############################################
############################################

import os
import sys
import glob
import CGATPipelines.Pipeline as P

from ruffus import *

P.getParameters(
    ["pipeline.ini"])


PARAMS = P.PARAMS


sample_map = {"201": "stool-day0-R1",
              "202": "stool-day0-R2",
              "203": "stool-day0-R3",
              "204": "stool-day0-R4",
              "205": "stool-day0-R5",
              "206": "stool-day0-R6",
              "207": "stool-day0-R7",
              "208": "stool-day0-R8",
              "209": "stool-day3-R1",
              "210": "stool-day3-R2",
              "211": "stool-day3-R3",
              "212": "stool-day3-R4",
              "213": "stool-day3-R5",
              "214": "stool-day3-R6",
              "215": "stool-day3-R7",
              "216": "stool-day3-R8",
              "217": "stool-day6-R1",
              "218": "stool-day6-R2",
              "219": "stool-day6-R3",
              "220": "stool-day6-R4",
              "221": "stool-day6-R5",
              "222": "stool-day6-R6",
              "223": "stool-day6-R7",
              "224": "stool-day6-R8",
              "225": "stool-day14-R1",
              "226": "stool-day14-R2",
              "227": "stool-day14-R3",
              "228": "stool-day14-R4",
              "229": "stool-day14-R5",
              "230": "stool-day14-R6",
              "231": "stool-day14-R7",
              "232": "stool-day14-R8",
              "233": "stool-day28-R1",
              "234": "stool-day28-R2",
              "235": "stool-day28-R3",
              "236": "stool-day28-R4",
              "237": "stool-day28-R5",
              "238": "stool-day28-R6",
              "239": "stool-day28-R7",
              "240": "stool-day28-R8"} 

basedir="/ifs/projects/proj029/backup/timecourse/RNA"
dirs=["run1", "run2", "run3", "run4"]

def combine(sample_map, read_in_pair=1):
    '''
    combine files per rad pair across lanes
    '''
    all_infiles = []
    outfiles = []
    for index in sample_map.keys():
        found = []
        if read_in_pair==1:
            outf = sample_map[index] + ".fastq.1.gz"
            suffix = "_1.fastq.gz"
        else:
            outf = sample_map[index] + ".fastq.2.gz"
            suffix = "_2.fastq.gz"
        for directory in dirs:
            direc = os.path.join(basedir, directory)
            pattern = os.path.join(direc, "*" + index + suffix)
            f = glob.glob(pattern)
            found.extend(f)
        infiles = " ".join(found)
        if os.path.exists(outf):
            continue
        all_infiles.extend(found)
    return all_infiles, outfiles

############################################
############################################
############################################

INFILES, OUTFILES = combine(sample_map, read_in_pair=1)

@collate(INFILES, 
         regex(r"/ifs/projects/proj029/backup/timecourse/RNA/run[0-9]/(.*)_([0-9]*)_([0-9]*)_1.fastq.gz"), 
         r"\3.fastq.2.gz")
def runCombineRead1(infiles, outfile):
    '''
    combine files
    '''
    infs = " ".join(list(infiles)) 
    statement = "zcat %(infs)s | gzip > %(outfile)s" % locals()
    P.run(statement)

############################################
############################################
############################################

INFILES, OUTFILES = combine(sample_map, read_in_pair=2)

@collate(INFILES, 
         regex(r"/ifs/projects/proj029/backup/timecourse/RNA/run[0-9]/(.*)_([0-9]*)_([0-9]*)_2.fastq.gz"), 
         r"\3.fastq.2.gz")
def runCombineRead2(infiles, outfile):
    '''
    combine files
    '''
    infs = " ".join(list(infiles)) 
    statement = "zcat %(infs)s | gzip > %(outfile)s" % locals()
    P.run()

############################################
############################################
############################################

# rename files

def getFiles(sample_map, read_in_pair=1):
    '''
    get input and output files
    '''
    all_files = []
    for ind, new_id in sample_map.iteritems():
        inf = ind + ".fastq." + str(read_in_pair) + ".gz"
        outf = new_id + ".fastq." + str(read_in_pair) + ".gz"
        if os.path.exists(outf):
            continue
        in_and_out = [inf, outf]
        all_files.append(in_and_out)
    return all_files

############################################
############################################
############################################

@follows(runCombineRead2)
@files(getFiles(sample_map, read_in_pair=2))
def renameRead2(infile, outfile):
    '''
    rename the files
    '''
    statement = '''mv %(infile)s %(outfile)s''' % locals()
    P.run()
           


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))






