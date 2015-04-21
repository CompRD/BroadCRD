#!/usr/bin/env python
###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################


import sys
import argparse
import fileinput 

def main( argv = [__name__]  ):

    redo_file = ""

    parser=argparse.ArgumentParser(description='Extracts jobs from a QGrind jobfile based on job id')
    parser.add_argument('jobfile', help='Original jobfile (default: stdin)', nargs = "?")
    parser.add_argument('--extract', help='File with list of jobs ids to extract (in first column)', metavar='FILENAME', default=None)
    parser.add_argument('--output', help='New jobfile (default: stdout)', metavar='FILENAME')
    args=parser.parse_args(argv[1:])



    # Load in list of jobs
    jobs = []
    for line in fileinput.input(args.jobfile):
        jobs.append(line.strip())

    # Load in list of jobs to extract
    with open(args.extract,'r') as fd:
        extract = map( str.strip, fd.readlines() )
    
    # List of jobs to redo
    redo_jobs = []

    for entry in extract:

        job_id = 0
        (left, tmp, right) = entry.partition(" ")
        if (left[0] != '!'):
            job_id = int(left)
        else :
            job_id = int(right.partition(" ")[0])

        redo_jobs.append(jobs[job_id])

    # Write redo list if requested
    if (args.output != None):
        print len(redo_jobs), "jobs to redo written to:", args.output
        with open(args.output, 'w') as redo_fd:
            for job in redo_jobs:
                redo_fd.write("%s\n" % job)
    else:
        for job in redo_jobs:
            print("%s" % job)

                
if __name__ == "__main__":
    sys.exit(main(sys.argv))
