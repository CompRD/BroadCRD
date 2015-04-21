#!/usr/bin/env python
###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################


import re
import sys
import argparse
import glob


def find_latest_log(search_dir):
    file_list = sorted(glob.glob(search_dir + "/log.*"))
    if (file_list in [None, []]):
        return ""
    else:
        return file_list[-1]

# Parse the LongProto log file    
def parse_log(log_file):
    if (log_file == ""):
        return ("missing", "")
    try:
        with open(log_file,'r') as log_fd:
            log = map( str.strip, log_fd.readlines() )
            last_log_line = ""
            for line in log :
                if (str.find(line, "Samtools failed") != -1) :
                    return ("ext_error", line)
                elif (str.find(line, "Cannot allocate memory") != -1) :
                    return ("ext_error", line)
                elif (str.find(line, "not found (required by LongProto") != -1) :
                    return ("ext_error", line)
                elif (str.find(line, ": done, time used = ") != -1) :
                    return ("ok", line)
                elif (str.find(line, "Killed.  Stopping.") != -1) :
                    if (str.find(last_log_line, "loading reads from picard cache") != -1) :
                        return ("ext_error", last_log_line)
                    else:
                        return ("killed", last_log_line)
                elif (len(line) != 0 ) :
                    last_log_line = line
            if (str.find(last_log_line, "loading reads from picard cache") != -1) :
                return ("ext_error", last_log_line)
            elif (str.find(last_log_line, "There are no reads, nothing to do.") != -1):
                return ("ok", line)
            else:
                return ("int_error", last_log_line)
    except IOError:
        return ("missing", "")

# Parse the qgrind job status file
def parse_status(status_file):
    try:
        with open(status_file,'r') as status_fd:
            status = map( str.strip, status_fd.readlines() )
            if (len(status) == 0):
                return "empty"
            else:
                return status[0]
    except IOError:
        return "missing"

# Parse the qgrind job info file
def parse_info(info_file):
    try:
        with open(info_file,'r') as status_fd:
            info = map( str.strip, status_fd.readlines() )
            if (len(info) == 0):
                return "????"
            else:
                return info[0].partition(" ")[2]
    except IOError:
        return "????"


def main( argv = [__name__]  ):

    redo_file = ""

    parser=argparse.ArgumentParser(description='Determines the status of QGrind jobs by Examining the qgrind status files and LongProto logs')
    parser.add_argument('--bad', help='Only display problem jobs', action='store_true')
    parser.add_argument('--longproto', help='Only examine LongProto logs - ignore QGrind status', action='store_true')
    parser.add_argument('--timeouts', help='Include timeout regions in new job list', action='store_true')
    parser.add_argument('--redo', help='Create new job list file containing failed jobs.', metavar='FILENAME', default=None)
    parser.add_argument('--jobfile', help='job list file. (default: jobfile.txt)', metavar='FILENAME', default="jobfile.txt")
    parser.add_argument('--path', help='Partial path to working directory: dataset/attempt', default=None)
    args=parser.parse_args(argv[1:])

    work_dir = "/wga/scr4/human_assemblies"
    if (args.path != None):
        work_dir = work_dir + "/" + args.path;

    # Obtain most recent active job number
    with open("nextjob.txt",'r') as fd:
        job_count = int(map( str.strip, fd.readlines() )[0])

    # Load in list of jobs
    with open(args.jobfile,'r') as fd:
        jobs = map( str.strip, fd.readlines() )
    
    if (args.longproto):
        job_count = len(jobs)

    # List of jobs to redo
    redo_jobs = []

    job_id = 0
    for job in jobs:

        job_cmd = re.split(" ", job)
        chr_number, temp, base_range =job_cmd[2].partition(":")

        region_dir = work_dir + "/" + chr_number + "/" + base_range
        status_file = "tmp/job." + str(job_id) +"/status"
        info_file = "tmp/job." + str(job_id) +"/info"

        log_file = find_latest_log(region_dir)

        # Parse the Qgrind status, info  and LongProto logs
        log_status, last_log_line = parse_log(log_file)
        qgrind_status = parse_status(status_file) if (args.longproto == False) else "skip"
        host = parse_info(info_file)

        if (qgrind_status == "missing") :
            if (log_status != "missing"):
                print "!MISSING:",
            redo_jobs.append(job)
        elif (qgrind_status == "empty") :
            if (log_status != "missing"):
                print "!EMPTY:",
            redo_jobs.append(job)
        elif (qgrind_status == "1") :
            if (log_status == "ok"):
                print "!ERROR:",
            redo_jobs.append(job)
        elif (qgrind_status == "-999") :
            if (log_status == "ok"):
                print "!TIMEOUT:",
            elif (log_status in ["int_error", "killed"]):
                log_status = "timeout"
            if (log_status in ["ext_error", "missing"]):
                redo_jobs.append(job)
            elif (args.timeouts):
                redo_jobs.append(job)
        elif (qgrind_status == "-888") :
            if (log_status == "ok"):
                print "!EXCEPTION:",
            redo_jobs.append(job)
        elif (qgrind_status == "0") :
            if (log_status != "ok"):
                print "!OK:",
        elif (qgrind_status != "skip"):
            print "!UNKNOWN (", qgrind_status, ")",

        if ( (qgrind_status != "0") | ( (args.bad == False) | (log_status != "ok") ) ):
            print job_id, log_status, chr_number, base_range, host, last_log_line

        job_id = job_id + 1
            
        # Halt if reach the last active job
        if (job_id == job_count):
            break

    # Write redo list if requested
    if (args.redo != None):
        print len(redo_jobs), "jobs to redo written to:", args.redo 
        with open(args.redo, 'w') as redo_fd:
            for job in redo_jobs:
                redo_fd.write("%s\n" % job)
                
if __name__ == "__main__":
    sys.exit(main(sys.argv))
