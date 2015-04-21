#!/usr/bin/env python3.1

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Original author: Michael G. Ross <mgross@broadinstitute.org>

# This program converts an SFF file to a FASTQ. If a barcode is specified,
# only the reads matching the barcode will be written to the FASTQ, and
# the barcode will be removed. Allows the specification of an argument to
# MakeSffInfo - which is useful when the default clip in the SFF file is
# broken and needs to be over-ridden.

import os
import re
import subprocess
import sys

# process arguments
if len(sys.argv) < 2:
    print('Usage: ' + os.path.basename(sys.argv[0])
          + ' PROJECT_NAME|SFF_FILE [FASTQ_NAME] [BARCODE] [MAKESFFINFO_ARG]',
          file=sys.stderr)
    exit(1)

if sys.argv[1][-4:] == '.sff':
    proj_name = sys.argv[1][:-4]
    sff_file = sys.argv[1]
else:    
    proj_name = sys.argv[1]
    lagoon_root = '/seq/lagoon_results/analysis/output/Broad/'
    proj_dirs_list = os.listdir(lagoon_root)
    proj_dir = []

    for d in proj_dirs_list:
        if re.match('^' + proj_name, d):
            proj_dir.append(d)
    
    if len(proj_dir) == 0:
        print('failed to match project in ' + lagoon_root, file=sys.stderr)
        exit(1)
    elif len(proj_dir) > 1:
        print('more than one match, specify project more specifically',
              file=sys.stderr)
        for d in proj_dir:
            print(d, file=sys.stderr)
        exit(1)

    proj_name = proj_dir[0]
    sff_file = os.path.join(lagoon_root, proj_name, 'rawlib.sff')

sff_file = os.path.abspath(sff_file)
sff_file_link = os.path.basename(sff_file)
sff_info_head = sff_file_link[0:-4]
sff_fasta = sff_info_head + '.fasta'

if len(sys.argv) > 2:
    fastq_name = sys.argv[2]
else:
    fastq_name = sff_info_head + '.fastq'
if len(sys.argv) > 3:
    sample_barcode = sys.argv[3]
else:
    sample_barcode = ''
if len(sys.argv) > 4:
    makesffinfo_arg = sys.argv[4]
else:
    makesffinfo_arg = None

# create directory to hold all the output
os.mkdir(proj_name)
os.chdir(proj_name)

logfilename = 'command.log'
# make a soft link to the original SFF
subprocess.call(['ln', '-s', sff_file, sff_file_link])

# run MakeSffInfo
makesffinfo_cmd = ['MakeSffInfo', 'SFF=' + sff_file_link, 'FASTA=' + sff_fasta]
if makesffinfo_arg:
    makesffinfo_cmd.append(makesffinfo_arg)
log_file = open(logfilename, 'w')
msi_ret = subprocess.call(makesffinfo_cmd, stdout=log_file, stderr=log_file)

# try to remove the too-short read
if msi_ret != 0:
    print('attempting to remove too-short reads', file=log_file)
    log_file.close()

    # extract the names of the short reads
    log_file = open(logfilename, 'r')
    for ts_line in log_file:
        tsmatch = re.match('^TOO SHORT READS:(.+)', ts_line)
        if tsmatch:
            break
    log_file.close()
    if not tsmatch:
        print('error in MakeSffInfo, apparently not a too-short read'
              ' see command.log for more info', file=sys.stderr)
        exit(1)

    # write the exclusion file for sfffile
    short_reads = tsmatch.group(1).split(',')
    short_filename = sff_info_head + '_too_short.txt'
    short_file = open(short_filename, 'w')
    for s in short_reads:
        print(s, file=short_file)
    short_file.close()

    # rename the original SFF
    sff_file_orig = sff_file_link + '.orig'
    log_file = open(logfilename, 'a')
    print('renaming ' + sff_file_link + ' to ' + sff_file_orig, file=log_file)
    os.rename(sff_file_link, sff_file_orig)

    # use sfffile to make a new SFF without the troublesome short reads
    sfffile_cmd = ['sfffile', '-e', short_filename, '-o', sff_file_link,
                   sff_file_orig]
    print('executing ' + ' '.join(sfffile_cmd), file=log_file)
    if (subprocess.call(sfffile_cmd, stdout=log_file, stderr=log_file) != 0):
        print('sfffile failed to remove the too-short reads, see '
              + logfilename, file=sys.stderr)
        exit(1)

    # re-run MakeSffInfo
    print('too-short reads stripped, retrying MakeSffInfo...', file=log_file)
    if subprocess.call(makesffinfo_cmd, stdout=log_file, stderr=log_file) != 0:
        print('MakeSffInfo failed, see ' + logfilename, file=sys.stderr)

# translate SffInfo files into a FASTQ
if (subprocess.call(['SffInfo2Fastq', 'SFFINFO=' + sff_info_head,
                     'FASTQ=' + fastq_name,
                     'BARCODE=' + sample_barcode], stdout=log_file,
                    stderr=log_file) != 0):
    print('SffInfo2Fastq failed, see ' + logfilename, file=sys.stderr)
    exit(1)

subprocess.call(['cat', logfilename])

exit(0)
