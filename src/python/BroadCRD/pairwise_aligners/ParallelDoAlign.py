#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Original author: Michael G. Ross <mgross@broadinstitute.org>

# This script parallelizes DoAlign. It splits the first FASTA file specified
# into approximately num_processes parts, runs DoAlign on each part in parallel
# (using the second FASTA as F2), merges the resulting QLTOUT files and
# renumbers their "sequence" and "QUERY" numbers to be consistent with the
# original FASTA.

import math
import os
import re
import subprocess
import sys
from BroadCRD.util import fasta

if len(sys.argv) < 5:
    print >>sys.stderr, ('Usage: ' + sys.argv[0] + ' file1.fasta file2.fasta ' 
        'num_processes output.qltout [doalign_param1 doalign_param2 ...]')
    exit(1)

fasta_filename = sys.argv[1]
fasta_filename2 = sys.argv[2]
num_processes = int(sys.argv[3])
qltout_combo_filename = sys.argv[4]
doalign_params = ['F2=' + fasta_filename2] + sys.argv[5:]

if fasta_filename[-6:] != '.fasta' or fasta_filename2[-6:] != '.fasta':
    print >>sys.stderr, ('first two arguments should be FASTA files with names '
        'ending in ".fasta"')
    exit(1)

# chop up the FASTA into chunks
fasta_reads = fasta.readFasta(fasta_filename)
chunk_size = int(math.ceil(len(fasta_reads) / num_processes))
subfasta_filenames = []
subfasta_qltouts = []
range_starts = []
range_stops = []
for p in range(num_processes):
    range_starts.append(min(p * chunk_size, len(fasta_reads)))
    range_stops.append(min((p + 1) * chunk_size, len(fasta_reads)))
    subfasta_filenames.append(os.path.basename(fasta_filename) + '_' + str(p))
    subfasta_qltouts.append(subfasta_filenames[p] + '.qltout')

    if os.path.exists(subfasta_filenames[p]):
        print >>sys.stderr, ('It looks like you already have a file named '
              + subfasta_filenames[p] + ' you should move or delete it before '
              + 'continuing.')
        exit(1)

    if os.path.exists(subfasta_qltouts[p]):
        print >>sys.stderr, ('It looks like you already have a file named '
              + subfasta_qltouts[p] + ' you should move or delete it before '
              + 'continuing.')
        exit(1)

for p in range(num_processes):
    fasta.writeFasta(fasta_reads, subfasta_filenames[p],
                     range(range_starts[p], range_stops[p]))

# align the chunks with many parallel DoAlign processes
DOALIGN_COMMAND = ['DoAlign'] + doalign_params
doalign_processes = []
for p in range(num_processes):
    doalign_processes.append(subprocess.Popen(DOALIGN_COMMAND +
                                              ['F1=' + subfasta_filenames[p]],
                                              stdout=open(subfasta_qltouts[p],
                                                          'w'),
                                              stderr=sys.stderr))

for proc, subfasta in zip(doalign_processes, subfasta_filenames):
    if proc.wait() != 0:
        print >>sys.stderr, 'DoAlign process for ' + subfasta + ' failed'
        exit(1)
    else:
        os.remove(subfasta)


# merge together the outputs, rewriting the query indices to be consistent
# with the original FASTA
sequence_pattern = re.compile('^(sequence\s+)(\d+)(.*)', re.S)
query_pattern = re.compile('^(QUERY\s+)(\d+)(.*)', re.S)
with open(qltout_combo_filename, 'w') as qltout_combo:
    for qltout_filename, increment in zip(subfasta_qltouts, range_starts):
        with open(qltout_filename, 'r') as qltout:
            for qltout_line in qltout:
                sequence_match = re.match(sequence_pattern, qltout_line)
                query_match = re.match(query_pattern, qltout_line)
                
                if sequence_match:
                    qltout_line = (sequence_match.group(1) + 
                                   str(int(sequence_match.group(2))
                                       + increment) +
                                   sequence_match.group(3))
                if query_match:
                    qltout_line = (query_match.group(1) + 
                                   str(int(query_match.group(2))
                                       + increment) +
                                   query_match.group(3))
                print >>qltout_combo, qltout_line,
        os.remove(qltout_filename)

exit(0)
                                                            
