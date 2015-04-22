#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Original author: Michael G. Ross <mgross@broadinstitute.org>

# A program for randomly splitting the lines of a textfile into two output
# files. Useful for constructing training and testing sets of reads.

import os
import random
import sys

try:
    source_filename = sys.argv[1]
except:
    print >>sys.stderr, ('Usage: {0} source.txt'
        .format(os.path.basename(sys.argv[0])))
    exit(1)

source_lines = 0
with open(source_filename, 'r') as source_file:
    for line in source_file:
        source_lines += 1

set_1 = random.sample(range(source_lines), source_lines / 2)
set_1.sort()

source_lines = 0
set_1_pos = 0
with open(source_filename, 'r') as source_file:
    with open(source_filename + '_1', 'w') as out_file_1:
        with open(source_filename + '_2', 'w') as out_file_2:
            for line in source_file:
                if set_1_pos < len(set_1) and source_lines == set_1[set_1_pos]:
                    out_file_1.write(line)
                    set_1_pos += 1
                else:
                    out_file_2.write(line)
                source_lines += 1

exit(0)
