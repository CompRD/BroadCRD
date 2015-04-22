#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program for converting one-based intervals to zero-based intervals.

# Original author: Michael G. Ross <mgross@broadinstitute.org>


import os
import re
import sys

try:
    input = sys.argv[1]
    output = sys.argv[2]
except:
    print >>sys.stderr, ('Usage: {0} input.intervals output.intervals'
        .format(os.paths.basename(sys.argv[0])))
    exit(1)

with open(input, 'r') as input_intervals:
    with open(output, 'w') as output_intervals:
        for i in input_intervals:
            (contig, remainder) = i.split(':')
            (start, stop) = remainder.split('-')
            start = str(int(start) - 1)
            stop = str(int(stop) - 1)
            print >>output_intervals, '{0}:{1}-{2}'.format(contig, start, stop)

exit(0)
            
