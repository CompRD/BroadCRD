#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program for translating BadCoverage's intervals format into the BED
# format required by liftOver. The zerobase argument is set to true if
# the intervals are based on zero-based indexing - it will shift them to
# produce a one-based indexed BED file.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import sys

try:
    input = sys.argv[1]

    if len(sys.argv) > 2:
        output = sys.argv[2]
    elif re.match('.*.intervals$', input):
        output = re.sub('.intervals$', '.bed', input)
    else:
        output = input + '.bed'
    
    if len(sys.argv) > 3:
        shiftcoords = int(sys.argv[3])
    else:
        shiftcoords = 0
except:
    print >>sys.stderr, ('Usage: {0} input.intervals output.bed '   
        '[shift_coordinates_offset]'.format(os.path.basename(sys.argv[0])))
    exit(1)

with open(input, 'r') as intervals:
    with open(output, 'w') as bedfile:
        print >>bedfile, 'track name=' + input
        for i in intervals:
            if re.match('^.*:\d+-\d+\s*', i):
                (contig, remainder) = i.split('\t')[0].rstrip().split(':')
                (start, stop) = remainder.split('-')
                start = str(int(start) + shiftcoords)
                stop = str(int(stop) + shiftcoords)
                print >>bedfile, contig + '\t' + start + '\t' + stop

exit(0)
