#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program for translating the first three columns of a BED file to 
# BadCoverage's intervals format.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import sys

try:
    input = sys.argv[1]
    output = sys.argv[2]
    
    if len(sys.argv) > 3:
        coordshift = int(sys.argv[3])
    else:
        coordshift = 0
except:
    print >>sys.stderr, ('Usage: {0} input.bed output.intervals [coordshift]'
        .format(os.path.basename(sys.argv[0])))
    exit(1)
    
with open(input, 'r') as beds:
    with open(output, 'w') as intervalsfile:
        for b in beds:
            bedfields = b.rstrip().split('\t')
            print >>intervalsfile, '{0}:{1}-{2}'.format(bedfields[0],
                 int(bedfields[1]) + coordshift, int(bedfields[2]) + coordshift)

exit(0)

