#!/usr/bin/env python2.6

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program that generates a specified number of random reads of a
# specified length.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import random
import sys

if len(sys.argv) != 3:
    print >> sys.stderr, ('Usage: ' + os.path.basename(sys.argv[0]) +
        ' num_reads read_length')
    exit(1)

num_reads = int(sys.argv[1])
read_length = int(sys.argv[2])

bases = 'ACGT'

for r in xrange(num_reads):
    print '>READ_' + str(r)
    for rl in xrange(read_length):
        sys.stdout.write(bases[random.randint(0,3)])
    sys.stdout.write('\n')

exit(0)