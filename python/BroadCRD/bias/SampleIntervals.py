#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Randomly select a user-specified number of intervals at each order of
# magnitude.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import math
import random
import re
import sys

intervals_file = sys.argv[1]
num_samples = int(sys.argv[2])
logsize = int(sys.argv[3])

num_intervals = 0
tot_lengths = 0

with open(intervals_file, 'r') as intervals:
    for i in intervals:
        if re.match('.*:(\d+)-(\d+)', i):
            itoks = i.split('\t')
            (contig, irange) = itoks[0].split(':')
            (start, finish) = irange.split('-')
            ilen = int(finish) - int(start)
            if math.floor(math.log10(ilen)) == logsize:
                tot_lengths += ilen
                num_intervals += 1
print >>sys.stderr, '{0} intervals totaling {1} bases'.format(num_intervals,
    tot_lengths)

sampled_intervals = 0
sampled_lengths = 0
sample_list = []
with open(intervals_file, 'r') as intervals:
    for i in intervals:
        if re.match('.*:(\d+)-(\d+)', i):
            itoks = i.split('\t')
            (contig, irange) = itoks[0].split(':')
            (start, finish) = irange.split('-')
            ilen = int(finish) - int(start)
            if math.floor(math.log10(ilen)) == logsize:
                sample_list.append(i)

for r in random.sample(xrange(len(sample_list)), num_samples):
    print sample_list[r],

exit(0)
