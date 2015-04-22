#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program for sampling the output of SequenceFilterCoverageBins.py to
# produce 100 bins at 0-9% relative coverage, 100 bins between 10-19%
# relative coverage, up until some reasonable 10N-(10N+9)% level,
# and ending with 100 bins between 10(N+1)-infinity% coverage.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import math
import random
import sys

bins_file = sys.argv[1]
filt_file = sys.argv[2]

cov_levels = {}
total_bins = 0
with open(bins_file, 'r') as bins:
    for bin in bins:
        cov = float(bin.split('\t')[2].rstrip()[0:-1])
        level = math.floor(cov / 10) * 10
        if level not in cov_levels:
            cov_levels[level] = []
        cov_levels[level].append(bin)
        total_bins += 1

cls = cov_levels.keys()
cls.sort()
total_lens = 0
break_point = -1
for c in cls:
    total_lens += len(cov_levels[c])
    cum_amt = float(total_lens) / float(total_bins)
    print str(c) + ' ' + str(len(cov_levels[c])) + ' ' + str(cum_amt)
    if break_point == -1 and cum_amt >= 0.99:
        break_point = c

with open(filt_file, 'w') as ff:
    remainder_bucket = []
    for c in cls:
        if c > break_point:
            remainder_bucket = remainder_bucket + cov_levels[c]
            del cov_levels[c]
        else:
            num_bins = len(cov_levels[c])
            kept_recs = random.sample(cov_levels[c], min(100, num_bins))
            kept_recs.sort()
            print '{0} {1}'.format(str(c), len(kept_recs))
            for k in kept_recs:
                print >>ff, k,
    num_bins = len(remainder_bucket)
    kept_recs = random.sample(remainder_bucket, min(100, num_bins))
    kept_recs.sort()
    print 'remainder {0}'.format(len(kept_recs))
    for k in kept_recs:
        print >>ff, k,
    
exit(0)