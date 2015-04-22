#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Sort intervals output by BadCoverage by length.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import sys

ifiles = sys.argv[1:]
for intervals_file in ifiles:
    (fname, ext) = os.path.splitext(os.path.basename(intervals_file))
    output_file = fname + '_sorted' + ext

    print 'sorting {0}'.format(intervals_file)
    
    ints_lens_list = []
    with open(intervals_file, 'r') as intervals:
        for i in intervals:
            itoks = i.split('\t')
            start_finish_match = re.match('.*:(\d+)-(\d+)', itoks[0])
            ilen = (int(start_finish_match.group(2))
                - int(start_finish_match.group(1)))
            ints_lens_list.append((i, ilen))
    
    ints_lens_list.sort(cmp=lambda x,y: cmp(y[1], x[1]))
    
    with open(output_file, 'w') as outervals:
        for i in ints_lens_list:
            print >>outervals, i[0],

exit(0)