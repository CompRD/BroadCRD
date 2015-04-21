#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Select output of BadCoverage's DETAILED_MOTIF_OUTPUT=True mode that is
# "bad" - i.e. relative coverage < 0.1.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import sys

intervals_file = sys.argv[1]

with open(intervals_file, 'r') as intervals_list:
    for interval in intervals_list:
        fields = interval.rstrip().split('\t')
        if len(fields) > 1:
            try:
                if float(fields[1]) < 0.1:
                    print interval,
            except:
                print interval,
    
exit(0)
