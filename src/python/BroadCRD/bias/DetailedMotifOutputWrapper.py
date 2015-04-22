#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# This program wraps the output of BadCoverage's DETAILED_MOTIF_OUTPUT=True
# option to make it easier to read.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import textwrap
import sys

intervals_filename = sys.argv[1]

with open(intervals_filename, 'r') as intervals:
    for i in intervals:
        i_fields = i.split('\t')
        
        seq = i_fields[7]
        cov = i_fields[8]
        qual = i_fields[9]
        
        seq_lines = textwrap.wrap(seq, 78)
        cov_lines = textwrap.wrap(cov, 78)
        qual_lines = textwrap.wrap(qual, 78)
        
        print '-' * 80
        print ' '.join(i_fields[0:7])
        print '\n'
        for s, c, q in zip(seq_lines, cov_lines, qual_lines):
            print 's:' + s
            print 'c:' + c
            print 'q:' + q
            print '\n'
        print '\n'.join(textwrap.wrap(' '.join(i_fields[10:]), 80))
        print '-' * 80
        print '\n'
            
exit(0)
