#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A simple program for selecting the intervals output by BadCoverage by
# regular expression.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import optparse
import re
import sys

oparser = optparse.OptionParser()
oparser.add_option('--regexp', action='store', type='string')
oparser.add_option('--min_length', action='store', type='int')
(options, args) = oparser.parse_args()

intervals_file = args[0]

ints_lens_list = []
with open(intervals_file, 'r') as intervals:
    for i in intervals:
        kept = True
        if options.regexp and not re.match(regexp, i):
            kept = False
        if options.min_length:
            data_match = re.match('.*:(\d+)-(\d+)', i)
            if (int(data_match.group(2)) - int(data_match.group(1)) <
                options.min_length):
                kept = False
        if kept:
            print i,
                        
exit(0)
