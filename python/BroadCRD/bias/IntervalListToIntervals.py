#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2013) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Converts the .interval_list files used in the Picard exome pipeline to specify
# the coordinates of the baits and targets to a more standard format friendly
# to our tools.

import sys

interval_list_file = sys.argv[1]

with open(interval_list_file, 'r') as interval_list:
    for interval in interval_list:
        if interval[0] != '@':
            interval_elements = interval.split('\t')[0:3]
            print '{0}:{1}-{2}'.format(interval_elements[0], 
                int(interval_elements[1]) - 1, interval_elements[2])
            
exit(0)
