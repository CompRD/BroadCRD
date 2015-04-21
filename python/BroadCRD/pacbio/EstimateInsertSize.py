#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program that attempts to estimate insert sizes for PacBio wells
# by averaging the subreads in each well, excluding the first and last
# one, which might be misleading because these subreads start and stop
# at random positions between the adapters. Of course, this could
# introduce a bias towards underestimation because long reads might
# not make the required three loops. The program should also correct
# for the high insertion rate, which will make raw subread lengths
# overestimate the insert length.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import sys

import numpy

from BroadCRD.util.fasta import readFasta

try:
    subreads_fasta = sys.argv[1]
except:
    print >>sys.stderr, ('Usage: ' + os.path.basename(sys.argv[0])
        + 'subreads.fasta')
    exit(1)

# gather all the reads in each well together
names_reads = readFasta(subreads_fasta)
well_subread_starts = {}
well_subread_ends = {}
for nr in names_reads:
    read_name_match = re.match('(.*)/(\d+)_(\d+)', nr[0])
    well_name = read_name_match.group(1)
    subread_start = int(read_name_match.group(2))
    subread_end = int(read_name_match.group(3))
    
    if subread_end <= subread_start:
        print >>sys.stderr, 'subreads should not end before they begin'
        exit(1)

    if well_name not in well_subread_starts:
        well_subread_starts[well_name] = []
        well_subread_ends[well_name] = []
    elif well_subread_ends[well_name][-1] >= subread_start:
        print >>sys.stderr, 'subreads should be in read-position order'
        exit(1)
        
    well_subread_starts[well_name].append(subread_start)
    well_subread_ends[well_name].append(subread_end)

# compute the average subread size for each well, discarding the first
# and last subreads, since they could be starting or quitting in some
# random location between SMRT bell adapters

insert_estimates = []
for well_name in well_subread_starts:
    num_complete_subreads = len(well_subread_starts[well_name]) - 2
    if num_complete_subreads > 0:
        avg_length = 0
        for s in range(1, len(well_subread_starts[well_name]) - 1):
            avg_length += (well_subread_ends[well_name][s]
                - well_subread_starts[well_name][s])
        avg_length /= num_complete_subreads
        print well_name + "\t" + str(avg_length)
        insert_estimates.append(avg_length)
        

insert_mean = numpy.mean(insert_estimates)
insert_var = numpy.var(insert_estimates)
insert_stddev = numpy.std(insert_estimates)

print 'mean estimated insert size = {0}'.format(insert_mean)
print 'var estimated insert size = {0}'.format(insert_var)
print 'stddev estimated insert size = {0}'.format(insert_stddev)
print 'number of wells = {0}'.format(len(insert_estimates))

exit(0)
