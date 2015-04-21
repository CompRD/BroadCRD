#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A quick and dirty program for measuring the approximate distance between
# independently aligned Ion jumps, with some equally quick and dirty filtering
# to report some quick and dirty statistics.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import math
import re
import sys
import numpy

from BroadCRD.util.samtools import openbam

sam1 = sys.argv[1]
sam2 = sys.argv[2]

read1start = {}
read1contig = {}
read1len = {}
read1names = {}

with openbam(sam1) as sam1reads:
    for read in sam1reads:
        if read[0] != '@':
            read_fields = read.split('\t')
            read_id = re.match('^.*:(.*:.*)', read_fields[0]).group(1)
            if not (int(read_fields[1]) & 0x4):
                read1start[read_id] = int(read_fields[3])
                read1contig[read_id] = read_fields[2]
                read1len[read_id] = len(read_fields[9])
                read1names[read_id] = read_fields[0]

dists = []
one_mate_unaligned = 0
wrong_contig = 0
with openbam(sam2) as sam2reads:
    for read in sam2reads:
        if read[0] != '@':
            read_fields = read.split('\t')
            read_id = re.match('^.*:(.*:.*)', read_fields[0]).group(1)
            if not (int(read_fields[1]) & 0x4):
                if read_id not in read1start:
                    one_mate_unaligned += 1
                elif read1contig[read_id] != read_fields[2]:
                    wrong_contig += 1
                else:
                    if read1start[read_id] < int(read_fields[3]):
                        d = (int(read_fields[3]) + len(read_fields[9])) \
                            - read1start[read_id]
                    else:
                        d = (read1start[read_id] + read1len[read_id]) \
                            - int(read_fields[3])
                    print '{0} {1}:{2} | {3} {4}:{5}'.format \
                        (read1names[read_id], read1contig[read_id],
                        read1start[read_id], read_fields[0], read_fields[2],
                        read_fields[3])
                    if d < 100000:
                        dists.append(d)
                    else:
                        print 'discarding distance of {0}'.format(d)

print 'correctly paired: {0}'.format(len(dists))
print 'mean separation: {0} (stddev: {1})'.format(numpy.mean(dists),
    numpy.std(dists))
print dists
print 'unmapped mates: {1} {0}'.format(one_mate_unaligned, len(read1start) - len(dists))
print 'wrong contig: {0}'.format(wrong_contig)
                

exit(0)