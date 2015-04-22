#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program for filtering BadCoverage's BAD_BIN_OUT output, removing bins
# with non-unique or RepeatMasked content.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import sys

bins_file = sys.argv[1]
rm_filename = sys.argv[2]
filt_bins_file = sys.argv[3]

bin_sequences = {}
with open(bins_file, 'r') as bins:
    for bin in bins:
        seq = bin.split('\t')[1]
        if seq != 'CONTENT':
            if seq not in bin_sequences:
                bin_sequences[seq] = 1
            else:
                bin_sequences[seq] += 1

kept_count = 0
discarded_count = 0
with open(bins_file, 'r') as bins:
    temp_filt_fname = filt_bins_file + '_temp1'
    with open(temp_filt_fname, 'w') as out:
        for bin in bins:
            seq = bin.split('\t')[1]
            if seq != 'CONTENT':
                if seq not in bin_sequences:
                    print >>sys.stderr, 'previously unseen sequence'
                    exit(1)
                if bin_sequences[seq] == 1:
                    kept_count += 1
                    print >>out, bin,
                else:
                    discarded_count += 1
bin_sequences = {}
print '{0} bins kept, {1} bins discarded'.format(kept_count, discarded_count)

print 'reading RepeatMasker data'
with open(rm_filename, 'r') as rmout_file:
    rminfo = {}
    rmlens = {}
    for rmout_line in rmout_file:
        rmline_match = re.match('^\s*(\d+\.*\d+\s+){4}(\S+)\s+(\d+)\s+(\d+)\s+' 
                                + '\\((\d+)\\)', rmout_line)
        if rmline_match:
            rmseq = rmline_match.group(2)
            rmstart = int(rmline_match.group(3)) - 1
            rmend = int(rmline_match.group(4)) - 1
            rmremain = int(rmline_match.group(5))
            if rmseq not in rminfo:
                print 'Reading exclusions from ' + rmseq
                rminfo[rmseq] = list()
                rmlens[rmseq] = 0
            rminfo[rmseq].append((rmstart, rmend))
            if rmlens[rmseq] == 0:
                rmlens[rmseq] = rmend + rmremain + 1
            elif rmlens[rmseq] != rmend + rmremain + 1:
                print('size mismatch in ' + rmseq)     

print 'constructing RepeatMasker mask arrays'
bad_areas = {}
for rmseq, maxlen in rmlens.items():
    bad_areas[rmseq] = bytearray(maxlen)
    for (bad_start, bad_end) in rminfo[rmseq]:
        bad_areas[rmseq][bad_start:(bad_end + 1)] = \
            (1 for i in range(bad_end - bad_start + 1))
    if len(bad_areas[rmseq]) != maxlen:
        print 'whoops size mismatch: ' + rmseq
        print len(bad_areas[rmseq])
        print maxlen
        exit(1)
    else:
        tot_bases = len(bad_areas[rmseq])
        tot_bad = sum(bad_areas[rmseq])
        print (str(tot_bad) + ' out of ' + str(tot_bases)
              + ' ({0:.2f})'.format(float(tot_bad) / float(tot_bases))
              + ' base pairs are bad in ' + rmseq)

print 'removing bins marked by RepeatMasker'
kept_count = 0
discarded_count = 0
with open(temp_filt_fname, 'r') as bins:
    with open(filt_bins_file, 'w') as out:
        for bin in bins:
            coords = bin.split('\t')[0]
            (bchr, endpts) = coords.split(':')
            (bstart, bend) = endpts.split('-')
            bstart = int(bstart)
            bend = int(bend)
            if (bchr not in bad_areas) or \
                (sum(bad_areas[bchr][bstart:(bend + 1)]) == 0):
                kept_count += 1
                print >>out, bin,
            else:
                discarded_count += 1
os.remove(temp_filt_fname)
print '{0} bins kept, {1} bins discarded'.format(kept_count, discarded_count)

exit(0)
