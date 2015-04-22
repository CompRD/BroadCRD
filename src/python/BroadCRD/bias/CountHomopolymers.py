#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Count up the homopolymer runs of various lengths (and bases) in a reference
# FASTA.

# Original author: Michael G. Ross <mgross@broadinstitute.org>


import sys

from BroadCRD.util.fasta import readFasta

reference = sys.argv[1]

ref_contigs = readFasta(reference)

hps = {'A':{}, 'C':{}, 'G':{}, 'T':{}}

for contig in ref_contigs:
    prev_base = None
    prev_base_count = 0
    for base in contig[1]:
        if prev_base and prev_base == base.upper():
            prev_base_count += 1
        else:
            if prev_base:
                if prev_base not in hps:
                    hps[prev_base] = {}

                if prev_base_count not in hps[prev_base]:
                    hps[prev_base][prev_base_count] = 0
                hps[prev_base][prev_base_count] += 1
                
            prev_base = base.upper()
            prev_base_count = 1

for base in hps.keys():
    print base
    counts = list(hps[base].keys())
    counts.sort()
    for c in counts:
        print '\t{0} - {1}'.format(c, hps[base][c])

exit(0)