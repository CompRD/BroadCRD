#!/usr/bin/env python2.6

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A quick and dirty program to find all the prefixes of a specified length
# in a FASTQ file, and display them in order of frequency (from lowest to
# highest). The intent is to help the user discover read barcodes.

# Original author: Michael Ross <mgross@broadinstitute.org>

from __future__ import print_function

import difflib
import operator
import sys

def compbase(base):
    if base == 'A':
        return 'T'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    elif base == 'T':
        return 'A'
    else:
        print('error - only A, C, G, or T are valid input', file=sys.stderr)
        exit(1)
    return

def revcomp(read):
    return ''.join([compbase(b) for b in read[::-1]])

def main():
    fastq = sys.argv[1]
    barcode_len = int(sys.argv[2])
    if len(sys.argv) >= 4:
        sim_cutoff = float(sys.argv[3])
    else:
        sim_cutoff = 1

    fastq_file = open(fastq, 'r')

    prefixes = {}
    prefix_keys = []

    fastq_line = fastq_file.readline()
    while fastq_line:
        if fastq_line[0] == '@':
            fastq_line = fastq_file.readline()
            prefix = fastq_line[0:min(len(fastq_line), barcode_len)].rstrip()

            if prefix in prefixes:
                prefixes[prefix] += 1
            elif sim_cutoff < 1:
                poss_matches = difflib.get_close_matches(prefix, prefix_keys,
                                                         n=1, cutoff=sim_cutoff)

                if poss_matches:
                    prefixes[poss_matches[0]] += 1
                else:
                    prefixes[prefix] = 1
                    prefix_keys.append(prefix)
                    print('.', end='')
                    sys.stdout.flush()

            fastq_file.readline()
            fastq_file.readline()
        fastq_line = fastq_file.readline()
    sorted_prefixes = sorted(prefixes.items(), key=operator.itemgetter(1))
    print()
    for k in sorted_prefixes:
        print(k)

if __name__ == '__main__':
    exit(main())
