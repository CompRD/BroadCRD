#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Look up the genes associated with a list of intervals in RefSeq.

import re
import sys

from BroadCRD.bias.CombineIntervals import Interval, intersect

refseq_file = sys.argv[1]
intervals_file = sys.argv[2]

records = {}
with open(refseq_file, 'r') as refseq:
    for linenum, data in enumerate(refseq):
        if linenum == 0 and data[0] == '#':
            fieldnames = data[1:].rstrip().split('\t')
        else:
            fields = data.rstrip().split('\t')
            rec = dict(zip(fieldnames, fields))
            if re.match('^chr[1-9][0-9]{0,1}$|^chr[XY]$', rec['chrom']):
                rec['chrom'] = re.sub('^chr', '', rec['chrom'])
                if rec['chrom'] not in records:
                    records[rec['chrom']] = []
                records[rec['chrom']].append(Interval(rec['chrom'],
                    int(rec['txStart']), int(rec['txEnd']), rec['name2']))

with open(intervals_file, 'r') as intervals:
    for inter_text in intervals:
        (chr, range) = inter_text.rstrip().split(':')
        (start, stop) = range.split('-')
        inter = Interval(chr, int(start), int(stop))
        gene_names = set()
        for rec in records[chr]:
            if intersect(inter, rec):
                gene_names.add(rec.name)
        inter.name = ' '.join(gene_names)
        print inter


exit(0)
