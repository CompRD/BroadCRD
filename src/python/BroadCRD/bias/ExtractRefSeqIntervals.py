#!/usr/bin/env python

# Extracts intervals from the RefSeq database - outputting a file for genes
# and for exons.

import re
import sys

refseq_file = sys.argv[1]
genes_file = sys.argv[2]
exons_file = sys.argv[3]

records = []
with open(refseq_file, 'r') as refseq:
    for linenum, data in enumerate(refseq):
        if linenum == 0 and data[0] == '#':
            fieldnames = data[1:].rstrip().split('\t')
        else:
            fields = data.rstrip().split('\t')
            rec = dict(zip(fieldnames, fields))
            if re.match('^chr[1-9][0-9]{0,1}$|^chr[XY]$', rec['chrom']):
                rec['chrom'] = re.sub('^chr', '', rec['chrom'])
                records.append(rec)

with open(genes_file, 'w') as genes:
    for rec in records:
        print >>genes, '{0}:{1}-{2}'.format(rec['chrom'], rec['txStart'], 
            rec['txEnd'])

with open(exons_file, 'w') as exons:
    for rec in records:
        exon_starts = rec['exonStarts'].split(',')
        exon_ends = rec['exonEnds'].split(',')
        exon_count = int(rec['exonCount'])
        for e in range(exon_count):
            print >>exons, '{0}:{1}-{2}'.format(rec['chrom'], exon_starts[e],
                exon_ends[e])

exit(0)
