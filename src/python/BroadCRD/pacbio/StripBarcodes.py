#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# This program strips barcodes from a FASTA or FASTQ of sequences using the 
# output of PacBio's PacBioBarcodeCCSDist software.

import re
import sys

from BroadCRD.util.fasta import readFastaDict, writeFasta
from BroadCRD.util.Fastq import readFastqDict, writeFastq

# note that the coordinates are relative to the version of the read that
# matched the barcode (so the left coords for an rc-matched read refer to
# the right side of the original read)

read_file = sys.argv[1]
barcode_results_file = sys.argv[2]
output_file = sys.argv[3]

fastq_data = False
if re.match('.*\.fastq$', read_file.lower()):
    if not re.match('.*\.fastq$', output_file.lower()):
        print >>sys.stderr, 'input is a FASTQ, output should be too'
        exit(1)
    input_reads = readFastqDict(read_file)
    fastq_data = True
else:
    input_reads = readFastaDict(read_file)

output_reads = []
fieldnames = None
with open(barcode_results_file, 'r') as barcode_results:
    for barcode_data in barcode_results:
        fields = barcode_data.rstrip().split('\t')
        if not fieldnames:
            fieldnames = fields
        else:
            data = dict(zip(fieldnames, fields))
            score_all = float(data['scoreAll'])
            score_right = float(data['scoreRight'])

            if (score_all > 30 and score_right < score_all * 0.6 and
                score_right > score_all * 0.4):
                read_name = data['fastaid']
                barcode = data['barcode']
                read_name_parse = re.match('(.*) \[revcomp\]$', read_name)
                left_seq_start = int(data['leftExtentEnd'])
                right_seq_end = int(data['rightExtentBegin']) - 1
                rc_match = False
                if read_name_parse:
                    rc_match = True
                    read_name = read_name_parse.group(1)
            
                if fastq_data:
                    (seq, qual) = input_reads[read_name]
                else:
                    seq = input_reads[read_name]
                if rc_match:
                    strip_seq = seq[-right_seq_end:-left_seq_start]
                    if fastq_data:
                        strip_qual = qual[-right_seq_end:-left_seq_start]
                else:
                    strip_seq = seq[left_seq_start:right_seq_end]
                    if fastq_data:
                        strip_qual = qual[left_seq_start:right_seq_end]
            
                if fastq_data:
                    output_reads.append(['{0}_{1}'.format(read_name, barcode),
                        strip_seq, strip_qual])
                else:
                    output_reads.append(['{0}_{1}'.format(read_name, barcode), 
                        strip_seq])

if fastq_data:
    writeFastq(output_reads, output_file)
else:
    writeFasta(output_reads, output_file)

exit(0)
