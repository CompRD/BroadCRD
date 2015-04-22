#!/usr/bin/env python3.1

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program that splits all reads in a FASTA in half, producing two
# FASTAs.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import sys
from BroadCRD.util.fasta import readFasta, writeFasta

if len(sys.argv) != 2:
    print('Usage ' + os.path.basename(sys.argv[0]) +
        ' input.fasta', file=sys.stderr)
    exit(1)

input_filename = sys.argv[1]

if input_filename[-6:] != '.fasta':
    print('file does not appear to be a FASTA', file=sys.stderr)

output_filename_first = os.path.basename(input_filename)[:-6] + '_first.fasta'
output_filename_second = os.path.basename(input_filename)[:-6] + '_second.fasta'

fasta_data = readFasta(input_filename)

fasta_data_first = []
fasta_data_second = []
for fasta_el in fasta_data:
    read_len = len(fasta_el[1])
    fasta_data_first.append((fasta_el[0], fasta_el[1][0:(read_len // 2)]))
    fasta_data_second.append((fasta_el[0], fasta_el[1][(read_len // 2):]))

writeFasta(fasta_data_first, output_filename_first, range(len(fasta_data)))
writeFasta(fasta_data_second, output_filename_second, range(len(fasta_data)))

exit(0)
