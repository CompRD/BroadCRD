#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program that trims all the reads in a FASTA file to a specified
# maximum length.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import sys
from BroadCRD.util.fasta import readFasta, writeFasta

if len(sys.argv) != 4:
    print >>sys.stderr, ('Usage ' + os.path.basename(sys.argv[0]) +
        ' input.fasta max_length output.fasta')
    exit(1)

input_filename = sys.argv[1]
max_length = int(sys.argv[2])
output_filename = sys.argv[3]

fasta_data = readFasta(input_filename)

for fasta_el in fasta_data:
    fasta_el[1] = fasta_el[1][0:min(max_length, len(fasta_el[1]))]

writeFasta(fasta_data, output_filename, range(len(fasta_data)))

exit(0)
