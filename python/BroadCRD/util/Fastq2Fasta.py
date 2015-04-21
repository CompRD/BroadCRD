#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2013) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Takes as input a FASTQ file and prints FASTA of the reads to standard out.

import textwrap
import sys
from BroadCRD.util.Fastq import readFastq

fastq_data = readFastq(sys.argv[1])

for fastq_rec in fastq_data:
    print '>{0}\n{1}'.format(fastq_rec[0], textwrap.fill(fastq_rec[1], 80))

exit(0)
