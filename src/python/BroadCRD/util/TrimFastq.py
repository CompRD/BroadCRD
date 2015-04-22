#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program that removes a specified number of bases from the beginning
# of all the reads in the FASTQ file.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import sys

try:
    fastq_file = sys.argv[1]
    trim = int(sys.argv[2])
except:
    print >>sys.stderr, ('Usage: {0} reads.fastq '
        'num_bases').format(os.path.basename(sys.argv[0]))
    exit(1)

with open(fastq_file, 'r') as fastq_file:
    cur_block = ''
    block_valid = False
    file_line = 0
    for fastq_line in fastq_file:
        fastq_line = fastq_line.rstrip()
        
        if ((file_line % 4 == 0 and fastq_line[0] != '@') or
            (file_line % 4 == 2 and fastq_line[0] != '+')):
            print >>sys.stderr, 'FASTQ formatting error'
            exit(1)

        if file_line % 4 == 0:
            if block_valid:
                block_lines = cur_block.split('\n')

                if len(block_lines[1]) != len(block_lines[3]):
                    print >>sys.stderr, 'Output error'
                    print >>sys.stderr, cur_block
                    exit(1)
                print cur_block,
            cur_block = ''
            block_valid = True
        elif file_line % 4 != 2:
            if len(fastq_line) > trim:
                fastq_line = fastq_line[trim:]
            else:
                block_valid = False

        cur_block += fastq_line + '\n'
        file_line += 1

    if block_valid:
        print cur_block,

exit(0)     
