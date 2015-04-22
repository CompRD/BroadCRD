#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# This program splits a SAM file, paired reads go to the same output file.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import math
import os
import random
import sys

from BroadCRD.util.Samtools import openbam

if len(sys.argv) != 4:
	print >>sys.stderr, ('Usage: ' + os.path.basename(sys.argv[0]) +
		' input.sam output1.sam output2.sam')
	exit(1)

samfile = sys.argv[1]
outfile1 = sys.argv[2]
outfile2 = sys.argv[3]

read_names1 = set()
read_names2 = set()
with openbam(samfile, header=True) as sam, open(outfile1, 'w') as out1, \
    open(outfile2, 'w') as out2:
    for rec in sam:
        if rec[0] != '@':
            rname = rec.split('\t')[0]
            if rname in read_names1:
                print >>out1, rec,
            elif rname in read_names2:
                print >>out2, rec,
            else:
                if random.random() < 0.5:
                    read_names1.add(rname)
                    print >>out1, rec,
                else:
                    read_names2.add(rname)
                    print >>out2, rec,
        else:
            print >>out1, rec,
            print >>out2, rec,

exit(0)
