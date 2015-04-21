#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# This program downsamples a SAM file, maintaining subread affiliations
# (all subreads from a well are both present or both absent in the output).

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import random
import re
import sys

from BroadCRD.util.Samtools import openbam

try:
    samfile = sys.argv[1]
    fraction = float(sys.argv[2])
    output = sys.argv[3]
except:
	print >>sys.stderr, ('Usage: ' + os.path.basename(sys.argv[0]) +
		' input.sam fraction output.sam')
	exit(1)

excluded_read_names = set()
kept_read_names = set()
with openbam(samfile, header=True) as sam, openbam(output, 'w') as outfile:
    for rec in sam:
        if rec[0] != '@':
            rname = rec.split('\t')[0]
            rmatch = re.match('^(.*/.*)/.*$', rname)
            if rmatch:
                rname = rmatch.group(1)

            if rname in kept_read_names:
                print >>outfile, rec,
            elif rname not in excluded_read_names:
                if random.random() < fraction:
                    kept_read_names.add(rname)
                    print >>outfile, rec,
                else:
                    excluded_read_names.add(rname)
        else:
            print >>outfile, rec,
            

exit(0)
