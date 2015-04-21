#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2013) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Count the bases in your BAM/SAM records.

import sys
from BroadCRD.util.Samtools import openbam

bamfile = sys.argv[1]

total_bases = 0

with openbam(bamfile, 'r') as bamdata:
    for samrec in bamdata:
        total_bases += len(samrec.split('\t')[9])
        
print 'total bases = {0}'.format(total_bases)

exit(0)
