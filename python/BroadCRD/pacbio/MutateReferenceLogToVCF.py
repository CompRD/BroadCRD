#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program that mutates the stdout stream produced by MutateReference.py
# to the VCF-formatted list of variants we would expect to get from reads
# aligned to that reference. Currently only works for SNPs.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import sys

try:
    log_file = sys.argv[1]
    vcf_file = sys.argv[2]
except:
    print >>sys.stderr, ('Usage: ' + os.path.basename(sys.argv[0])
        + 'input.log output.vcf')
    exit(1)

sub_re = re.compile('(.+):(\d+) substitute (.)->(.)')
del_re = re.compile('(.+):(\d+) delete')
ins_re = re.compile('(.+):(\d+) insert')
with open(log_file, 'r') as log:
    with open(vcf_file, 'w') as vcf:
        for line in log:
            sub_match = re.match(sub_re, line)
            del_match = re.match(del_re, line)
            ins_match = re.match(ins_re, line)
            
            if sub_match:
                chr = sub_match.group(1)
                pos = int(sub_match.group(2)) + 1
                ref = sub_match.group(4)
                alt = sub_match.group(3)
                print >>vcf, '{0}\t{1}\t.\t{2}\t{3}'.format(chr, pos, ref, alt)
            if del_match:
                print >>sys.stderr, 'deletions are unsupported'
                exit(1)
            if ins_match:
                print >>sys.stderr, 'insertions are unsuppored'
                exit(1)

exit(0)