#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# An unholy mix of processing required to get contigs to align via BWASW.
# Contigs longer than 100 kbases are split and ambiguous base codes are
# set to N (that's more of a Picard issue).
# All hail Cthulu!

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import optparse
import re
import sys

from BroadCRD.util.fasta import readFasta, writeFasta

optparser = optparse.OptionParser()
optparser.add_option('--max_size', type='int', default=100000)
(options, args) = optparser.parse_args()

fasta_data = readFasta(args[0])
new_data = []

for c, contig in enumerate(fasta_data):
    new_seq = re.sub('[^ACGTNacgtn]', 'N', contig[1])
    split_list = []
    if len(new_seq) > options.max_size:
        remainder = new_seq
        while len(remainder) > options.max_size:
            new_contig_name = '{0}_{1}'.format(contig[0], len(split_list))
            if len(remainder) > 2 * options.max_size:
                split_end = options.max_size
            else:
                split_end = len(remainder) / 2
            new_contig_seq = remainder[0:split_end]
            remainder = remainder[split_end:]
            split_list.append((new_contig_name, new_contig_seq))
        split_list.append(('{0}_{1}'.format(contig[0], len(split_list)),
            remainder))
    else:
        split_list = [(contig[0], new_seq)]
    new_data += split_list

writeFasta(new_data, args[1])

exit(0)
