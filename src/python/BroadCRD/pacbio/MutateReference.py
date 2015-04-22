#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

import os
import random
import sys

from BroadCRD.util.fasta import readFasta, writeFasta

try:
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    sub_rate = float(sys.argv[3])
    ins_rate = float(sys.argv[4])
    del_rate = float(sys.argv[5])
except:
    print >>sys.stderr, ('Usage: ' + os.path.basename(sys.argv[0]) +
        ' input.fasta output.fasta substitution_rate insertion_rate' +
        ' deletion_rate')
    exit(1)

print 'reading ' + input_fasta
ref = readFasta(input_fasta)

max_indel_size = 100
nucleotides = 'ACGT'
nuc_mutations = {'A':'CGT', 'C':'AGT', 'G':'ACT', 'T':'ACG'}

for contig in ref:
    print 'mutating contig ' + contig[0]
    new_contig_string = []
    base = 0
    while base < len(contig[1]):
        # insertion?
        if random.random() < ins_rate:
            ins_seq_len = random.randint(1, max_indel_size)
            ins_seq = []
            for i in xrange(ins_seq_len):
                ins_seq.append(nucleotides[random.randint(0,
                    len(nucleotides) - 1)])
            print '{0}:{1} insert '.format(contig[0], base) + ''.join(ins_seq)
            new_contig_string += ins_seq
            
        # deletion?
        if random.random() < del_rate:
            skip = min(random.randint(1, max_indel_size), len(contig[1]) - base)
            print ('{0}:{1} delete '.format(contig[0], base)
                + ''.join(contig[1][base:base + skip]))
            base += skip
            
        # substitution?
        if base < len(contig[1]):
            if random.random() < sub_rate:
                old_base = contig[1][base]
                substitutions = nuc_mutations[old_base]
                new_base = substitutions[random.randint(0,
                    len(substitutions) - 1)]
                print ('{0}:{1} substitute '.format(contig[0], base)
                    + '{0}->{1}'.format(old_base, new_base))
                new_contig_string.append(new_base)
            else:
                new_contig_string.append(contig[1][base])
            
        base += 1
    contig[1] = ''.join(new_contig_string)

print 'writing ' + output_fasta
writeFasta(ref, output_fasta)

exit(0)
