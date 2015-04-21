#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Score the accuracy of PacBio's identification of the barcodes and simulated
# amplicons generated by MakeFakeMultiplexedBarcodes.py.

import sys

results_file = sys.argv[1]

correct_answer_count = 0
incorrect_answer_count = 0

first_line_skipped = False
with open(results_file, 'r') as results:
    for result_line in results:
        if first_line_skipped:
            fields = result_line.split('\t')
            correct_answer = ''.join(fields[0].split('_')[0:2])
            inferred_answer = fields[1]
            if correct_answer == inferred_answer:
                correct_answer_count += 1
            else:
                incorrect_answer_count += 1
        else:
            first_line_skipped = True

total_reads_barcoded = correct_answer_count + incorrect_answer_count

print 'total reads barcoded = {0}'.format(total_reads_barcoded)
print 'correct answers = {0} ({1:.2f})'.format(correct_answer_count,
    float(correct_answer_count) / total_reads_barcoded)
print 'incorrect answers = {0} ({1:.2f})'.format(incorrect_answer_count,
    float(incorrect_answer_count) / total_reads_barcoded)

exit(0)
