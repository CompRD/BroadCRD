#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A messy program for taking the output of PacBio's alpha barcoding software
# (run on CCS reads) and using it to chop up the relevant
# filtered_subreads.fast{a,q} files.

import re
import string
import sys

from BroadCRD.util.fasta import readFasta, writeFasta
from BroadCRD.util.Fastq import readFastq, writeFastq

def read_data(bc_file):
    bc_reads = {}
    bc_counts = {}
    bc_score_totals = {}
    bc_revcomps = {}
    
    with open(bc_file, 'r') as bc_data:
        skip_first = False
        for bc_line in bc_data:
            if not skip_first:
                skip_first = True
            else:
                bc_fields = bc_line.split('\t')
                score = float(bc_fields[3])
                barcode = int(re.match('F(\d+)', bc_fields[1]).group(1))
                read = re.sub(' \[revcomp\]$', '', bc_fields[0])
                if read != bc_fields[0]:
                    bc_revcomps[read] = True
                else:
                    bc_revcomps[read] = False
                bc_reads[read] = barcode
                if barcode not in bc_counts:
                    bc_counts[barcode] = 0
                    bc_score_totals[barcode] = 0
                bc_counts[barcode] += 1
                bc_score_totals[barcode] += score
    
    bc_list = list(bc_counts.keys())
    bc_list.sort()
    print 'barcode histogram:'
    sum = 0
    for bc in bc_list:
        print '{0}\t{1}\t{2}'.format(bc, bc_counts[bc], bc_score_totals[bc] /
            bc_counts[bc])
        sum += bc_counts[bc]
    print '{0} barcoded CCS reads'.format(sum)
    print 'mean count per barcode = {0}'.format(float(sum) / len(bc_counts))
    
    return (bc_reads, bc_revcomps)

def summarize_counts(subread_count):
    rounds_hist = {}
    for c in subread_count.keys():
        if subread_count[c] not in rounds_hist:
            rounds_hist[subread_count[c]] = 0
        rounds_hist[subread_count[c]] += 1
    rounds_list = list(rounds_hist.keys())
    rounds_list.sort()
    print 'subread count histogram:'
    for r in rounds_list:
        print '{0}\t{1}'.format(r, rounds_hist[r])
    return

def assign_subreads(subreads_file, bc_reads, bc_revcomps, output_prefix,
    cleanup_list):
    if subreads_file[-6:] == '.fasta':
        subreads = readFasta(subreads_file)
    elif subreads_file[-6:] == '.fastq':
        subreads = readFastq(subreads_file)
    else:
        print >>sys.stderr, 'bailing out, unknown subreads format'
        exit(1)
    
    if cleanup_list:
        correctly_barcoded = set(map(lambda x: x.rstrip(),
            open(cleanup_list, 'r')))
    
    split_subreads = {}
    total_subreads = 0
    unmatched_subreads = 0
    unmatched_subread_count = {}
    matched_subread_count = {}
    for s in subreads:
        total_subreads += 1        
        ccs_name = re.sub('\d+_\d+$', 'ccs', s[0])
        if (ccs_name in bc_reads and
            (not cleanup_list or ccs_name in correctly_barcoded)):
            if bc_revcomps[ccs_name]:
                s[0] = s[0] + ' [revcomp]'
                s[1] = s[1][::-1].translate(string.maketrans
                    ('ACTGactg', 'TGACtgac'))
                if len(s) > 2:
                    s[2] = s[2][::-1]
            bc = bc_reads[ccs_name]
            if bc not in split_subreads:
                split_subreads[bc] = []
            split_subreads[bc].append(s)
            if ccs_name not in matched_subread_count:
                matched_subread_count[ccs_name] = 0
            matched_subread_count[ccs_name] += 1
        else:
            unmatched_subreads += 1
            if ccs_name not in unmatched_subread_count:
                unmatched_subread_count[ccs_name] = 0
            unmatched_subread_count[ccs_name] += 1
            

    if output_prefix:
        for bc in split_subreads.keys():
            if subreads_file[-6:] == '.fasta':
                writeFasta(split_subreads[bc],
                    '{0}_{1}.fasta'.format(output_prefix, bc))
            elif subreads_file[-6:] == '.fastq':
                writeFastq(split_subreads[bc],
                    '{0}_{1}.fastq'.format(output_prefix, bc))
            else:
                print >>sys.stderr, ('bailing out, should have bailed out '
                    'earlier!')
                exit(1)

    print 'total subreads: {0}'.format(total_subreads)
    print 'unmatched subreads: {0}'.format(unmatched_subreads)
    summarize_counts(unmatched_subread_count)
    print 'matched subreads: {0}'.format(total_subreads - unmatched_subreads)
    summarize_counts(matched_subread_count)
    
    return        

def main():
    bc_file = sys.argv[1]

    subreads_fasta = None
    output_prefix = None
    cleanup_list = None
    
    if len(sys.argv) > 2:
        subreads_fasta = sys.argv[2]    
    if len(sys.argv) > 3:
        output_prefix = sys.argv[3]
    if len(sys.argv) > 4:
        cleanup_list = sys.argv[4]
        
    (bc_reads, bc_revcomps) = read_data(bc_file)
    
    if subreads_fasta:
        assign_subreads(subreads_fasta, bc_reads, bc_revcomps, output_prefix,
            cleanup_list)
    
    return 0

if __name__ == '__main__':
    exit(main())
