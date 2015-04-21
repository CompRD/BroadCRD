#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A simple script for running BWA on a FASTQ file in a Picard-esque way.
# Note that it relies on having the "picard" script from
# crdutil/scripts/ in the path, and it explicitly uses the Picard versions
# of BWA rather than what might be in the user's path.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import optparse
import os
import re
import shlex
import subprocess
import sys

from BroadCRD.util.Picard import picard_exec
from BroadCRD.util.Samtools import openbam

# Remove multiple alignments introduced by BWASW to represent
# read chimerism - currently just keep the first alignment for
# each record.
def remove_multiple_assignments(samfile):                
    querysortedbamfile = re.sub('.sam$', '.querysorted.bam', samfile)
    picard_exec('SortSam', 'I=' + samfile, 'O=' + querysortedbamfile,
        'SO=queryname')
    with open(samfile, 'w') as fsamfile:
        with openbam(querysortedbamfile) as bam_records:
            current_read = None
            current_records = []
            for rec in bam_records:
                fields = rec.split('\t')
                if current_read != fields[0]:
                    if current_read:
                        print >>fsamfile, current_records[0],
                    current_read = fields[0]
                    current_records = [rec]
                else:   
                    current_records.append(rec)
    os.remove(querysortedbamfile)
    return

def main():
    oparser = optparse.OptionParser('usage: %prog [options] '
        'reads.fastq|fastq_list.txt reference.fasta num_threads '
        '[output_prefix]')
    oparser.add_option('--bwasw', action='store_true', default=False)
    (options, args) = oparser.parse_args()
    
    readfile = args[0]
    reference = args[1]
    threads = int(args[2])

    if readfile[-4:] == '.txt':
        fastqfile_list = open(readfile, 'r')
    else:
        fastqfile_list = [readfile]
    
    for fastqfile in fastqfile_list:
        fastqfile = fastqfile.rstrip()
        if len(args) > 3 and len(fastqfile_list) == 1:
            output_prefix = args[3]
        else:
            output_prefix = (re.match('(.*).fastq$',
                os.path.basename(fastqfile)).group(1))
        
        if options.bwasw:
            output_prefix += '.bwasw'
        else:
            output_prefix += '.bwashort'
    
        outsamfile = output_prefix + '.aligned.sam'
        mergedbamfile = re.sub('.sam$', '.merged.bam', outsamfile)
        sortedbamfile = re.sub('.bam$', '.sorted.bam', mergedbamfile)
        dupmarkedbamfile = re.sub('.bam$', '.duplicates_marked.bam', 
            sortedbamfile)
        
        bwa = '/seq/software/picard/current/3rd_party/bwa/bwa'
        
        if options.bwasw:
            bwasw_cmd = (bwa + ' bwasw -t {0} -f {1} {2} {3}').format(threads,
                outsamfile, reference, fastqfile)
            subprocess.check_call(shlex.split(bwasw_cmd))
            
            remove_multiple_assignments(outsamfile)
        else:
            saifile = output_prefix + '.sai'
            bwaaln_cmd = bwa + ' aln {0} -q 5 -l 32 -k 2 -t {3} -o 1 -f {1} {2}'
            bwaaln_cmd = bwaaln_cmd.format(reference, saifile, fastqfile, 
                threads)
            subprocess.check_call(shlex.split(bwaaln_cmd))
         
            bwasampe_cmd = (bwa + ' samse -f {0} {1} {2} {3}'
                .format(outsamfile, reference, saifile, fastqfile))
            subprocess.check_call(shlex.split(bwasampe_cmd))
            os.remove(saifile)
            
        unmappedbamfile = output_prefix + '.unmapped.bam'
        picard_exec('FastqToSam', 'FASTQ=' + fastqfile, 
            'QUALITY_FORMAT=Standard',
            'SAMPLE_NAME=' + fastqfile, 'O=' + unmappedbamfile)
               
        picard_exec('MergeBamAlignment', 'R=' + reference,
            'UNMAPPED=' + unmappedbamfile, 'ALIGNED=' + outsamfile,
            'OUTPUT=' + mergedbamfile, 'PAIRED_RUN=False',
            'VALIDATION_STRINGENCY=SILENT')
        os.remove(outsamfile)
        os.remove(unmappedbamfile)
           
        picard_exec('SortSam', 'I=' + mergedbamfile, 'O=' + sortedbamfile,
            'SO=coordinate')
        os.remove(mergedbamfile)
        
        dupmarkedmetricfile = re.sub('.bam$', '.duplicate_metrics', 
            dupmarkedbamfile)
        picard_exec('MarkDuplicates', 'I=' + sortedbamfile,
            'O=' + dupmarkedbamfile,
            'M=' + dupmarkedmetricfile)
        os.remove(sortedbamfile)
        
        gcbiasfile = re.sub('.bam$', '.gc_bias.pdf', dupmarkedbamfile)
        gcbiasmetricfile = re.sub('.pdf$', '.txt', gcbiasfile)
        picard_exec('CollectGcBiasMetrics', 'I=' + dupmarkedbamfile,
            'O=' + gcbiasmetricfile, 'CHART=' + gcbiasfile, 'R=' + reference)
    
    return 0

if __name__ == '__main__':
    exit(main())

