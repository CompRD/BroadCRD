#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A simple script for running BWA on a BAM file. Parameters inspired by
# Picard. Note that it relies on having the "picard" script from
# crdutil/scripts/ in the path, and it explicitly uses the Picard versions
# of BWA rather than what might be in the user's path.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import shlex
import subprocess
import sys

from BroadCRD.util.Picard import picard_exec

def main(argv):
    try:
        bamfile = argv[1]
        reference = argv[2]
        paired = argv[3]
        threads = int(argv[4])
        if paired != 'true' and paired != 'false':
            print >>sys.stderr, ('The only valid paired values are "true" or '
                '"false"')
            return 1
        if len(argv) > 5:
            output_prefix = argv[5]
        else:
            output_prefix = (re.match('(.*).bam', os.path.basename(bamfile))
                .group(1))
            output_prefix = re.sub('\.aligned', '', output_prefix)
            output_prefix = re.sub('\.duplicates_marked', '', output_prefix)
            
    except:
        print >>sys.stderr, ('Usage: {0} reads.bam reference.fasta '
            'paired threads [tag]'.format(os.path.basename(argv[0])))
        return 1
    
    unmappedbamfile = output_prefix + '.unmapped.bam'
    outsamfile = output_prefix + '.aligned.sam'
    mergedbamfile = re.sub('.sam$', '.merged.bam', outsamfile)
    sortedbamfile = re.sub('.bam$', '.sorted.bam', mergedbamfile)
    dupmarkedbamfile = re.sub('.bam$', '.duplicates_marked.bam', sortedbamfile)
    
    if os.path.exists(dupmarkedbamfile):
        print >>sys.stderr, ('{0} exists - should I skip {1}?'
            .format(dupmarkedbamfile, bamfile))
        skip = 'x'
        while skip != 'y' and skip != 'n' and skip != '':
            skip = raw_input('(y/n) [y]: ')
        if skip == 'y' or skip == '':
            return 0
    
    picard_exec('RevertSam', 'I=' + bamfile, 'O=' + unmappedbamfile,
        'RESTORE_ORIGINAL_QUALITIES=true', 'REMOVE_DUPLICATE_INFORMATION=true',
        'REMOVE_ALIGNMENT_INFORMATION=true')
    
    fastq1file = output_prefix + '.1.fastq'
    fastq2file = output_prefix + '.2.fastq'
    picard_exec('SamToFastq', 'I=' + unmappedbamfile, 'F=' + fastq1file,
        'F2=' + fastq2file)
    
    bwa = '/seq/software/picard/current/3rd_party/bwa/bwa'
    sai1file = output_prefix + '.1.sai'
    sai2file = output_prefix + '.2.sai'
    bwaaln_cmd = bwa + ' aln {0} -q 5 -l 32 -k 2 -t {3} -o 1 -f {1} {2}'
    bwaaln1_cmd = bwaaln_cmd.format(reference, sai1file, fastq1file, threads)
    bwaaln2_cmd = bwaaln_cmd.format(reference, sai2file, fastq2file, threads)
    subprocess.check_call(shlex.split(bwaaln1_cmd))
    subprocess.check_call(shlex.split(bwaaln2_cmd))
 
    bwasampe_cmd = (bwa + ' sampe -P -f {0} {1} {2} {3} {4} {5}'
        .format(outsamfile, reference, sai1file, sai2file, fastq1file,
        fastq2file))
    subprocess.check_call(shlex.split(bwasampe_cmd))
    os.remove(sai1file)
    os.remove(sai2file)
    os.remove(fastq1file)
    os.remove(fastq2file)
       
    picard_exec('MergeBamAlignment', 'R=' + reference,
        'UNMAPPED=' + unmappedbamfile, 'ALIGNED=' + outsamfile,
        'OUTPUT=' + mergedbamfile, 'PAIRED_RUN=' + paired,
        'VALIDATION_STRINGENCY=SILENT')
    os.remove(outsamfile)
    os.remove(unmappedbamfile)
       
    picard_exec('SortSam', 'I=' + mergedbamfile, 'O=' + sortedbamfile,
        'SO=coordinate')
    os.remove(mergedbamfile)
    
    dupmarkedmetricfile = re.sub('.bam$', '.duplicate_metrics', 
        dupmarkedbamfile)
    picard_exec('MarkDuplicates', 'I=' + sortedbamfile, 'O=' + dupmarkedbamfile,
        'M=' + dupmarkedmetricfile)
    os.remove(sortedbamfile)
    
    return 0

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1][-4:] == '.txt':
        with open(sys.argv[1], 'r') as bam_list:
            for bam in bam_list:
                bam = bam.rstrip()
                print 'aligning ' + bam 
                ret = main([sys.argv[0], bam] + sys.argv[2:])
                if ret != 0:
                    exit(ret)
        exit(0)
    else:
        exit(main(sys.argv))
