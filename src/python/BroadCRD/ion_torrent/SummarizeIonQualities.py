#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program to summarize the quality scores in an Ion Torrent run and
# their reliability. It accepts the analysis name as its one and only
# argument.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import getpass
import glob
import optparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

def locate_analysis_dir(analysis_name):
    analysis_root = '/seq/lagoon_results/analysis/output/Broad'
    exact = os.path.join(analysis_root, analysis_name)
    if os.path.exists(exact):
        all_possibilities = [exact]
    else:
        all_possibilities = []

    all_possibilities += glob.glob(exact + '_*')

    alt_analysis_root = '/seq/ion_analysis_archive/TorrentServer/analysis'
    alt_exact = os.path.join(alt_analysis_root, analysis_name)
    alt_auto = os.path.join(alt_analysis_root, 'Auto_' + analysis_name)

    if os.path.exists(alt_exact):
        all_possibilities += [alt_exact]
    if os.path.exists(alt_auto):
        all_possibilities += [alt_auto]

    all_possibilities = (all_possibilities + glob.glob(alt_exact + '_*') +
        glob.glob(alt_auto + '_*'))

    if len(all_possibilities) > 1:
        all_possibilities.sort()
        print >>sys.stderr, ('Please be more specific or use the '
            '--analysis_dir option to specify one of these possible analyses:'
            '\n' + '\n'.join(all_possibilities))
        return None
    elif len(all_possibilities) == 0:
        print >> sys.stderr, ('Unable to locate analysis directory for '
            + analysis_name)
        return None
    else:
        return all_possibilities[0]

def extract_read_data(analysis_dir, ionroot_dir, ref):
    user_temp_root = '/broad/hptmp/' + getpass.getuser()
    if not os.path.exists(user_temp_root):
        os.mkdir(user_temp_root)
    temp_dir = tempfile.mkdtemp(dir=user_temp_root)

    all_sam_files = glob.glob(os.path.join(analysis_dir, '*.sam'))
    all_sam_files += glob.glob(os.path.join(analysis_dir, '*.bam'))
    if len(all_sam_files) > 1:
        print >> sys.stderr, ('Found {0} SAM/BAM files in {1}, rather than 1'
            .format(len(all_sam_files), analysis_dir))
        exit(1)
    sam_file = all_sam_files[0]

    if not os.path.exists(sam_file):
        print >> sys.stderr, 'No SAM or BAM file found, exiting.'
        exit(1)

    # figure out which root directory we're using
    if not ionroot_dir:
        if re.match('^/seq/ion_analysis_archive/TorrentServer', analysis_dir):
            ionroot_dir = '/seq/ion_analysis_archive/TorrentServer'
        elif re.match('^/seq/lagoon_results', analysis_dir):
            ionroot_dir = '/seq/lagoon_results'
        else:
            print >>sys.stderr, ('Unknown analysis root directory ' +
                analysis_dir)
            exit(1)

    # find the reference FASTA
    if not ref:
        bfastout_file = os.path.join(analysis_dir, 'bfastAnalysis_err.txt')
        if os.path.exists(bfastout_file):
            with open(bfastout_file, 'r') as bfastout:
                for bfline in bfastout:
                    fastam = re.match('Validating fastaFileName (.*)\.', bfline)
                    if fastam:
                        ref = os.path.join(ionroot_dir, fastam.group(1)[1:])
                        break
        else:
            samtools_cmd = ['samtools', 'view', '-H']
            if sam_file[-4:] == '.sam':
                samtools_cmd += ['-S']
            samtools_cmd += [sam_file]

            refregexp = re.compile('^@PG\s+ID:tmap.*-f'
                '\s+/results/(referenceLibrary/\S+)')

            sam_header = (subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE)
                .communicate()[0].split('\n'))

            for sam_line in sam_header:
                if sam_line[0] != '@':
                    break
                else:
                    regresult = re.match(refregexp, sam_line)
                    if regresult:
                        ref = os.path.join(ionroot_dir, regresult.group(1))
                        break

    if not ref or not os.path.exists(ref):
        print >> sys.stderr, 'Failed to extract valid reference FASTA'
        exit(1)

    # convert to CRD formats
    subprocess.check_call(['SAM2CRDDump', 'SAM=' + sam_file, 'REF_FASTA=' + ref,
        'OUT_HEAD=' + os.path.join(temp_dir, 'SIQ'), 'NO_HEADER=True',
        'WRITE_PAIRS=False', 'WRITE_NAMES=False'])
    reads = os.path.join(temp_dir, 'SIQ.fastb')
    quals = os.path.join(temp_dir, 'SIQ.qualb')
    aligns = os.path.join(temp_dir, 'SIQ.qltout')
    return (reads, ref, aligns, quals, temp_dir)

def print_stats(q_stats, prefix, output_jira=False, rsquared=None):
    if not output_jira:
        print(prefix + ' bases min quality = ' + q_stats[0])
        print(prefix + ' bases lower quartile quality = ' + q_stats[1])
        print(prefix + ' bases median quality = ' + q_stats[2])
        print(prefix + ' bases upper quartile quality = ' + q_stats[3])
        print(prefix + ' bases max quality = ' + q_stats[4])
    else:
        print ('||R^2 between qualities and errors' +
            '||' + prefix + ' bases min quality' +
            '||' + prefix + ' bases median quality' +
            '||' + prefix + ' bases max quality||')
        print ('|' + rsquared +
            '|' + q_stats[0] +
            '|' + q_stats[2] +
            '|' + q_stats[4] + '|')
    return

def summarize_quals(reads, ref, aligns, quals, print_table, read_pos_report,
    output_jira):
    evq_pipe = subprocess.Popen(['EvaluateQuals', 'READS=' + reads,
                                 'REF=' + ref,
                                 'ALIGNS=' + aligns, 'QUALS=' + quals,
                                 'NO_HEADER=True',
                                 'MAX_ALIGN_ERROR_RATE=0.4',
                                 'READ_POS_REPORT=' + str(read_pos_report)],
                                 stdout=subprocess.PIPE)
    evq_out = evq_pipe.communicate()[0]
    rsquared = None
    all_stats = None
    aligned_stats = None
    unaligned_stats = None
    for evq_line in evq_out.splitlines():
        if print_table:
            print evq_line
        rsquared_match = re.match('^R2ideal\s+(\S+)', evq_line)
        all_stats_match = re.match('^Predicted quality of all.*:\s+'
                                   '(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$',
                                   evq_line)
        aligned_stats_match = re.match('^Predicted quality of aligned.*:\s+'
                                       '(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$',
                                       evq_line)
        unaligned_stats_match = re.match('^Predicted quality of unaligned.*:\s'
                                         '(\d+)\s+(\d+)\s+(\d+)\s+(\d+)'
                                         '\s+(\d+)$', evq_line)
        if rsquared_match:
            rsquared = rsquared_match.group(1)
        elif all_stats_match:
            all_stats = all_stats_match.groups()
        elif aligned_stats_match:
            aligned_stats = aligned_stats_match.groups()
        elif unaligned_stats_match:
            unaligned_stats = unaligned_stats_match.groups()

    if rsquared and not output_jira:
        print('R^2 between qualities and errors = ' + rsquared)
    if all_stats:
        print_stats(all_stats, 'all', output_jira, rsquared)
    if aligned_stats and not output_jira:
        print_stats(aligned_stats, 'aligned')
    if unaligned_stats and not output_jira:
        print_stats(unaligned_stats, 'unaligned')

    return

def main(argv):
    oparser = optparse.OptionParser(usage='%prog [options] IonAnalysisName')
    oparser.add_option('--full_table', action='store_true', default=False)
    oparser.add_option('--read_pos_report', action='store_true', default=False)
    oparser.add_option('--analysis_dir', action='store', type='string',
        default=None)
    oparser.add_option('--ionroot_dir', action='store', type='string',
        default=None)
    oparser.add_option('--output_jira', action='store_true', default=False)
    oparser.add_option('--reference', action='store', type='string',
        default=None)
    (options, argv) = oparser.parse_args(argv)

    if not options.analysis_dir:
        analysis_name = argv[1]
        analysis_dir = locate_analysis_dir(analysis_name)
    else:
        analysis_dir = options.analysis_dir

    if not analysis_dir:
        exit(1)

    print 'Processing ' + analysis_dir

    # Extract information
    (reads, ref, aligns, quals, temp_dir) = extract_read_data(analysis_dir,
        options.ionroot_dir, options.reference)

    # Summarize quality statistics
    summarize_quals(reads, ref, aligns, quals, options.full_table,
        options.read_pos_report, options.output_jira)

    # Remove temp files
    shutil.rmtree(temp_dir)

    # Profit
    exit(0)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
