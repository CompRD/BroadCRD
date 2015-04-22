#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# This program reorders the @SQ fields in a BAM file header to match the
# order of sequences in a reference FASTA file. Then it resorts the BAM
# reads in coordinate order (which will follow the new sequence order).
# Its main purpose is to enable BadCoverage to run on BAMs produced
# by the Pacific Biosciences pipeline.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 4:
    print >>sys.stderr, ('Usage: {0} reads.bam ref.fasta '
    'output.bam').format(os.path.basename(sys.argv[0]))
    exit(1)
    
bam_file = sys.argv[1]
ref_file = sys.argv[2]
new_sorted_bam_file = sys.argv[3]

if new_sorted_bam_file[-4:] == '.bam':
    new_sorted_bam_file = new_sorted_bam_file[:-4]
else:
    print >>sys.stderr, 'Output file name must end with ".bam"'
    exit(1)

# read in headers from the BAM
bam_headers = subprocess.Popen(['samtools', 'view', '-H', bam_file],
    stdout=subprocess.PIPE).communicate()[0].split('\n')

seq_headers = filter(lambda x: x[0:3] == '@SQ', bam_headers)

# read in sequence names from the reference
ref_names = subprocess.Popen(['grep', '^>', ref_file],
    stdout=subprocess.PIPE).communicate()[0].rstrip().split('\n')
ref_names = map(lambda x: re.match('^>(\S+)', x).group(1), ref_names)

# record locations of the sequences in the reference
ref_name_locations = {}
for loc, rn in enumerate(ref_names):
    ref_name_locations[rn] = loc

# reorder BAM sequence headers
def seq_cmp(x, y):
    x_name = re.match('^.*SN\:(\S+)', x).group(1)
    y_name = re.match('^.*SN\:(\S+)', y).group(1)
    return cmp(ref_name_locations[x_name], ref_name_locations[y_name])

seq_headers.sort(seq_cmp)

# insert the reorderd headers back in place
seq_head_start = 0
while (seq_head_start < len(bam_headers) and
    bam_headers[seq_head_start][0:3] != '@SQ'):
    seq_head_start += 1
if seq_head_start == len(bam_headers):
    print >> sys.stderr, 'Failed to find start of @SQ block'
    exit(1)
    
seq_head_end = seq_head_start + 1
while (seq_head_end < len(bam_headers) and
    bam_headers[seq_head_end][0:3] == '@SQ'):
    seq_head_end += 1
    
remainder = seq_head_end + 1
while (remainder < len(bam_headers) and bam_headers[remainder] != '@SQ'):
    remainder += 1
if remainder != len(bam_headers):
    print >>sys.stderr, '@SQ elements not located in a contiguous block'
    exit(1)

bam_headers = (bam_headers[0:seq_head_start] + seq_headers +
    bam_headers[seq_head_end:])

# create the output BAM file
temp_dir = tempfile.mkdtemp()
header_file = os.path.join(temp_dir, re.sub('.bam$', '.reorderedheader',
    os.path.basename(bam_file)))
with open(header_file, 'w') as header_filestream:
    print >>header_filestream, '\n'.join(bam_headers),

sam_file = os.path.join(temp_dir, re.sub('.bam$', '.sam',
    os.path.basename(bam_file)))
with open(sam_file, 'w') as sam_filestream:
    subprocess.check_call(['samtools', 'view', bam_file], stdout=sam_filestream)

sam_header_file = os.path.join(temp_dir, re.sub('.bam$', '.reorderedheader.sam',
    os.path.basename(bam_file)))
with open(sam_header_file, 'w') as sam_header_filestream:
    subprocess.check_call(['cat', header_file, sam_file],
        stdout=sam_header_filestream)

new_bam_file = os.path.join(temp_dir, re.sub('.bam$', '.reorderedheader.bam',
    os.path.basename(bam_file)))
with open(new_bam_file, 'w') as new_bam_filestream:
    subprocess.check_call(['samtools', 'view', '-bS', sam_header_file],
        stdout=new_bam_filestream)

subprocess.check_call(['samtools', 'sort', new_bam_file, new_sorted_bam_file])

os.remove(header_file)
os.remove(sam_file)
os.remove(sam_header_file)
os.remove(new_bam_file)
os.rmdir(temp_dir)

exit(0)
