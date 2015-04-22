#!/usr/bin/tcsh

# script assemble data using platanus
# Expects pairs of fastq files

set threads = 48
set mem = 1000
set K = 32

set readA = orig_data/frag_reads_split.A.fastq
set readB = orig_data/frag_reads_split.B.fastq

set bin_dir = ../bin
set working_dir = working
set read_list = reads2corr.lst
set out_head = $working_dir/na12878

$bin_dir/platanus assemble -t $threads -m $mem -k $K -o $out_head -f $readA $readB  > assemble_log.txt

$bin_dir/platanus scaffold -t $threads -o $out_head -c ${out_head}_contig.fa -b ${out_head}_contigBubble.fa -IP2 $readA $readB > scaffold_log.txt

$bin_dir/platanus gap_close -t $threads -o $out_head -c ${out_head}_scaffold.fa -IP2 $readA $readB > gap_log.txt
