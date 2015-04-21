#!/usr/bin/tcsh

# script to error correct reads using BGI soap code.
# Expects pairs of fastq files, listed in the file reads2corr.lst

set threads = 48
set K = 27
set read_length = 250

set bin_dir = ../bin
set working_dir = working
set read_list = reads2corr.lst
set out_head = $working_dir/na12878

# Create kmer spectrum
$bin_dir/KmerFreq_HA -k $K -l $read_list -p $out_head -t $threads -L $read_length >& kmerfreq.log 

# Error correct reads
$bin_dir/Corrector_HA -k $K -l 3 -a 0 -e 0 -w 0 -q 35 -t $threads -j 1 $out_head.freq.gz $read_list >& correct.log

