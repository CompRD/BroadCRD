#!/usr/bin/tcsh

# script to assemble error corrected reads using BGI soap code.

set threads = 48
set K = 63

set bin_dir = ../bin
set working_dir = working
set config = na12878_corr.cfg
set out_head = $working_dir/na12878_k$K

# Create kmer spectrum

$bin_dir/SOAPdenovo-63mer pregraph -s $config -d 2 -p $threads -K $K -o $out_head >& pregraph.log

$bin_dir/SOAPdenovo-63mer contig -g $out_head -p $threads >& contig.log

$bin_dir/SOAPdenovo-63mer map -s $config -g $out_head -p $threads >& map.log

$bin_dir/SOAPdenovo-63mer scaff -g $out_head -p $threads >& scaff.log
