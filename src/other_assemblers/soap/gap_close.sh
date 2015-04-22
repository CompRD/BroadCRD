#!/usr/bin/tcsh

# script to close gaps in a soap de novo assembly using BGI code.

set threads = 48
set K = 63

set bin_dir = ../bin
set working_dir = working
set config = na12878_corr.cfg
set out_head = $working_dir/na12878_k$K

# Close gaps

$bin_dir/GapCloser -a $out_head.scafSeq -b $config -o $out_head.scafSeq.gap_closed -l 155 -p 25 -t $threads >& gap_closer.log
