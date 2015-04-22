#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Print out size statistics from a collection of intervals.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import optparse
import math
import re
import sys

import numpy

from BroadCRD.util.fasta import readFastaDict

def n50(size_list):
    size_list.sort()
    halfway = len(size_list) / 2
    if len(size_list) % 2 == 1:
        return size_list[halfway]
    else:
        return (size_list[halfway - 1] + size_list[halfway]) / 2

def main():
    oparser = optparse.OptionParser()
    oparser.add_option('--reference', action='store', type='string')
    oparser.add_option('--detail', action='store_true', default=False)
    (options, args) = oparser.parse_args()
    
    ref = None
    if options.reference:
        ref = readFastaDict(options.reference)
    
    for intervals_file in args:
        print '--' + intervals_file + '--'
        size_count = []
        base_count = []
        at_count = []
        gc_count = []
        hp_N50 = []
        total_bases = 0
        size_list = []
        hp_size_list = {}
    
        if options.detail:
            print 'Interval,GC fraction,N50 homopolymer length'
    
        with open(intervals_file, 'r') as intervals:
            for i in intervals:
                i_match = re.match('^(.*):(\d+)-(\d+)', i)
                if i_match:
                    contig = i_match.group(1)
                    start = int(i_match.group(2))
                    stop = int(i_match.group(3))
                    if start < stop:
                        log_size = int(math.floor(math.log10(stop - start)))
                        if log_size >= len(size_count):
                            ext_size = (log_size - len(size_count) + 1)
                            size_count.extend([0]*ext_size)
                            base_count.extend([0]*ext_size)
                            at_count.extend([0]*ext_size)
                            gc_count.extend([0]*ext_size)
                            hp_N50.extend([0]*ext_size)
                        size_count[log_size] += 1
                        base_count[log_size] += stop - start
                        size_list += [stop - start]*(stop-start)
    
                        local_at_count = None
                        local_gc_count = None
                        local_hp_N50 = None
    
                        total_bases += stop - start
    
                        if ref:
                            current_hp_length = 0
                            current_hp_base = ref[contig][start]
                            int_hps = []
    
                            local_at_count = 0
                            local_gc_count = 0
    
                            for b in ref[contig][start:stop].lower():
                                if b == current_hp_base:
                                    current_hp_length += 1
                                else:
                                    int_hps += ([current_hp_length] *
                                        current_hp_length)
                                    current_hp_base = b
                                    current_hp_length = 1
    
                                if b == 'a' or b == 't':
                                    local_at_count += 1
                                elif b == 'g' or b == 'c':
                                    local_gc_count += 1
                            int_hps += [current_hp_length] * current_hp_length
                            local_hp_N50 = numpy.median(int_hps)
                        
                            if log_size not in hp_size_list:
                                hp_size_list[log_size] = []
                            hp_size_list[log_size] += int_hps

                            at_count[log_size] += local_at_count
                            gc_count[log_size] += local_gc_count
                            hp_N50[log_size] += local_hp_N50
    
                        if options.detail:
                            local_gc_percent = (local_gc_count /
                                float(stop - start))
                            print '{0}:{1}-{2},{3},{4}'.format(contig, start, 
                                stop, local_gc_percent, local_hp_N50)
                    else:
                        print >>sys.stderr, 'skipping {0}'.format(i.rstrip())

    
        for s, c in enumerate(size_count):
            print '>={0}\t-\t{1}\t{2:0.0f}%\t[{3}]'.format(10**s, c,
                float(base_count[s]) / total_bases * 100, base_count[s]),
            if ref and c > 0:
                print ('AT={0:0.0f}%\tGC={1:0.0f}%\tHPN50={2:.1f}\t'
                    'HPMEAN={3:.1f}'.format(
                    float(at_count[s]) / base_count[s] * 100,
                    float(gc_count[s]) / base_count[s] * 100,
                    n50(hp_size_list[s]), numpy.mean(hp_size_list[s])))
            else:
                print
        print 'Total Bases = {0}'.format(total_bases)            
        print 'Intervals N50 = {0}'.format(n50(size_list))
        print 'Intervals Max = {0}'.format(max(size_list))
        print 'Intervals Min = {0}'.format(min(size_list))
        if ref:
            full_hp_size_list = reduce(lambda x, y: x +
                hp_size_list[y], hp_size_list, [])
            print 'Homopolymer N50 = {0}'.format(n50(full_hp_size_list))
    
    return 0

if __name__ == '__main__':
    exit(main())

