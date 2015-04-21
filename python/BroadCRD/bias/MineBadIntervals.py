#!/usr/bin/env python

import re
import sys

from BroadCRD.util.fasta import readFastaDict

def N50(lengths):
    lengths.sort()
    total = sum(lengths)
    half = 0
    for i, ln in enumerate(lengths):
        half += ln
        if (2 * half == total) and (i < len(lengths) - 1):
            return (ln + lengths[i + 1]) / 2
        if (2 * half >= total):
            return ln

def gc_count(sequence):
    count = 0
    for s in sequence:
        if s.upper() == 'G' or s.upper() == 'C':
            count += 1
    return count

def unamb_refsize(refdict):
    size = 0
    for seq in refdict.itervalues():
        for base in seq:
            if (base == 'A' or base == 'a' or base == 'C' or base == 'c' or
                base == 'G' or base == 'g' or base == 'T' or base == 't'):
                size += 1    
    return size

def main(intervals_filename, reference):
    #refdict = readFastaDict(reference)
    #ref_bases = unamb_refsize(refdict)
    ref_bases = 2870624073
    #Y_size = unamb_refsize({'chrY':refdict['chrY']})
    Y_size = 25652954
    lengths = []
    gc_content = []
    chrY_lengths = []
    with open(intervals_filename, 'r') as intervals_file:
        for interval in intervals_file:
            (loc, content) = interval.split('\t')
            start = int(re.search(':(\d+)-', loc).group(1))
            finish = int(re.search(':\d+-(\d+)', loc).group(1))
            if re.match('chrY', loc):
                chrY_lengths.append(finish - start)
            else:
                lengths.append(finish - start)
                gc_content.append(gc_count(content))
    total_bases_wY = sum(lengths + chrY_lengths)
    percent_bases_wY = float(total_bases_wY) / ref_bases
    print 'all - total={0} [{1}] N50={2}'.format(sum(lengths + chrY_lengths),
        percent_bases_wY, N50(lengths + chrY_lengths))
    total_bases = sum(lengths)
    percent_bases = float(total_bases) / (ref_bases - Y_size)
    total_gc = sum(gc_content)
    print 'sans chrY - total={0} [{1}] N50={2}'.format(total_bases,
        percent_bases, N50(lengths))
    print '\tGC={0}'.format(float(total_gc) / total_bases)
    return 0


if __name__ == '__main__':
    exit(main(sys.argv[1], sys.argv[2]))
