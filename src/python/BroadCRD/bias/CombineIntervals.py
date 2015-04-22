#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

import optparse
import re
import sys

class Interval:
    def __init__(self, contig, start, finish, name=None):
        self.contig = contig
        self.start = start
        self.finish = finish
        self.name = name
        
    def __cmp__(self, other):
        if self.contig == other.contig:
            return cmp(self.start, other.start)
        else:
            return cmp(self.contig, other.contig)
        
    def __str__(self):
        rep_str = '{0}:{1}-{2}'.format(self.contig, self.start, self.finish)
        if self.name:
            rep_str += ' #' + self.name 
        return rep_str
        
    def size(self):
        return self.finish - self.start

def overlap(a, b):
    return ((a.contig == b.contig) and
        ((a.start <= b.start and b.start < a.finish) or
        (b.start <= a.start and a.start < b.finish)))

def intersect(a, b):
    if overlap(a, b):
        return Interval(a.contig, max(a.start, b.start),
            min(a.finish, b.finish))
    else:
        return None

def union(a, b):
    if overlap(a, b):
        return Interval(a.contig, min(a.start, b.start),
            max(a.finish, b.finish))
    else:
        return None

def subtract(a, b):
    if overlap(a, b):
        res = []
        if a.start < b.start:
            res += [Interval(a.contig, a.start, b.start)]
        if b.finish < a.finish:
            res += [Interval(a.contig, b.finish, a.finish)]
        return res
    else:
        return None

def load_intervals(intervals_fname):
    intervals = []
    with open(intervals_fname, 'r') as interval_lines:
        for i in interval_lines:
            itext = i.rstrip().split('\t')[0]
            if re.match('^\S+:\d+-\d+$', itext):
                (contig, span) = i.rstrip().split('\t')[0].split(':')
                (start, finish) = span.split('-')
                intervals.append(Interval(contig, int(start), int(finish)))
    intervals.sort()
    return intervals

def merge_overlaps(intervals):
    i = 0
    while i < len(intervals) - 1:
        u = union(intervals[i], intervals[i + 1])
        if u:
            intervals[i:(i + 2)] = [u]
        else:
            i += 1
    return
    
def intersect_overlaps(intervals):
    intersections = []
    i = 0
    while i < len(intervals) - 1:
        j = i + 1
        while j < len(intervals) and overlap(intervals[i], intervals[j]):
            intersections += [intersect(intervals[i], intervals[j])]
            j += 1
        i += 1
    return intersections

def subtract_overlaps(intervals_a, intervals_b):
    results = list(intervals_a)
    r = 0
    b = 0
    while r < len(results):
        while (b < len(intervals_b) and
            (intervals_b[b].contig < results[r].contig or
            (intervals_b[b].contig == results[r].contig and
            intervals_b[b].finish <= results[r].start))):
            b += 1
        
        bp = b
        while (bp < len(intervals_b) and
            intervals_b[bp].contig == results[r].contig and
            intervals_b[bp].start < results[r].finish and
            not overlap(results[r], intervals_b[bp])):
            bp += 1

        if bp < len(intervals_b) and overlap(results[r], intervals_b[bp]):
            results[r:(r + 1)] = subtract(results[r], intervals_b[bp])
        else:
            r += 1
    return results
        
def count_spans(intervals):
    return reduce(lambda x, y: x + y.size(), intervals, 0)

def main():
    oparser = optparse.OptionParser('usage: %prog [options] first.intervals '
        'second.intervals')
    oparser.add_option('--union', action='store_true')
    oparser.add_option('--intersect', action='store_true')
    oparser.add_option('--subtract', action='store_true')
    oparser.add_option('--summarize_only', action='store_true')
    (options, args) = oparser.parse_args()
    
    if ((options.union and (options.intersect or options.subtract)) or
        (options.intersect and (options.union or options.subtract)) or
        (options.subtract and (options.intersect or options.union))):
        print >>sys.stderr, ('specify exactly one of --union, --intersect, or '
            'subtract')
        return 1
    
    intervals_file_a = args[0]
    intervals_file_b = args[1]
    
    intervals_a = load_intervals(intervals_file_a)
    intervals_b = load_intervals(intervals_file_b)
    
    merge_overlaps(intervals_a)
    merge_overlaps(intervals_b)
    
    if options.subtract:
        result = subtract_overlaps(intervals_a, intervals_b)
    elif options.union:
        result = intervals_a + intervals_b
        result.sort()
        merge_overlaps(result)
    elif options.intersect:
        intervals_c = intervals_a + intervals_b
        intervals_c.sort()
        result = intersect_overlaps(intervals_c)
    else:
        print >>sys.stderr, ('somehow we got to this point without a valid '
            'operator')
        return 1
    
    if options.summarize_only:
        print '{0} bases in {1}'.format(count_spans(intervals_a), 
            intervals_file_a)
        print '{0} bases in {1}'.format(count_spans(intervals_b),
            intervals_file_b)
        print ('{0} bases in result'.format
            (count_spans(result), intervals_file_a, intervals_file_b))
    else:
        print '\n'.join(map(lambda x: str(x), result))
    
    return 0
    
if __name__ == '__main__':
    exit(main())
