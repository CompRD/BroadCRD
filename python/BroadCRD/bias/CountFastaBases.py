#!/usr/bin/env python

import sys

a = 0
c = 0
g = 0
t = 0

for fasta_filename in sys.argv[1:]:
    with open(fasta_filename, 'r') as fasta_file:
        print 'Processing ' + fasta_filename
        for line in fasta_file:
            if len(line) > 0 and line[0] != '>':
                for base in line:
                    if base == 'a' or base =='A':
                        a += 1
                    elif base == 'c' or base =='C':
                        c += 1
                    elif base == 'g' or base =='G':
                        g += 1
                    elif base == 't' or base =='T':
                        t += 1

known_bases = a + c + g + t

print('total known bases: {0}'.format(known_bases))
print('A: {0} ({1}%)'.format(a, 100 * a / known_bases))
print('C: {0} ({1}%)'.format(c, 100 * c / known_bases))
print('G: {0} ({1}%)'.format(g, 100 * g / known_bases))
print('T: {0} ({1}%)'.format(t, 100 * t / known_bases))

exit(0)

