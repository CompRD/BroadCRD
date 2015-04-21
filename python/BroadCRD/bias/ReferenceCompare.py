#!/usr/bin/env python

from BroadCRD.util.fasta import readFastaDict

ucsc = readFastaDict('hg19.fa')
broad = readFastaDict('/seq/references/Homo_sapiens_assembly19/v1/'
    'Homo_sapiens_assembly19.fasta')

ucsc_lens = {}
for c in ucsc:
    ucsc_lens[len(ucsc[c])] = c

unmatched = []
for c in broad:
    if (len(broad[c]) in ucsc_lens and
        broad[c].upper() == ucsc[ucsc_lens[len(broad[c])]].upper()):
        print '{0} {1}'.format(c, ucsc_lens[len(broad[c])])
    else:
        unmatched.append(c)

print '--UNMATCHED--'
for u in unmatched:
    print u
    
exit(0)
