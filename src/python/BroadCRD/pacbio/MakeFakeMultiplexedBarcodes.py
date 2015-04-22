#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Generate new PacBio barcode combinations from a FASTA of unique front and rear
# barcodes, and similarly chop and mix barcoded amplicons to make a test data 
# set.

import random
import string
import sys

from BroadCRD.util.fasta import readFasta, writeFasta

def multiplex(front_halves, back_halves, multiplex_label):
    output_length = min(len(front_halves), len(back_halves))
    output = []
    for o in range(output_length):
        front = front_halves[o][1][0:(len(front_halves[o][1]) / 2)]
        back = back_halves[o][1][(len(back_halves[o][1]) / 2):]
        seq = front + back
        if random.random() > 0.5:
            seq = seq[::-1].translate(string.maketrans('ACTGactg', 'TGACtgac'))
        output.append((multiplex_label + '_{0}'.format(o), seq))
    return output

def split_barcodes(barcodes):
    front_barcodes = [None]*(len(barcodes) / 2)
    back_barcodes = [None]*(len(barcodes) / 2)
    
    for b in barcodes:
        if b[0][0] == 'F':
            front_barcodes[int(b[0][1:])] = b
        elif b[0][0] == 'R':
            back_barcodes[int(b[0][1:])] = b
        else:
            assert False

    return (front_barcodes, back_barcodes)

def main():
    barcode_fasta = sys.argv[1]
    amplicon_fastas = sys.argv[2:]
    
    barcodes = readFasta(barcode_fasta)
    
    assert len(barcodes) == len(amplicon_fastas) * 2
    (front_barcodes, back_barcodes) = split_barcodes(barcodes)
    
    amplicon_fastas.sort()
    print 'Processing, in order: ' + ', '.join(amplicon_fastas)
    amplicons = [readFasta(afasta) for afasta in amplicon_fastas]
    
    multiplex_barcodes = []
    with open('multiplex_decode.txt', 'w') as multiplex_decode:
        d = 0
        for a, afasta in enumerate(amplicons):
            for b in range(a, a + 12):
                bind = b % len(amplicons)
                writeFasta(multiplex(amplicons[a], amplicons[bind],
                    'F_{0}'.format(d)), 
                    'multiplex_amplicons_{0}.fasta'.format(d))
                multiplex_barcodes.append(('F{0}'.format(d), 
                    front_barcodes[a][1]))
                multiplex_barcodes.append(('R{0}'.format(d),
                    back_barcodes[bind][1]))
                print >>multiplex_decode, 'F{0} = {1} {2}'.format(d, a, bind)
                d += 1
    
    writeFasta(multiplex_barcodes, 'multiplex_barcodes.fasta')

    return 0

if __name__ == '__main__':
    exit(main())
