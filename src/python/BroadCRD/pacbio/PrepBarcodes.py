#!/usr/bin/env python

import string
import sys
from BroadCRD.util.fasta import readFasta, writeFasta

# This program, which might not be generally useful, preps barcodes from
# the format provided by the lab, into the orientation and naming convention
# expected by PacBio's PacBioBarcodeCCSDist software.

barcodes = readFasta(sys.argv[1])

for b in range(len(barcodes)):
    if barcodes[b][0][-1] == 'F':
        barcodes[b][0] = 'F{0}'.format(b / 2)
    else:
        barcodes[b][0] = 'R{0}'.format(b / 2)
        barcodes[b][1] = barcodes[b][1][::-1].translate(string.maketrans
            ('ACTGactg', 'TGACtgac'))

writeFasta(barcodes, sys.argv[2])

exit(0)
