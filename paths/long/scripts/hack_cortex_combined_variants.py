#!/usr/bin/env python

#
# WARNING: this is really ugly code with a really specific purpose -- it
# will turn certain combined subs/indel variants in a cortex output file
# into a subs/indel pair of variant lines.  This only works for lines
# with a single-base substitution combined with an indel, but this
# accounts for approximately 7100 of 7200 such combined events in the
# CORTEX file.
#
# The whole point of this was just to see if we were unnecessarily
# penalizing CORTEX.
#


import sys
import argparse
import math
import os
import copy
from VCF import *


def check_vcf( input ):
    v = VCF( input )

    debug = False

    print "\n".join(v.metadata)

    for line in v.lines():
        r = line.ref
        a = line.alt_list

        if len(r) > 1:
            if len(a) > 1:
                raise Exception("WARNING: multi-allelic change not coded for")

            line1 = copy.deepcopy( line )
            line2 = copy.deepcopy( line )

            if len(a[0]) > 1:
                if len(r) < len(a[0]):          # insertion
                    if len(r) == 2:
                        line1.ref = r[1]
                        line1.alt = r[1]+a[0][2:] 
                        line1.pos += 1
                        line2.ref = r[1]
                        line2.alt = a[0][1] 
                        line2.pos += 1
                    if debug: print "===="
                    print line1
                    print line2
                    if debug: print line
                elif len(a[0]) < len(r):
                    if len(a[0]) == 2:
                        line1.ref = r[1:]
                        line1.alt = r[1]
                        line1.pos += 1
                        line2.ref = r[1]
                        line2.alt = a[0][1]
                        line2.pos += 1
                    if debug: print "===="
                    print line1
                    print line2
                    if debug: print line
                else:
                    print line
            else:
                print line
        else:
            print line


def main( argv=[__name__] ):
    parser=argparse.ArgumentParser(description='rewrite VCF file in reference contig order')
    parser.add_argument( 'input_vcf', help='VCF input file' )
    args = parser.parse_args(argv[1:])

    return check_vcf( args.input_vcf )

if __name__ == "__main__":
    sys.exit(main(sys.argv))
