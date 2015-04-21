#!/usr/bin/env python
###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################

import argparse
import sys
import re
import os

def find_region(dir, chr, pos):
    chr_dir=os.path.join(dir, chr)
    for root, dirs, files in os.walk( chr_dir, topdown=True ):
        if root != chr_dir: break
        for dir in dirs:
            match=re.match('(\d+)-(\d+)', dir)
            if match:
                start=int(match.groups()[0])
                end=int(match.groups()[1])
                if pos >= start and pos <= end:
                    print "{}/{}".format(root,dir)



################################################################################
#
#    M A I N
#
################################################################################
def main( argv = [__name__] ):
   parser=argparse.ArgumentParser(description="take a coordinate as chromosome and position and find the corresponding region dir")
   parser.add_argument('chromosome', nargs=1)
   parser.add_argument('position', nargs=1, type=int)
   parser.add_argument('directory', help='top level, above the chromosome', nargs=1)
   args=parser.parse_args(argv[1:])

   return find_region( args.directory[0], args.chromosome[0], args.position[0] )


if __name__ == '__main__':
    sys.exit(main(sys.argv))
