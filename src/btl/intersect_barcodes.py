#!/usr/bin/env python
###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################

import sys
import re
import copy
import argparse
import os
import math
import datetime

class analyze_barcode:
    def __init__(self,nBins):
        self.barcode_counts={}
        self.nBins=nBins

    def report(self):
        print "report"
        if self.nBins == 2 :
            n1=0
            n2=0
            n12=0

            for value in self.barcode_counts.itervalues():
                if value[0]>0 and value[1]>0:
                    n12+=1
                elif value[0]>0:
                    n1+=1
                elif value[1]>0:
                    n2+=1
                else:
                    print >>sys.stderr, "implementation error"
                    sys.exit(1)
            print "category    #barcodes:"
            print "file1-only", n1
            print "file2-only", n2
            print "both files", n12
        else:
            print "not implemented"

    def analyze(self,file_path,iBin):
        if not os.path.exists(file_path):
            print >>sys.stderr, "file ", file_path, " does not exist."
            sys.exit(1)
        print "analyzing ", file_path
        count=0
        with open(file_path,'r') as fin:
            for line in fin:
                if len(line)==0 or line[0]!='>': continue
                count+=1

                barcode = line.split('_')[1]

                self.barcode_counts.setdefault(barcode,[0]*self.nBins)[iBin]+=1
                if count%1000000==0:
                    sys.stderr.write('.')
        sys.stderr.write('\n')

if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='analyze sis file')
    parser.add_argument( 'barcode1', help='fasta file generated by Sante\'s code' )
    parser.add_argument( 'barcode2', help='fasta file generated by Sante\'s code' )
    args = parser.parse_args(sys.argv[1:])

    analysis = analyze_barcode(2)

    analysis.analyze(args.barcode1,0)
    analysis.analyze(args.barcode2,1)

    analysis.report()

    sys.exit(0)
