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
    def __init__(self):
        self.nBarcode=0
        self.nUnpaired  =0
        self.nUnpairedMulti =0
        self.nUnpairedMultiChrom =0
        self.nAlignments=0

        #number of IVT alignment is sum_(barcode,chromosome) (num_fw_alignment -1 + num_rc_alignment-1)
        self.nIVT=0
        self.nMultiChrom=0

        self.nPositives=0
        self.nNegatives=0

        #number of IVT alignment is sum_(barcode,chromosome) abs(num_fw_alignment - num_rc_alignment)
        self.nAbsFwdRev=0

        self.singletDist={}

    def report(self):
        print 'number of alignments:                                                ', self.nAlignments
        print 'number of IVT alignments:                                            ', self.nIVT
        print 'abs of number of forward alignments minus reverse alignments:        ', self.nAbsFwdRev
        print
        print 'number of barcodes with alignments:                                  ', self.nBarcode
        print 'number of barcodes with unpaired alignments:                         ', self.nUnpaired
        print 'number of barcodes with unpaired alignments with multiple alignments:', self.nUnpairedMulti
        print 'number of barcodes with multi-chrom alignment:                       ', self.nMultiChrom
        print 'number of barcodes with unpaired, multi-chrom alignment:             ', self.nUnpairedMultiChrom
        print
        print 'number of positive jumps:                                            ', self.nPositives
        print 'number of negative jumps:                                            ', self.nNegatives

        print "statistics for unpaired barcode"
        print "per-chrom abs barcode"
        aa=0
        bb=0
        print "chromosome_index abs(num_barcode_with_fw_alignment - num_barcode_with_rc_alignment) (num_barcode_with_fw_alignment + num_barcode_with_rc_alignment)"
        for chrom in self.singletDist:
            aa+=abs(self.singletDist[chrom][0]-self.singletDist[chrom][1])
            bb+=abs(self.singletDist[chrom][0]+self.singletDist[chrom][1])
            print chrom, abs(self.singletDist[chrom][0]-self.singletDist[chrom][1]),self.singletDist[chrom][0]+self.singletDist[chrom][1]
        print "sum(chrom) abs(num_barcode_with_fw_alignment - num_barcode_with_rc_alignment) (num_barcode_with_fw_alignment + num_barcode_with_rc_alignment)",aa,bb
    def add_collection(self,collection):
        # a collection of data per barcode,
        # collection[chrom][0] is the (front,back) of a FW alignment of the gDNA to chromosome index 'chrom', according to Sante's code
        # collection[chrom][1] is the (front,back) of a RC alignment of the gDNA to chromosome index 'chrom', according to Sante's code
        if len(collection) == 0 : return

        self.nBarcode+=1

        bPaired=False

        locNAlign=0

        locChromDist={}

        # chromosome-by-chromosome
        for chrom in collection:
            # locN[0] is the number of FW alignment of this barcode to this chromosome
            # locN[1] is the number of RC alignment of this barcode to this chromosome
            locN=[0,0]

            #cluster up the alignment -- there can be overlaping alignment because of IVT duplciation
            #clusters[0] - fw
            #clusters[1] - rc

            clusters=[[],[]]
            for rev in range(2):
                ranges=collection[chrom][rev]
                ranges.sort()
                last_back=-1
                for entry in ranges:
                    self.nAlignments+=1
                    locN[rev]+=1
                    locNAlign+=1
                    front = entry[0]
                    back  = entry[1]
                    if front > last_back: clusters[rev].append([])
                    clusters[rev][-1].append((front,back));
                    last_back=max(last_back,back)
            if len(clusters[0])>0 and len(clusters[1])>0 :
                # if there are both fw and rc clusters, there is a jump
                # let's see if the jump of this barcode of this chromosome is positive or negative
                bPaired=True

                p_distances=[]
                n_distances=[]

                for c0 in clusters[0]:
                    for c1 in clusters[1]:
                        for r0 in c0:
                            for r1 in c1:
                                distance = r0[1] -r1[0]
                                if distance >0:
                                    p_distances.append(distance)
#                                    self.nPositives+=1
                                else:
                                    n_distances.append(distance)
#                                    self.nNegatives+=1
                p_distances.sort()
                n_distances.sort()

                #let's do a majority vote for each barcode and each chromosome instead of counting all jumps
                if len(p_distances) > len(n_distances):
                    self.nPositives+=1
                else:
                    self.nNegatives+=1

            self.nAbsFwdRev += abs(locN[0]-locN[1])

            if locN[0]>1: self.nIVT+= locN[0] -1
            if locN[1]>1: self.nIVT+= locN[1] -1

            if len(clusters[0])>0 or len(clusters[1])>0 :
                locChromDist[chrom]=locN
        if not bPaired:
            self.nUnpaired+=1
            if locNAlign> 1: self.nUnpairedMulti+=1
            if len(locChromDist)> 1:
                self.nUnpairedMultiChrom+=1
            elif len(locChromDist) == 1:
                chrom = locChromDist.keys()[0]
                if chrom not in self.singletDist: self.singletDist[chrom]=[0,0]
                assert( (locChromDist[chrom][0]==0) != (locChromDist[chrom][1]==0) )
                if (locChromDist[chrom][1]==0):
                    self.singletDist[chrom][0]+=1
                else:
                    self.singletDist[chrom][1]+=1




        if len(locChromDist) >1: self.nMultiChrom+=1


    
def analyze_sis(records):
    last_barcode=-1

    # a collection of data per barcode,
    # collection[chrom][0] is the (front,back) of a FW alignment of the gDNA to chromosome index 'chrom', according to Sante's code
    # collection[chrom][1] is the (front,back) of a RC alignment of the gDNA to chromosome index 'chrom', according to Sante's code
    collection={}

    #initiate analyzer
    analyzer=analyze_barcode()

    #line counter of sis file
    cc=0

    #fwrev[chrom][0] is the total number of gDNA (across all chromosome and all barcode) having FW alignment to that chromosome
    #fwrev[chrom][1] is the total number of gDNA (across all chromosome and all barcode) having RC alignment to that chromosome
    fwrev={}

    for line in records:
        buffer=line.split()
        # each line has 4 fields, 
        # barcode_index chromosome_index, front of alignment, back of alignment
        # if chromosome_index is negative, -chrom-1 is the chromosome index
        assert(len(buffer)==4)

        barcode=int(buffer[0])
        if barcode != last_barcode:
            #update analyzer with information collected with this barcode
            analyzer.add_collection(collection)
            collection={}
        last_barcode=barcode

        chrom=int(buffer[1])
        front=int(buffer[2])
        back =int(buffer[3])
        rev=0
        if chrom < 0 :
            rev=1
            chrom = -(chrom+1)
        assert(chrom>=0)

        if chrom not in fwrev: fwrev[chrom]=[0,0]
        fwrev[chrom][rev]+=1

        if chrom not in collection: collection[chrom]=[ [], [] ]
        collection[chrom][rev].append( (front,back) )
        cc+=1
        if cc%1000000==0 :
            sys.stdout.write( '.')
            sys.stdout.flush()

    aa=0
    aaa=0
    print "chromosome_index abs(num_fw_alignment - num_rc_alignment) (num_fw_alignment + num_rc_alignment)"
    for chrom in fwrev:
        aa+=abs(fwrev[chrom][0]-fwrev[chrom][1])
        aaa+=fwrev[chrom][0]+fwrev[chrom][1]
        print chrom, abs(fwrev[chrom][0]-fwrev[chrom][1]),fwrev[chrom][0]+fwrev[chrom][1]
    print "sum_chromosome abs(num_fw_alignment - num_rc_alignment) total", aa, aaa
    sys.stdout.write('\n')

    #update analyzer with information collected with the last barcode
    analyzer.add_collection(collection)
    analyzer.report()

        
    

if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='analyze sis file')
    parser.add_argument( 'sis_in', help='sis file generated by Sante\'s code' )
    args = parser.parse_args(sys.argv[1:])


    if not os.path.exists( args.sis_in ):
        print >>sys.stderr, "file ", args.sis_in, " does not exist."
        sys.exit(1)
    print str(datetime.datetime.now()), "analyzing"
    with open(args.sis_in,'r') as fin:
        analyze_sis(fin)
    print str(datetime.datetime.now()), "done"
    sys.exit(0)
