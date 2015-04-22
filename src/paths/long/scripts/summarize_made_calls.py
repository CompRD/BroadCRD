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
import locus_masks

class MadeCallsReport:
    def __init__(self,masks):
        self.masks=masks
        self.nMasks=len(self.masks)

        self.snp=0
        self.ids=0

        self.unmasked_snp=0
        self.unmasked_ids=0

        self.masked_snp=[0]*len(masks)
        self.masked_ids=[0]*len(masks)
        self.nLines=0;
    def report(self):
        print 'total number of snp:    ', self.snp
        print 'total number of ids:    ', self.ids
        print 'number of unmasked snp: ', self.unmasked_snp
        print 'number of unmasked ids: ', self.unmasked_ids
        print 'masks:                  ', [ i for i in map(lambda f: f.name(),self.masks) ]
        print 'masked snp:             ', self.masked_snp
        print 'masked ids:             ', self.masked_ids
        print 'number of lines:        ', self.nLines


    def analyze(self,vcf_file):
        c_idx=0;
        p_idx=1;
        r_idx=3;
        a_idx=4;
        f_idx=6;
        g_idx=9;
        with open(vcf_file,'r') as fin:
            for line in fin:
                if self.nLines % 100000 ==0: self.report()
                self.nLines+=1
                if len(line)>0 and line[0] != '#':
                    buffer=line.split();
                    assert len(buffer)==10

                    if buffer[c_idx]=='.': continue;
                    assert buffer[c_idx][:3] == 'chr'
                    chrom = buffer[c_idx][3:]

                    pos = buffer[p_idx]
                    ref = buffer[r_idx]
                    alt = buffer[a_idx]
                    fil = buffer[f_idx]
                    gt  = buffer[g_idx]

                    if pos=='.' or fil!='PASS' or ref=='.' or alt=='.' or gt=='.': continue

                    pos=int(pos)

                    gt = gt.split(':')[0]
                    if gt=='.': continue

                    gt = gt.replace('|','/').split('/')
                    if ( (gt[0]=='0'or gt[0]=='.')  and (gt[1]=='0' or gt[1]=='.')): continue;


                    if len(ref)==len(alt):
                        if ref != alt:
                            self.snp+=len(ref)
                            bMasked=False
                            for ii in range(self.nMasks):
                                if self.masks[ii].isMasked(chrom,pos,pos+len(ref)-1):
                                    self.masked_snp[ii]+=len(ref)
                                    bMasked=True
                            if not bMasked: self.unmasked_snp+=len(ref)
                    else:
                        self.ids+=1
                        bMasked=False
                        for ii in range(self.nMasks):
                            if self.masks[ii].isMasked(chrom,pos,pos+len(ref)-1):
                                self.masked_ids[ii]+=1
                                bMasked=True
                        if not bMasked: self.unmasked_ids+=1


if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='categorize missed variant calls')
    parser.add_argument( 'vcf_file', help='VCF file' )
    args = parser.parse_args(sys.argv[1:])


    if not os.path.exists( args.vcf_file ):
        print >>sys.stderr, "CompareVars directory ", args.vcf_file, " does not exist."
        sys.exit(1)

    print str(datetime.datetime.now()), "reading dust masks"
    fDust="/wga/scr4/blau/dustmasker/hg19/ofn"
    dust_masks=locus_masks.dustmask(fDust)

    print str(datetime.datetime.now()), "reading seg dup masks"
    fSegDup="/wga/scr4/blau/segdup/GRCh37GenomicSuperDup.tab"
    segdup_masks=locus_masks.segdupmask(fSegDup)

    print str(datetime.datetime.now()), "reading gc-85 masks"
    fGC85="/wga/scr4/blau/mross/gc85.intervals"
    gc85_masks=locus_masks.crd_mask(fGC85)
    gc85_masks.setName('GC85')

    print str(datetime.datetime.now()), "reading gc-75 masks"
    fGC75="/wga/scr4/blau/mross/gc75.intervals"
    gc75_masks=locus_masks.crd_mask(fGC75)
    gc75_masks.setName('GC75')

#    print 'segdup masked bases: ', segdup_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
#    print 'DUST   masked bases: ', dust_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
#    print 'gc85   masked bases: ', gc85_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
#    print 'gc75   masked bases: ', gc75_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
#    print str(datetime.datetime.now()), "done reading masks"
#
#    union_masks = locus_masks.unionmask([segdup_masks,dust_masks])
#    print 'UNION  masked bases: ', union_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")

    report = MadeCallsReport([dust_masks,segdup_masks,gc85_masks,gc75_masks])

    print str(datetime.datetime.now()), "analyzing"

    report.analyze(args.vcf_file);

    print 'final report:'
    report.report();

    sys.exit(0)
