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
import zlib
import datetime
from Bio import SeqIO
import locus_masks

def zlib_ratio_DNA( data ):
    if ( not data ) or len(data)==0 : return 1
    return float(len(zlib.compress(data)))/len(data)

def Entropy_Shannon_DNA( data ):
    entropy=0.0
    if ( not data ) or len(data)==0 : return entropy
    loc_data=data.upper()
    length=float(len(data))
    for base in ('A','C','G','T'):
        p_x=float(data.count(base))/length
        if p_x>0:
          entropy -= p_x * math.log(p_x,2)
    return entropy



dna2bit = [-100000]*256
dna2bit[ord('A')]=0
dna2bit[ord('a')]=0
dna2bit[ord('C')]=1
dna2bit[ord('c')]=1
dna2bit[ord('G')]=2
dna2bit[ord('g')]=2
dna2bit[ord('T')]=3
dna2bit[ord('t')]=3

def Entropy_Shannon_DNA_2mer( data ):
    entropy=0.0
    if ( not data ) or len(data)<2 : return entropy

    n2mer = len(data)-1

    count2mer = [0]*16

    for offset in range(n2mer):
        value = dna2bit[ord(data[offset])]+ 4*dna2bit[ord(data[offset+1])]
        count2mer[value]+=1

    for count in count2mer:
        p_x=float(count)/float(n2mer)
        if p_x>0: entropy -= p_x * math.log(p_x,2)

    return entropy

def Entropy_Shannon_DNA_3mer( data ):
    entropy=0.0
    if ( not data ) or len(data)<3 : return entropy

    n3mer = len(data)-2

    count3mer = [0]*64

    for offset in range(n3mer):
        value = dna2bit[ord(data[offset])]+ 4*dna2bit[ord(data[offset+1])] + 16*dna2bit[ord(data[offset+2])]
        count3mer[value]+=1

    for count in count3mer:
        p_x=float(count)/float(n3mer)
        if p_x>0: entropy -= p_x * math.log(p_x,2)

    return entropy

def Entropy_Shannon_DNA_4mer( data ):
    entropy=0.0
    if ( not data ) or len(data)<4 : return entropy

    n4mer = len(data)-3

    count4mer = [0]*64*4

    for offset in range(n4mer):
        value =    dna2bit[ord(data[offset])] \
               + 4*dna2bit[ord(data[offset+1])] \
               + 16*dna2bit[ord(data[offset+2])] \
               + 64*dna2bit[ord(data[offset+3])] 
        count4mer[value]+=1

    for count in count4mer:
        p_x=float(count)/float(n4mer)
        if p_x>0: entropy -= p_x * math.log(p_x,2)

    return entropy



class MissedCallsReport:
    def __init__(self,caller,masks):
        self.caller=caller
        self.masks=masks
        self.clusters=[]

        self.nSingle=0
        self.nMultiple=0
        self.nFMasked=0
        self.nFUnMasked=0

        self.nFosmid=[0,0]
        self.nMiss=[0,0]
        self.nMissSNP=[0,0]
        self.nMissIND=[0,0]
        self.nMissCPX=[0,0]


        self.nMissMasked=[0,0]
        self.nMissUnMasked=[0,0]

        self.compMiss=[]
        self.compMade=[]

        self.f_comp=re.compile("Fosmid")
        self.c_comp=re.compile(caller)
        self.a_comp=re.compile("(Fosmid|"+caller+")")

    def loadClusters(self,clusters):
        for cluster in clusters: self.processCluster(cluster)

    def processCluster(self,cluster):
        assert( len(cluster) != 0 )
        bin_idx = -1
        if len(cluster)==1:
            self.nSingle+=1
            bin_idx=0
        else:
            self.nMultiple+=1
            bin_idx=1

        for line in cluster:
            chrom_idx=0
            front_idx=1
            back_idx=2
            ref_idx=3
            alt_idx=4
            info_idx=5
            comp_idx=6

            if self.f_comp.search( line[info_idx] ) is None: continue

            bIsMasked=False
            for mask in self.masks:
                bIsMasked = bIsMasked or mask.isMasked(line[chrom_idx],line[front_idx],line[back_idx])

#            bIsMasked = self.masks.isMasked(line[chrom_idx],line[front_idx],line[back_idx])
            if bIsMasked:
                self.nFMasked+=1
            else:
                self.nFUnMasked+=1

            self.nFosmid[bin_idx]+=1
            complexity=line[comp_idx]
            if self.c_comp.search( line[info_idx] ) is None:
                self.nMiss[bin_idx]+=1
                self.compMiss.append(complexity)

                if bIsMasked:
                    self.nMissMasked[bin_idx]+=1
                else:
                    self.nMissUnMasked[bin_idx]+=1

                if( len(line[ref_idx]) != len(line[alt_idx])):
                    self.nMissIND[bin_idx]+=1
                elif len(line[ref_idx]) == 1:
                    self.nMissSNP[bin_idx]+=1
                else:
                    self.nMissCPX[bin_idx]+=1
            else:
                self.compMade.append(complexity)

    def report(self,out_dir=None):
        print
        print "Caller:                             ",self.caller#,self.nSingle,self.nMultiple
        print "fosmid calls:                       ",self.nFosmid
        print "Missed calls:                       ",self.nMiss
        print "Missed SNP:                         ",self.nMissSNP
        print "Missed IND:                         ",self.nMissIND
        print "Missed CPX:                         ",self.nMissCPX
        print
        print "fosmid calls(masked,unmasked):      ",self.nFMasked,self.nFUnMasked
        print "Missed Masked:                      ",sum(self.nMissMasked)
        print "fraction missed in masked fosmid:   ",float(sum(self.nMissMasked))/self.nFMasked
        print "Missed UnMasked:                    ",sum(self.nMissUnMasked)
        print "fraction missed in unmasked fosmid: ",float(sum(self.nMissUnMasked))/self.nFUnMasked
        print

        if out_dir:
            fCompMade=out_dir+"/"+self.caller+".compMade"
            fCompMiss=out_dir+"/"+self.caller+".compMiss"
#            print "writing to: ", fCompMade, fCompMiss
            with open(fCompMade,'w') as fout:
                for entry in self.compMade: fout.write("{0:e}\n".format(entry))
            with open(fCompMiss,'w') as fout:
                for entry in self.compMiss: fout.write("{0:e}\n".format(entry))
              


def parse_vall_line( line ):
    fid_idx=0
    flank_idx=1
    loc_pos_idx=2
    pos_idx=3
    ref_idx=4
    alt_idx=5
    info_idx=6
    buffer=line.split()

    (chrom,front)=buffer[pos_idx].split(':')
    front = int(front)
    ref=buffer[ref_idx]
    back = front+len(ref)-1

    alt=buffer[alt_idx]
    info="".join(buffer[info_idx:])
    flanks = buffer[flank_idx].split('/')
    assert(len(flanks)==2)
#    if len(ref)==1:
#        region=flanks[0]+alt+flanks[1]
#    else:
#        shift = len(ref)-1
#        region=flanks[0]+alt+flanks[1][shift:]
    region=flanks[0]+alt+flanks[1]
#    region=alt.join( buffer[flank_idx].split('/') )
#    return [chrom,front,back,ref,alt,info,zlib_ratio_DNA(region)]
    return [chrom,front,back,ref,alt,info,Entropy_Shannon_DNA_2mer(region),flanks]

def vall_cluster(vall_log,comp):
    vall_log_loc=sorted(vall_log)
    clusters=[]
    last_back=-1
    last_chrom=""

    chrom_idx=0
    front_idx=1
    back_idx=2
    info_idx=5

    for entry in vall_log_loc:
        if comp.search( entry[info_idx] ) is None: continue

        chrom = entry[chrom_idx]
        front = entry[front_idx]
        back  = entry[back_idx]
        if chrom != last_chrom or front > last_back :
            clusters.append([])
        clusters[-1].append(entry)
        last_chrom=chrom
        last_back = max(last_back,back)
    return clusters

def variant_tag(entry):
    r_idx=3
    f_idx=7
    ref=entry[r_idx]
    flanks=entry[f_idx]
    ref_seq = flanks[0] + ref + flanks[1]
    max_homo_len=0

    start=0
    for ii in range(1,len(ref_seq)):
        if ref_seq[ii]!=ref_seq[ii-1]:
            homo_len= ii-1-start+1
            max_homo_len = max(max_homo_len,homo_len)
            start = ii
    max_homo_len = max(max_homo_len,len(ref_seq)-1-start+1)
    if(max_homo_len>=5):
        return "HOMO5+"
    else:
        return ""


def analyze_fosmid(fdir,reports,out_vall=None):
    vall_file = fdir+"/variants.all"
    vall_log = []
    if not os.path.exists(vall_file): return
    with open(vall_file, 'r') as fvall:
        for line in fvall:
            vall_log.append( parse_vall_line(line))
            outstr=line[:-1]
            if out_vall is not None:
                for mask in reports[0].masks:
                    if mask.isMasked( vall_log[-1][0], vall_log[-1][1], vall_log[-1][2] ):
                        outstr+=","+mask.name()
                tmp=variant_tag(vall_log[-1])
                if len(tmp)>0:
                    outstr+=","+tmp
                out_vall.write(outstr+'\n')

    for report in reports:
        report.loadClusters( vall_cluster(vall_log,report.a_comp ) )

def analyze_all( root_dir , reports , out_vall=None, max_fid=106 ):
    for fid in range(max_fid+1):
        fdir = root_dir+"/"+str(fid)
        if not os.path.exists(fdir): continue
        analyze_fosmid(fdir , reports , out_vall)
    

if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='categorize missed variant calls')
    parser.add_argument( 'cv_dir', help='CompareVars directory' )
    parser.add_argument( 'out_dir', help='output directory' )
    args = parser.parse_args(sys.argv[1:])


    if not os.path.exists( args.cv_dir ):
        print >>sys.stderr, "CompareVars directory ", args.cv_dir, " does not exist."
        sys.exit(1)

    if os.path.exists( args.out_dir ):
        print >>sys.stderr, "Output directory ", args.out_dir, " already exists."
        sys.exit(1)
    print str(datetime.datetime.now()), "reading masks"
    fDust="/home/unix/blau/wga/dustmasker/hg19/ofn"
    dust_masks=locus_masks.dustmask(fDust)
    print str(datetime.datetime.now()), "reading masks"
    fSegDup="/home/unix/blau/wga/segdup/GRCh37GenomicSuperDup.tab"
    segdup_masks=locus_masks.segdupmask(fSegDup)
#    print 'segdup masked bases: ', segdup_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
#    print 'DUST   masked bases: ', dust_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
    print str(datetime.datetime.now()), "done reading masks"
#    union_masks = locus_masks.unionmask([segdup_masks,dust_masks])
#    print 'UNION  masked bases: ', union_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")

    callers = ("GATK-100","GATK-250","CORTEX","DISCOVAR")
    reports = [ MissedCallsReport(caller,[dust_masks,segdup_masks]) for caller in callers ]

    print str(datetime.datetime.now()), "analyzing"
    os.makedirs(args.out_dir)
    with open(args.out_dir+'/variant.all.marked','w') as vall_out:
        analyze_all(args.cv_dir,reports,vall_out)
    print str(datetime.datetime.now()), "done analyzing"

    for report in reports:
        report.report(args.out_dir)

    sys.exit(0)
