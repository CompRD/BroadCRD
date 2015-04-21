###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################

from Bio import SeqIO

# example usage:
#
#    print str(datetime.datetime.now()), "reading masks"
#    fDust="/home/unix/blau/wga/dustmasker/hg19/ofn"
#    dust_masks=locus_masks.dustmask(fDust)
#    print str(datetime.datetime.now()), "reading masks"
#    fSegDup="/home/unix/blau/wga/segdup/GRCh37GenomicSuperDup.tab"
#    segdup_masks=locus_masks.segdupmask(fSegDup)
#    print 'segdup masked bases: ', segdup_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
#    print 'DUST   masked bases: ', dust_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")
#    print str(datetime.datetime.now()), "done reading masks"
#    union_masks = locus_masks.unionmask([segdup_masks,dust_masks])
#    print 'UNION  masked bases: ', union_masks.nMasked("/wga/scr4/bigrefs/human19/genome.fasta")

class mask_base:
    def load_mask(self,filename):
        pass
    def __init__(self,filename):
        self.load_mask(filename)

    def nMasked(self,fName):
        count=0
        for record in SeqIO.parse( fName, "fasta"):
            loc_seq = str(record.seq)
            key =  record.id
            if key not in self.ranges: continue
            for interval in self.ranges[key]:
                sequence = str( loc_seq[max(0,interval[0]-1):interval[1]])
                count += len(sequence) - sequence.count('N') - sequence.count('n')
        return count

    def isMasked(self,chrom,front,back):
        if chrom not in self.ranges: return False
        loc_list = self.ranges[chrom]
        if len(loc_list)==0 : return False

        left_limit = self.find_left(loc_list,0,len(loc_list)-1,front)
        if ( left_limit >= 0 ) and ( front <= loc_list[left_limit][1] ):
            assert (loc_list[left_limit][0] <= front and front <= loc_list[left_limit][1] )
            return True

        left_limit = self.find_left(loc_list,0,len(loc_list)-1,back)
        if ( left_limit >= 0 ) and ( front <= loc_list[left_limit][1] ):
            assert (loc_list[left_limit][0] <= back and front <= loc_list[left_limit][1] )
            return True

        return False
    def find_left(self, list, lo, hi, val):
        if hi < lo or list[lo][0] > val : return -1
        if list[hi][0] <= val : return hi
        if list[lo][0] == val : return lo

        mid = (lo+hi)//2
        mid_val = list[mid][0]

        if mid_val == val:
            return mid
        elif mid_val < val:
            tmp = self.find_left(list,mid+1,hi,val)
            if tmp==-1: return mid
            else: return tmp
        else:
            return self.find_left(list,lo,mid,val)

class dustmask(mask_base):
    def name(self):
        return "DUST"
    def load_mask(self,filename):
        self.ranges={}
        with open(filename,'r') as fin:
            key=None
            for line in fin:
                if len(line) < 2: continue
                buffer = line.split()
                if buffer[0][0] == '>':
                    key = buffer[0][1:]
                    if key not in self.ranges:
                        self.ranges[key]=[]
                else :
                    assert( len(buffer)==3)
                    assert( key )
                    front = int(buffer[0])+1
                    back = int(buffer[2])
                    self.ranges[key].append((front,back))
        for key in self.ranges:
            self.ranges[key].sort()
            loc_list = self.ranges[key]
            nn=len(loc_list)
            if nn < 2 : continue
            for ii in range(1,nn):
                assert( loc_list[ii][0] > loc_list[ii-1][1])

class segdupmask(mask_base):
    def name(self):
        return "SEG_DUP"
    def load_mask(self,filename):
        self.ranges={}
        with open(filename,'r') as fin:
            for line in fin:
                if len(line) < 2: continue
                buffer = line.split()
                assert ( len(buffer) == 29 )
                if buffer[0] == 'chrom': continue

                chrom = buffer[0][3:]
                front = int(buffer[1])
                back  = int(buffer[2])
                if chrom not in self.ranges: self.ranges[chrom]=[]
                self.ranges[chrom].append((front,back))
        for key in self.ranges:
            self.ranges[key].sort()
            loc_list = self.ranges[key]
            nn=len(loc_list)
            if nn < 2 : continue

            new_list=[]

            cluster_front=-1
            cluster_back=-1
            for entry in loc_list:
                if entry[0] > cluster_back:
                    if cluster_front>=0:
                        new_list.append( (cluster_front,cluster_back) )
                    cluster_front = entry[0]

                cluster_back=max(cluster_back,entry[1])

            if cluster_front>=0:
                new_list.append( (cluster_front,cluster_back) )

            self.ranges[key]=new_list

            loc_list = self.ranges[key]
            nn=len(loc_list)
            if nn < 2 : continue
            for ii in range(1,nn):
                assert( loc_list[ii][0] > loc_list[ii-1][1])

class crd_mask(mask_base):
    def name(self):
        return self.mask_name
    def setName(self,new_name):
        self.mask_name=new_name
    def load_mask(self,filename):
        self.mask_name="CRD"

        self.ranges={}
        with open(filename,'r') as fin:
            for line in fin:
                if len(line) < 2: continue
                buffer = line.split(':')
                assert len(buffer)==2
                range_0_base = buffer[1].split('-')
                assert len(range_0_base)==2

                chrom = buffer[0]
                front = int(range_0_base[0])+1
                back  = int(range_0_base[1])

                if chrom not in self.ranges: self.ranges[chrom]=[]
                self.ranges[chrom].append((front,back))
        for key in self.ranges:
            self.ranges[key].sort()
            loc_list = self.ranges[key]
            nn=len(loc_list)
            if nn < 2 : continue

            new_list=[]

            cluster_front=-1
            cluster_back=-1
            for entry in loc_list:
                if entry[0] > cluster_back:
                    if cluster_front>=0:
                        new_list.append( (cluster_front,cluster_back) )
                    cluster_front = entry[0]

                cluster_back=max(cluster_back,entry[1])

            if cluster_front>=0:
                new_list.append( (cluster_front,cluster_back) )

            self.ranges[key]=new_list

            loc_list = self.ranges[key]
            nn=len(loc_list)
            if nn < 2 : continue
            for ii in range(1,nn):
                assert( loc_list[ii][0] > loc_list[ii-1][1])

class unionmask(mask_base):
    def name(self):
        return "UNION"
    def load_mask(self,masks):
        self.ranges={}
        for mask in masks:
            for chrom in mask.ranges:
                if chrom not in self.ranges: self.ranges[chrom]=[]
                self.ranges[chrom].extend(mask.ranges[chrom])
            
        for key in self.ranges:
            self.ranges[key].sort()
            loc_list = self.ranges[key]
            nn=len(loc_list)
            if nn < 2 : continue

            new_list=[]

            cluster_front=-1
            cluster_back=-1
            for entry in loc_list:
                if entry[0] > cluster_back:
                    if cluster_front>=0:
                        new_list.append( (cluster_front,cluster_back) )
                    cluster_front = entry[0]

                cluster_back=max(cluster_back,entry[1])

            if cluster_front>=0:
                new_list.append( (cluster_front,cluster_back) )

            self.ranges[key]=new_list

            loc_list = self.ranges[key]
            nn=len(loc_list)
            if nn < 2 : continue
            for ii in range(1,nn):
                assert( loc_list[ii][0] > loc_list[ii-1][1])

