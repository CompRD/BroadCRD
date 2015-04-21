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
import tempfile
import fnmatch

######################################
# Exceptions
######################################

class NotDiscovarVCFException(Exception):
    def __init__(self, file, msg):
        newmsg="file " + file + " is not a Discovar-generated VCF: "+msg
        Exception.__init__(self,newmsg)

class BadMetaDataException(Exception):
    def __init__(self, filename, expected, found):
        sys.stderr.write("In file {}\n".format(filename))
        sys.stderr.write("Expected the following metadata:\n")
        for e in expected: sys.stderr.write(e)
        sys.stderr.write("but FOUND the following metadata:\n")
        for f in found: sys.stderr.write(f)
        Exception.__init__(self, """
            Short story: you are mixing filtered and unfiltered VCFs or
            VCFs from different Discovar versions.  We do not have code
            to merge metadata, so the files being merged must have the
            same metadata headers.""")

class BadInputException(Exception):
    def __init__(self,msg):
        Exception.__init__(self,msg)


###################################################################################
# TempFileDeleteOnDestroy - named temporary file that deletes on garbage collection,
# rather than file close().  The file is closed immediately.
#
# built on top of tempfile.NamedTemporaryFile, takes the same arguments
# and forwards the name property.  The delete= arg still determines
# whether the file is deleted, but the *when* has changed.  It now happens on
# object destruction.
###################################################################################
class TempFileDeleteOnDestroy:
    def __init__(self, **kwargs):

        # default "delete" arg is True
        if not kwargs.has_key('delete'): self._delete = True
        else: self._delete = kwargs['delete']

        # forward False to NamedTemporaryFile because we'll handle it
        kwargs['delete']=False

        self._tempfile = tempfile.NamedTemporaryFile(**kwargs)
        self._tempfile.close()

    def __del__(self):
        if self._delete:
            self._tempfile.unlink(self._tempfile.name)

    @property
    def name(self):
        return self._tempfile.name


###################################################################################
# Region - integer coordinates tuple of (chr,start,end)
#
# This is here to encapsulate all of the ugliness associated with chromosome naming
#
# internally it's a tuple with numeric chromosomes
#   X=23, Y=24
# constructs and prints with "chr:start-end"
# probably should switch to using a string for "chr" and have it be an assembly
# contig name
###################################################################################
class Region(tuple):
    """Genome coordinates internally a tuple of (chr, start, end), but
    externally it's chr:start-end."""
    # this is a tuple, and therefore immutable, so we need to reimplement __new__
    def __new__(self,str):
        """construct a Region from the string chr:start-end"""
        m=re.match('([^:]+):([0-9]+)-([0-9]+)',str)
        if not m or len(m.groups()) != 3:
            raise self.RegionException("error parsing region: "+str)

        (chr, start, end) = m.groups()

        istart=int(start)
        iend=int(end)

        return super(Region,self).__new__(self,tuple((chr,istart,iend)))

    def __str__(self):
        return "{}:{}-{}".format(*self)

    def __add__(self,y):
        """Provides the union of two regions if the chromosome is the
        same.  Throws BinaryMathRegionException otherwise."""
        (self_chr,self_start,self_end) = self
        (y_chr,y_start,y_end) = y
        if self_chr != y_chr:
            raise self.BinaryMathRegionException(self,y,"not on same chromosome")
        start = min(self_start, y_start)
        end = max(self_end,y_end)
        return Region("{}:{}-{}".format(self_chr, start, end))

    @property
    def chr(self): return self[0]

    @property
    def start(self): return self[1]

    @property
    def end(self): return self[2]

    class RegionException(Exception):
        def __init__(self,msg):
            Exception.__init__(self, msg)

    class BinaryMathRegionException(Exception):
        def __init__(self,x,y,msg):
            Exception.__init__(self, "can't do math on regions {} and {}: {}".format(x,y,msg))


################################################################################
# find_discovar_regions - find DiscovarRegion lines in a file or throw exception
#
# looks for vcf file lines of the form (including ##):
#
#      ##DiscovarRegion=<CHR="X",START=32380000,END=32410000>
#
# return a list of Regions
################################################################################
def find_discovar_regions(input):
    target='##DiscovarRegion'
    regions=[]
    for s in open(input,'r'):
        if len(s) >= len(target) and s[:len(target)] == target:
            m=re.match('##DiscovarRegion=<CHR="([^"]+)",START=([0-9]+),END=([0-9]+)', s)
            if not m: raise NotDiscovarVCFException(input,'error parsing '+s)
            regions.append(Region("{}:{}-{}".format(*m.groups())))
        elif s.split('\t')[0] == '#CHROM':      # for efficiency, skip the data section of the file
            break


    if not regions:
        raise NotDiscovarVCFException(input, "no "+target+" line found")

    return regions

################################################################################
# grab_metadata - pull all metadata lines from an input VCF file
################################################################################
def grab_metadata(input):
    metalines=[]
    samplelines=[]
    end_tag='#CHROM'
    skip_tag='##DiscovarRegion'
    for s in open(input,'r'):
        tag=s.split('\t')[0].split('=')[0]
        if tag == end_tag: break
        elif tag == skip_tag: pass
        elif tag == '##SAMPLE': samplelines.append(s)
        else: metalines.append(s)

    return ( metalines, samplelines )


################################################################################
# scan_inputs - call find_discovar_regions for each file and sanity check and
# sort Regions
################################################################################
def scan_inputs(inputs):

    # find out the union of the samples available in all VCFs
    file_regions=[]
    sample_match=[True]
    (metadatazero,samplezero) = ( set(x) for x in grab_metadata(inputs[0]) )

    for input in inputs:
        # sanity check the metadata
        # we assume that it is EXACTLY the same, with the exception of the
        # regions.  We now tolerate differences in samples, but only to
        # exclude regions with an incomplete number of samples.  It's up
        # to the user to decide if that's gone totally awry.
        if input != inputs[0]:
            # first check the metadata for an exact match
            (metadata,sampledata)=( set(x) for x in grab_metadata(input) )
            if metadatazero != metadata:
                raise BadMetaDataException(input, metadatazero, metadata)
            # then we check to see if this VCF introduces a novel sample
            before=len(samplezero)
            samplezero.update( sampledata )     # add any new sample names
            after=len(samplezero)
            # and now we check to see if this VCF has a full complement
            # of samples
            if after > before:
                sample_match = [ False for x in sample_match ]

            # either way, judge the current sample collection
            # we may have introduced a new sample, but still be deficient
            sample_match.append( sampledata == samplezero )

        regions=find_discovar_regions(input)
        if len(regions) != 1: raise BadInputException('input file cannot contain more than one region: ' +input)
        file_regions.append((input,regions[0]))

    # toss out any regions without a full complement of samples and
    # print a warning message.
    clean_file_regions=[]
    for match, region in zip( sample_match, file_regions ):
        if match:
            clean_file_regions.append( region )
        else:
            print >>sys.stderr, "WARNING: VCF {} does not include a full complement of samples...excluding it.".format( region[0] )


    return sorted(clean_file_regions, key=lambda s: s[1])      # sort by 2nd member - region


def chr_compare(x,y):
    # find common prefix of x,y
    i=0
    while i < len(x) and i < len(y) and x[i] == y[i]:
        i += 1
    if len(x[i:]) > 0 and len(y[i:]) > 0:
        x=x[i:]
        y=y[i:]
    # try to sort numerically
    try:
        xi=int(x)
        yi=int(y)
        return xi-yi
    except:
        if x < y:
            return -1
        elif x == y:
            return 0
        else:
            return 1

    assert(0);  # never reached
    return 0;

################################################################################
# merge_vcfs - subject a list of VCF files to merging
################################################################################
def merge_vcfs(inputs,output,verbose=False,hierarchical=True):

    outdir=os.path.dirname(output)
    if not outdir: outdir='.'

    file_regions=scan_inputs(inputs)
    print "{} regions found in VCF files:".format(len(file_regions))
    if verbose:
        for r in file_regions: print r

    # make a per-chromosome list of file_regions in chr_file_regions
    chr_file_regions={}
    for r in file_regions:
        if not chr_file_regions.has_key(r[1].chr):
            chr_file_regions[r[1].chr] = list()
        chr_file_regions[r[1].chr].append(r)


    # work out each chromosome
    temp_files_d={}
    ndeleted = 0
    for chr in chr_file_regions:
        print "processing chromosome {}".format(chr)
        tempout=TempFileDeleteOnDestroy(dir=outdir,delete=True,prefix='tmp.chr{}.'.format(chr))

        if hierarchical: ndeleted += hierarchical_merge_vcfs(chr_file_regions[chr],tempout.name,verbose)
        else: ndeleted += linear_merge_vcfs(chr_file_regions[chr],tempout.name)

        temp_files_d[chr]=tempout

    # append chromosomes to output
    #
    # pull out list of files sorted by chromosome. Because we don't have
    # access to the reference here, we try to be smart about chromosome
    # names by removing common prefixes and sorting numerically, when
    # possible.
    temp_files = [ temp_files_d[k] for k in sorted(temp_files_d.keys(),cmp=chr_compare) ]
    bl=copy_without( open(temp_files[0].name,'r'), open(output,'w') )
    for tf in temp_files[1:]:
        bl=copy_without( open(tf.name), open(output,'a'), '#', bl )


    print "Discarded {} conflicting or redundant entries".format(ndeleted)

################################################################################
# hierarchical_merge_vcfs - merge pairs, then pairs of pairs, etc.
################################################################################
def hierarchical_merge_vcfs(file_regions,output,verbose=False):
    outdir=os.path.dirname(output)
    if not outdir: outdir='.'

    new_temp_files=[]
    temp_files=None
    ndeleted=0
    while len(file_regions) > 1:
        print "merging {} regions".format(len(file_regions))
        sys.stdout.flush()
        new_file_regions=[]

        while len(file_regions) > 1:
            fr2 = file_regions.pop()
            fr1 = file_regions.pop()

            ( fr1_f, fr1_r ) = fr1
            ( fr2_f, fr2_r ) = fr2

            new_r = fr1_r + fr2_r

            temp=TempFileDeleteOnDestroy(dir=outdir,delete=True,prefix='tmp.{}.'.format(new_r))

            ndeleted += merge_vcf_avoid_block( fr1, fr2, temp.name, verbose )

            new_temp_files.append(temp)
            new_file_regions.append( ( temp.name, new_r ) )

            if verbose: print "{} + {} -> {}".format(fr1_r, fr2_r, new_r )
        else:
            if len( file_regions ) > 0:
                new_file_regions.extend(file_regions)

        file_regions = sorted(new_file_regions, key=lambda s: s[1])

    # copy the final file_regions entry to the output
    open( output, 'w' ).writelines( open( file_regions[0][0], 'r').readlines() )

    # all new_temp_files with be deleted upon garbage colleciton
    return ndeleted


################################################################################
# linear_merge_vcfs - merge a and b, then c with (a,b), then d with (a,b,c)
################################################################################
def linear_merge_vcfs(file_regions,output,verbose=False):
    outdir=os.path.dirname(output)
    if not outdir: outdir='.'
    ndeleted=0
    tempf1 = TempFileDeleteOnDestroy(dir=outdir, delete=False)
    tempf2 = TempFileDeleteOnDestroy(dir=outdir, delete=False)
    tempf1_name = tempf1.name
    tempf2_name = tempf2.name
    last_fr = ( tempf1_name, None )

    for next_fr in file_regions:
        ndeleted += merge_vcf_avoid_block( last_fr, next_fr, tempf2_name, verbose )
        last_fr = ( tempf2_name, next_fr[1] )
        # swap temp2 and temp1 names
        tempf1_name,tempf2_name = tempf2_name,tempf1_name

    os.rename( tempf1_name, output )
    os.unlink( tempf2_name )

    return ndeleted


################################################################################
# append_blocks - pulls blocks and micro_blocks out of a VCF file
#
# a block is a variant block represented by a unique BL= tag numbers
# a microblock is the (reference) coordinate range covered by a single
#       line of the VCF file.
################################################################################
def append_blocks( vcf_file, blocks, micro_blocks, verbose = False):
    if verbose: print "processing VCF file {}".format(vcf_file)
    with open(vcf_file,'r') as f:
        start=-1
        end=-1

        block_by_number={}

        header=True
        lastpos = None

        for line in f:
            if header and line[0] == '#': continue

            header=False

            entries=line.strip('\n').split()
            if len(entries) < 8: raise BadInputException("Bad VCF data line? "+line)

            # unpack the line
            (chr,cpos,id,ref,alt,qual,filter,info) = entries[:8]
            pos = int(cpos)

            assert( lastpos <= pos )
            lastpos = pos

            # unpack info line into a dictionary
            info_d = info_to_dict(info)

            try: this_bl = int( info_d['BL'] )
            except: raise BadInputException('no BL tag in info field: ' + info)

            end=max(pos+len(ref)-1,end)
            micro_blocks.append( (pos,end) )        # every line is a micro_block

            if block_by_number.has_key(this_bl):    # if the block exists, we just extend the end
                (olds, olde) = block_by_number[this_bl]
                block_by_number[this_bl] = (olds, end)
            else:                                   # new block
                block_by_number[this_bl] = (pos, end)

    blocks.extend( block_by_number.values() )



################################################################################
# divider_from_blocks - pick closest valid partition coordinate within a
# block set.
################################################################################
def divider_from_blocks( blocks, divider_in ):
    divider=divider_in
    blocks.sort()
    merged_blocks=[ [blocks[0][0],blocks[0][1]] ];
    for b in blocks:
        if b[0] > merged_blocks[-1][1]+1 :
            merged_blocks.append( [ b[0] , b[1] ])
        elif b[1] > merged_blocks[-1][1]:
            merged_blocks[-1][1] = b[1]
    offending_block = -1
    for mb_idx in range( len(merged_blocks) ):
        if ( merged_blocks[mb_idx][0] <= divider and divider <= merged_blocks[mb_idx][1] ):
            offending_block = mb_idx
    if( offending_block >= 0):
        if ( divider-merged_blocks[offending_block][0] < merged_blocks[offending_block][1]-divider ):
            divider = merged_blocks[offending_block][0] - 1
        else:
            divider = merged_blocks[offending_block][1] + 1
    return divider


################################################################################
# info_to_dict - breaks a VCF info line into a dictionary of elements
################################################################################
def info_to_dict( info ):
    # abuse of python -- this turns an INFO line like this:
    # BL=2;TR=5;SF;KL=6
    # into a dictionary of { BL:2, TR:5, SF:None, KL:6 }
    #
    info_els = info.split(';')
    # split into pairs around equals sign and cast to dictionary
    # handle the case where there's no equals sign
    info_d = dict( [ x.split('=') if '=' in x else (x, None) for x in info_els ] )

    return info_d


################################################################################
# dict_to_info - takes a  dictionary of elements (presumably from
# info_to_dict, above) and makes a VCF "info" line.
################################################################################
def dict_to_info( info_d ):
    #
    # again, abuse of python -- this is the opposite of info_to_dict (see above)
    #
    info_els = [ "=".join(x) if x[1] else str(x[0]) for x in info_d.items() ]
    info = ";".join(info_els)
    return info


################################################################################
# block_line - return the BL= block number associated with a line
################################################################################
def block_line( line ):
    entries=line.strip('\n').split('\t')
    if len(entries) < 8: raise BadInputException("Bad VCF data line? "+line)

    # unpack the line
    (chr,pos,id,ref,alt,qual,filter,info) = entries[:8]

    # find the BL tag
    info_d = info_to_dict( info )

    try: bl = int( info_d['BL'] )
    except: raise BadInputException('no BL tag in info field: ' + info)

    return bl



################################################################################
# block_modify - change the BL= block number associated with a line
################################################################################
def block_modify( line, blockno ):
    entries = line.strip('\n').split('\t')
    if len(entries) < 8: raise BadInputException("Bad VCF data line? "+line)

    # unpack the line
    info = entries[7]

    info_d = info_to_dict( info )
    if info_d.has_key('BL'): info_d['BL'] = str(blockno)
    info = dict_to_info( info_d )

    entries[7] = info

    return "\t".join(entries)+"\n"




################################################################################
# copy_without - copy from one file object to another file object
# and skip lines starting with a certain string
################################################################################
def copy_without( infd, outfd, skip = None, next_blockno = 1 ):

    blockno_to_new = {}

    for line in infd:
        if not skip or len(line) < len(skip) or line[:len(skip)] != skip:
            if line[0] == '#':
                outfd.write(line)
            else:
                this_block = block_line( line )
                if not blockno_to_new.has_key( this_block ):
                    blockno_to_new[ this_block ] = next_blockno
                    next_blockno += 1
                outfd.write(block_modify(line, blockno_to_new[this_block]))

    return next_blockno        # always set up for a new block

################################################################################
# merge_vcf_avoid_block - workhorse function to merge the VCF properly
#
# We try to merge between variant blocks but, if necessary, will divide
# variant blocks.
################################################################################
def merge_vcf_avoid_block(lo,hi,out,verbose=False):
    if verbose: print "\nmerging {}+{}->{}\n".format(lo,hi,out)

    strip_tag = '##DiscovarRegion'

    nDeleted=0;
    lo_f, lo_r = lo         # files, regions
    hi_f, hi_r = hi

    #
    # check if this is disjoint and do the trivial merge if so
    #

    if lo_r == None:
        if verbose: print "single file; trivial copy"
        copy_without(open(hi_f,'r'), open(out, 'w'), strip_tag)
        return 0
    elif lo_r.chr != hi_r.chr or lo_r.end <= hi_r.start:
        # no overlap between VCF files -- straight merge
        if verbose: print "no overlap; trivial merge"
        bl=copy_without(open(lo_f,'r'), open(out, 'w'), strip_tag)     # strip only ##DiscovarRegion
        copy_without(open(hi_f,'r'), open(out, 'a'), '#', bl)          # strip all lines starting with #
        return 0

    #
    # formal merge
    #

    blocks=[];
    micro_blocks=[];
    assert( lo_r.end > hi_r.start )
    divider = (lo_r.end + hi_r.start) // 2              # start off in the middle of the overlap

    append_blocks( lo_f, blocks, micro_blocks, verbose )
    append_blocks( hi_f, blocks, micro_blocks, verbose )


    if len(blocks) > 0:
        divider = divider_from_blocks(blocks,divider)
        if( lo_r.start > divider or divider > lo_r.end or hi_r.start > divider or divider > hi_r.end):
            print("for "+lo_f+" and "+hi_f)
            print("\toverlapping region fully spanned by blocks")
            print("\tseeking divider in micro blocks")
            divider = (lo_r.end + hi_r.start)//2
            divider = divider_from_blocks(micro_blocks,divider)


        if( lo_r.start > divider or divider > lo_r.end or hi_r.start > divider or divider > hi_r.end):
            divider = (hi_r.end + hi_r.start)//2
            print("for "+lo_f+" and "+hi_f)
            print("\toverlapping region fully spanned by blocks and micro blocks")
            print("\tusing mid-point as divider")

        if not (lo_r.start <= divider and divider <= lo_r.end and \
                hi_r.start <= divider and divider <= hi_r.end):
            sys.stderr.write("""an error has occurred:
                    lo_r.start={}     lo_r.end={}
                    hi_r.start={}     hi_r.end={}
                    divider={}
            """.format( lo_r.start, lo_r.end, hi_r.start, hi_r.end, divider) )
            sys.exit(1)

    if verbose: print "\tdivider is {}".format(divider)


    with open(out,'w') as fout:
        new_block_no = 1
        block_no_map = {}
        with open(lo_f,'r') as fin:
            for line in fin:
                # skip blank lines
                if line=='\n': continue
                # pass comment lines, except ##DiscovarRegion (strip_tag)
                elif line[:len(strip_tag)] == strip_tag:  continue
                elif line[0]=='#':
                    fout.write(line)
                else:
                    tmp = line.split('\t')
                    assert( len(tmp) >=8 )
                    pos = int( tmp[1] )
                    this_block_no = block_line(line)

                    if pos <= divider:
                        if not block_no_map.has_key( this_block_no ):
                            block_no_map[this_block_no] = new_block_no
                            new_block_no += 1

                        fout.write(block_modify(line, block_no_map[this_block_no]) )
                    else:
                        nDeleted+=1;

        block_no_map = {}       # reset map, but continue the block numbering
        with open(hi_f,'r') as fin:
            for line in fin:
                if line=='\n' or line[0]=='#': continue
                else:
                    tmp = line.split('\t')
                    assert( len(tmp) >=8 )
                    pos = int( tmp[1] )
                    this_block_no = block_line(line)

                    if pos >= divider:
                        if not block_no_map.has_key(this_block_no):
                            block_no_map[this_block_no] = new_block_no
                            new_block_no += 1

                        fout.write(block_modify( line, block_no_map[this_block_no]))
                    else:
                        nDeleted+=1
    return nDeleted


################################################################################
# find_vcf_files - walk a directory tree looking for VCF files
################################################################################
def find_vcf_files( path, pattern ):

    print "scanning for VCF files matching pattern {}".format(pattern)

    vcf_files=[]

    for root, dirs, files in os.walk( path ):
        vcf_files.extend([ os.path.join(root,f) for f in files if \
                fnmatch.fnmatch( f, pattern ) ])

    return vcf_files


################################################################################
#
#    M A I N
#
################################################################################
def main( argv = [__name__] ):
   default_pattern = '*.variant.filtered.vcf'
   parser=argparse.ArgumentParser(description="Take vcf files generated by Discovar for disjoint or slightly-overlapping regions")
   parser.add_argument('--output', '-o', help='output file', nargs=1, required=True)
   parser.add_argument('--verbose', '-v', help='print verbose output', action='store_true')
   parser.add_argument('--pattern', '-p',
       help='pattern of VCF files to accept (default: {})'.format(default_pattern),
       default=default_pattern)
   parser.add_argument('inputdir', help='top-level directory in which to find VCF files produced by Discovar', nargs=1)
   args=parser.parse_args(argv[1:])

   print "scanning {} for VCF files...".format( args.inputdir[0] )
   sys.stdout.flush()
   inputFiles = find_vcf_files( args.inputdir[0], args.pattern )

   sys.stderr.flush()
   sys.stdout.flush()
   return merge_vcfs( inputFiles, args.output[0], args.verbose )


if __name__ == '__main__':
    sys.exit(main(sys.argv))
