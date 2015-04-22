#!/usr/bin/env python

import sys
import os
import subprocess
c_idx = 0
s_idx = 1
e_idx = 2
vcf_idx = 3
vs_idx = 4
ve_idx = 5

def append_blocks( blocks, micro_blocks, v_file):
    with open(v_file,'r') as f:
        start=-1
        end=-1

        for line in f:
            entries = line.strip('\n').split()
            if len(entries) < 1: continue;

            if entries[0]=="BLOCK":
                if start >= 0:
                    blocks.append( (start,end) )
                start=-1
                end=-1
            elif  entries[0]=="VAR":
                coor = int(entries[2].split(':')[1])
                ref = entries[3].strip()
                if start <0: start = coor
                end = max(coor + len(ref) - 1,end)
                micro_blocks.append( (coor,end) )
        if start >= 0:
            blocks.append( (start,end) )




def divider_from_blocks( blocks, divider_in):
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




def merge_vcf_avoid_block(lo,hi,out):
    nDeleted=0;
    assert( lo[e_idx] > hi[s_idx] )
    divider = (lo[e_idx] + hi[s_idx])//2

    blocks=[];
    micro_blocks=[];
    append_blocks( blocks, micro_blocks, lo[ve_idx] )
    append_blocks( blocks, micro_blocks, hi[vs_idx] )

    if( len(blocks) > 0):
        divider = divider_from_blocks(blocks,divider)
        if( lo[s_idx] > divider or divider > lo[e_idx] or hi[s_idx] > divider or divider > hi[e_idx]):
            sys.stderr.write("for "+lo[ve_idx] + " and " + hi[vs_idx] + ",\n")
            sys.stderr.write("WARNING: overlapping region fully spanned by blocks\n")
            sys.stderr.write("WARNING: seeking divider in micro blocks\n")
            divider = (lo[e_idx] + hi[s_idx])//2
            divider = divider_from_blocks(micro_blocks,divider)

        if( lo[s_idx] > divider or divider > lo[e_idx] or hi[s_idx] > divider or divider > hi[e_idx]):
            divider = (lo[e_idx] + hi[s_idx])//2
            sys.stderr.write("for "+lo[ve_idx] + " and " + hi[vs_idx] + ",\n")
            sys.stderr.write("WARNING: overlapping region fully spanned by blocks and micro blocks\n")
            sys.stderr.write("WARNING: using mid-point as divider\n")
        assert( lo[s_idx] < divider and divider < lo[e_idx] and hi[s_idx] < divider and divider < hi[e_idx])
    with open(out,'w') as fout:
        with open(lo[vcf_idx],'r') as fin:
            for line in fin:
                if line=='\n': continue
                elif line[0]=='#': fout.write(line)
                else:
                    tmp = line.split('\t')
                    assert( len(tmp) >=8 )
                    pos = int( tmp[1] )
                    if pos <= divider: fout.write(line)
                    else: nDeleted+=1;
        with open(hi[vcf_idx],'r') as fin:
            for line in fin:
                if line=='\n' or line[0]=='#': continue
                else:
                    tmp = line.split('\t')
                    assert( len(tmp) >=8 )
                    pos = int( tmp[1] )
                    if pos >= divider: fout.write(line)
                    else: nDeleted+=1
    return nDeleted

def merge_vcf(lo,hi,out):

    nDeleted=0

    assert( lo[s_idx] < hi[s_idx])
    if lo[e_idx] < hi[s_idx]:
        os.system("(cat "+ lo[vcf_idx] + " | grep ^# ; cat " 
                         + lo[vcf_idx] + " | grep -v ^# ; cat " 
                         + hi[vcf_idx] + " | grep -v ^# ) > " + out );
    else:
#        os.system("~/wga/src/BroadCRD/paths/long/distributed/merge_vcf.sh "+ lo[vcf_idx] + " " + hi[vcf_idx] + " " + out );
        nDeleted+=merge_vcf_avoid_block(lo,hi,out)
    return nDeleted;



if __name__ == '__main__':
    sys.stderr.write("This code is obsolete. Please use another merger.\n");
    exit (1)




    nArg = len(sys.argv) - 1

    if nArg !=2:
        sys.stderr.write("\n");
        sys.stderr.write(sys.argv[0] + " data_directory target_directory.\n");
        sys.stderr.write("\n");
        sys.stderr.write("Merge the VCF files in the data directory to one vcf file per chromosome in target_directory\n");
        sys.stderr.write("\n");
        sys.stderr.write("Assumptions:\n");
        sys.stderr.write("a) .variant and .variant.vcf files are stored in sub-directories data_directory/chromosome/start-end\n");
        sys.stderr.write("b) for the list regions defined by start-end, each region overlap with at most one other region to the left,\n")
        sys.stderr.write("   and at most one other region to the right\n")
        sys.stderr.write("c) if a file prefix.variant.vcf exists, a .variant file of the name prefix.variant exists\n");
        sys.stderr.write("\n");
        exit(1)

    d_dir = sys.argv[1]
    t_dir = sys.argv[2]

    if not os.path.isdir(d_dir ):
        sys.stderr.write( d_dir+ " does not exist as a directory.\n")
        exit(1)

    if os.path.exists(t_dir ):
        sys.stderr.write( t_dir+ " exists and cannot serve as the target directory.\n")
        exit(1)
    os.makedirs(t_dir)

    old_dir=t_dir+"/t0"
    new_dir=t_dir+"/t1"
    os.makedirs(old_dir)
    os.makedirs(new_dir)

#    out_file = t_dir+"/out.vcf"

    file_list_file=t_dir+"/in_files"

    sys.stderr.write("searching for vcf files in " + d_dir);
    os.system( "find " + d_dir + " -name '*vcf' > " + file_list_file );
    sys.stderr.write(" done\n")

    with open(file_list_file,'r') as f:
        file_list = f.readlines()


    data_list = []

    for entry in file_list:
        vcf_file = entry.strip('\n')

        slash_1 = vcf_file.find('/',len(d_dir))
        slash_2 = vcf_file.find('/',slash_1+1)
        slash_3 = vcf_file.find('/',slash_2+1)

        chromosome = vcf_file[slash_1+1:slash_2]
        bound = vcf_file[slash_2+1:slash_3].split('-')

        var_file = vcf_file[:-4]

        data_list.append( (chromosome, int(bound[0]), int(bound[1]), vcf_file, var_file, var_file ) )



    data_list.sort()
    for ii in range( 1, len(data_list) ):
        assert( data_list[ii][1] > data_list[ii-1][1] or data_list[ii-1][0] != data_list[ii][0])



    old_list = data_list;
    nDeleted=0

    nStepDone=0;

    while( len(old_list) > 1 ):
#        if(nStepDone == 4): break;
        nOld = len(old_list)
        sys.stderr.write( "merging " + str(nOld) + " entries.\n")

        new_list=[]
        os.system( "rm -f " + new_dir + "/*" );

        if nOld > 1 :
            bDone = True;
            for ii in range(0, nOld-1):
                if old_list[ii][0] == old_list[ii+1][0] :
                    bDone=False;
                    break;
            if bDone:
#                for ii in range(0, nOld):
#                    new_file = new_dir + "/c" + str(ii) + ".vcf"
#                    os.system( "cp " + old_list[ii][vcf_idx] + " " + new_file);
                break;
                

#        print "bayo old list:"
#        for entry in old_list: print entry
        ii=0
        while ( ii < nOld ):
            if ( ii + 1 == nOld):
                new_file = new_dir + "/" + str(ii) + ".vcf"
                new_list.append( ( old_list[ii][c_idx]
                                 , old_list[ii][s_idx]
                                 , old_list[ii][e_idx]
                                 , new_file
                                 , old_list[ii][vs_idx]
                                 , old_list[ii][ve_idx]
                                 )
                               )
                os.system( "cp " + old_list[ii][vcf_idx] + " " + new_file);
                ii+=1
            else:
                next = ii + 1
                if old_list[ii][c_idx] == old_list[next][c_idx]:
                    new_file = new_dir + "/" + str(ii) + "_" + str(next) + ".vcf"
                    new_list.append( ( old_list[ii][c_idx]
                                     , old_list[ii][s_idx]
                                     , old_list[next][e_idx]
                                     , new_file
                                     , old_list[ii][vs_idx]
                                     , old_list[next][ve_idx]
                                     )
                                   )
                    nDeleted+=merge_vcf( old_list[ii] , old_list[next] , new_file )
                    ii=next+1
                else:
                    new_file = new_dir + "/" + str(ii) + ".vcf"
                    new_list.append( ( old_list[ii][c_idx]
                                     , old_list[ii][s_idx]
                                     , old_list[ii][e_idx]
                                     , new_file
                                     , old_list[ii][vs_idx]
                                     , old_list[ii][ve_idx]
                                     )
                                   )
                    os.system( "cp " + old_list[ii][vcf_idx] + " " + new_file);
                    ii+=1

#        print "bayo old list:"
#        for entry in old_list: print entry
#        print "bayo new list:"
#        for entry in new_list: print entry
        old_list=new_list
        old_dir,new_dir=new_dir,old_dir
        nStepDone+=1;
    os.system( "cp " + old_dir + "/*vcf " + t_dir)
    sys.stderr.write( "number of discard entries: "+ str(nDeleted) + "\n");
    sys.stderr.write( "output file directory: "+ t_dir+"\n");

#    bCheck=False
#    if bCheck:
#        tmp_vcf = t_dir + "/tmp.vcf"
#        tmp_gz = t_dir + "/tmp.vcf.gz"
#        tmp_tbi = t_dir + "/tmp.vcf.gz.tbi"
#        sorted_gz = t_dir + "/sorted.vcf.gz"
#        os.system( "cat "+out_file+ " | vcf-sort | bgzip -c > " + sorted_gz)
#        os.system( "tabix " + sorted_gz )
#
#        for entry in data_list:
#            vcf_file=entry[vcf_idx]
#            print entry[s_idx],entry[e_idx],vcf_file
#            os.system("rm -f "+tmp_vcf)
#            os.system("rm -f "+tmp_gz)
#            os.system("rm -f "+tmp_tbi)
#            os.system("cp "+vcf_file+" "+tmp_vcf)
#            os.system("~/packages/tabix-0.2.6/bgzip -c "+tmp_vcf+" > "+tmp_gz)
#            os.system("tabix "+tmp_gz)
#            os.system("vcf-isec -f -c "+tmp_gz + " " + sorted_gz + " | grep -v ^#")
