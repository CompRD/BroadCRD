
import sys
import subprocess
from Bio import SeqIO



if __name__ == '__main__':

  nArg =  len(sys.argv) -1;

  if( nArg != 2 and nArg!=3):
      print nArg
      sys.stderr.write( "usage: " + sys.argv[0] + " contigs_fasta ref_head [fasta_tag]\n")
      sys.exit(1)


  contigsFN = sys.argv[1];
  ref_head = sys.argv[2];

  qlut_tgt = "SEQS="+contigsFN+" L="+ref_head+".lookup"

  qlut_cmd = "QueryLookupTable K=12 MM=12 MC=0.15 PARSEABLE=True SMITH_WAT=True QUERY_NAMING=from_record " + qlut_tgt + " | grep '^QUERY' "

  qlut_proc=subprocess.Popen(qlut_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
  qlut_proc.wait()


  idx_contig=1

  idx_src_begin=2
  idx_src_end=3
  idx_src_n=4

  idx_src_rc=5
  idx_des_begin=7
  idx_des_end=8
  idx_des_n=9
  idx_n_blk=10
  idx_blk_start=11

  alignments =[ ]

  appearances =dict()

  ref_len = -1


  for entry in qlut_proc.stdout:
      buffer=entry.split()

      nBlocks = int(buffer[idx_n_blk])
      contig_name=buffer[idx_contig]

      if contig_name in appearances:
          appearances[contig_name]+=1
      else:
          appearances[contig_name] =1

      if ref_len < 0:
          ref_len = int(buffer[idx_des_n])
      else:
          assert( ref_len == int(buffer[idx_des_n]) )
          
          

      n_snp = 0
      n_ids = 0
      for bb in range(nBlocks):
          n_snp += abs(int( buffer[idx_blk_start + 3*bb +2] ))
          n_ids += abs(int( buffer[idx_blk_start + 3*bb ] ))

      alignments.append( (int(buffer[idx_des_begin]) #0
                         ,int(buffer[idx_des_end])   #1
                         ,contig_name                #2
                         ,int(buffer[idx_src_rc])!=0 #3
                         ,int(buffer[idx_src_begin]) #4
                         ,int(buffer[idx_src_end])   #5
                         ,int(buffer[idx_src_n])     #6
                         ,n_snp                      #7
                         ,n_ids                      #8
                         )
                        )
  sorted_alignments=sorted(alignments)

  handle = open(contigsFN,"r")

  seq_dict = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))

  handle.close()

  end = 0 
  scaffold=""

  nSNP=0
  nIDS=0
  nLen=0
  totalN=0

  for entry in  sorted_alignments:
      loc_begin = entry[4]
      loc_end   = entry[5]
      if appearances[entry[2]]>1 or entry[1] <= entry[0] or loc_end <= loc_begin or entry[1] <= end:
          continue

      assert len(str(seq_dict[entry[2]].seq )) == entry[6] 

      if entry[0] > end:
          nNs = entry[0] - end
          scaffold+= 'N'*nNs
          totalN+=nNs
      else:
          loc_begin += end - entry[0]
      if end==0: loc_begin=0;
      assert loc_begin < loc_end
      if entry[3] :
          scaffold+=str(seq_dict[entry[2]].seq.reverse_complement())[loc_begin:loc_end]
      else:
          scaffold+=str(seq_dict[entry[2]].seq)[loc_begin:loc_end]
      nSNP+=entry[7]
      nIDS+=entry[8]
      nLen+=loc_end-loc_begin
      end = max(end,entry[1])
  nerr=nSNP+nIDS
  if ref_len > 0:
      nNs = ref_len - end
      scaffold+= 'N'*nNs
      totalN+=nNs
      penalty = totalN*10+nSNP+nIDS
  else:
      penalty=1000000000
  hint="_"+str(nLen-nerr)+'_'+str(ref_len)+'_'+str(nLen)+'_'+str(nIDS)+'_'+str(nSNP)+'_'+str(totalN)+'_'+str(penalty)
  if nArg==3:
      print ">"+sys.argv[3]+hint
  else:
      print '>scaffold_with_N_'+hint
  print scaffold
