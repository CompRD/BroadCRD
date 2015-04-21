
import sys
import subprocess
from Bio import SeqIO


def SummarizeAlignments(seqs,lookup_head,qlog=None,qlog_header=""):
  qlut_tgt = "SEQS="+seqs+" L="+lookup_head+".lookup"
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
  ref_len=-1
  for entry in qlut_proc.stdout:
      buffer=entry.split()

      if qlog is not None:
          qlog.write("%s "%qlog_header)
          for bb in range(idx_blk_start):
              qlog.write("%s "%buffer[bb])
          qlog.write("\n")


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

#  for entry in  sorted_alignments: print entry
  last_part=0
  last_end=0

  nSNP=0
  nIDS=0
  nContigBases=0
  nAlignedBases=0


  for entry in  sorted_alignments:
      idx = entry[2].rfind('part')
      assert idx>0 and idx<len(entry[2])
      part = int( entry[2][idx+4:] )
      if( part != last_part+1): continue;
      last_part=part
#      assert entry[0] >= last_end
      if entry[0] < last_end:
         nIDS+= last_end-entry[0]
         sys.stderr.write("\nwarning: overlaping contig\n")
      assert entry[1] >  entry[0]
      assert entry[4] == 0 or entry[0]==0
      assert entry[5] == entry[6] or entry[1]==ref_len
      last_end = entry[1]

      nSNP+=entry[7]
      nIDS+=entry[8]
      nContigBases+=entry[6]
      nAlignedBases+=entry[5]-entry[4]
  return (last_part-1,last_part,nAlignedBases,nContigBases,nSNP,nIDS)
      

if __name__ == '__main__':

  nArg =  len(sys.argv) -1;

  if  nArg != 2 :
      sys.stderr.write( "usage: " + sys.argv[0] + " before after \n")
      sys.exit(1)

  with open(sys.argv[1],"r") as handle:
      before_dict = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
  with open(sys.argv[2],"r") as handle:
      after_dict = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))

  ref_dir="/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin"

  ignored_fosmids=(16,33,35,36,38,67)

  log=[]

  qlog = open("query.log","w")

  for after_tag in after_dict:
      assert after_tag[0]=='F'
      idx = after_tag.find('_')
      assert idx >1 and idx < len(after_tag)
      fosmid_id = int(after_tag[1:idx])

      if fosmid_id in ignored_fosmids: continue

      sys.stderr.write("working on fosmid %d: "%fosmid_id)

      idx = after_tag.rfind('|')
      assert idx >=0 and idx < len(after_tag)
      before_tag=after_tag[:idx]
      assert before_tag in before_dict
      SeqIO.write( before_dict[before_tag], "before.fasta", "fasta")
      SeqIO.write( after_dict[after_tag], "after.fasta", "fasta")

      bN = str(before_dict[before_tag].seq).strip('N').count('N')
      aN = str(after_dict[after_tag].seq).strip('N').count('N')

      subprocess.check_call("BreakAtNs MIN_TO_BREAK=1 IN=%s.fasta OUT=%s.noN.fasta"%('before','before')
                           ,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
      subprocess.check_call("BreakAtNs MIN_TO_BREAK=1 IN=%s.fasta OUT=%s.noN.fasta"%('after','after')
                           ,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
      subprocess.check_call("MakeLookupTable SOURCE=%s/fos.%d.fasta OUT_HEAD=./fos.%d LO=True"%(ref_dir,fosmid_id,fosmid_id)
                           ,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
      loc_log=[fosmid_id,]
      for entry in ("before.noN.fasta","after.noN.fasta"):
          sys.stderr.write("%s "%entry)
          header="%3d %s\t"%(fosmid_id,entry[:entry.find('.')])
          loc_log.append(SummarizeAlignments(entry,"fos.%d"%fosmid_id,qlog,header))
      loc_log.append(bN)
      loc_log.append(aN)
      sys.stderr.write("\n")
      log.append(loc_log)
  log=sorted(log)
  print "# fosmid stage nGap nContig nAlignedBases nContigBases nSNP nIDS nN"
  for entry in log:
      sys.stdout.write( "%3d "%entry[0]+" before " )
      for value in entry[1]: sys.stdout.write( str(value)+"\t" )
      sys.stdout.write("%d\n"%entry[3])
      sys.stdout.write( "%3d "%entry[0]+" after  " )
      for value in entry[2]: sys.stdout.write( str(value)+"\t" )
      sys.stdout.write("%d\n"%entry[4])



#  return (last_part-1,last_part,nAlignedBases,nContigBases,nSNP,nIDS)
  qlog.close()
  exit(0)
