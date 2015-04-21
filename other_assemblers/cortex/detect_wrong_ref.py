
import sys
from Bio import SeqIO
import re

if __name__ == '__main__':

  nArg =  len(sys.argv) -1;


  if nArg <1 or nArg>5:
      sys.stderr.write( "cmd vcf reference\n")
      exit(1)

  vName = sys.argv[1];
  rName = sys.argv[2];

  with open(rName,"r") as handle: ref_dict = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))

  chr_idx=0
  pos_idx=1
  id_idx=2
  ref_idx=3

  nErr=0
  nTot=0

  nErrBC61=0
  nErrBC31=0
  nErrPD31=0

  nTotBC61=0
  nTotBC31=0
  nTotPD31=0

  with open(vName,"r") as handle:
      for line in handle:
          if len(line)==0 or line[0]=='#': continue
          nTot+=1
          buffer = line.split()


          bBC61 = re.match('UNION_BC_k61',buffer[id_idx]) is not None
          bBC31 = re.match('UNION_BC_k31',buffer[id_idx]) is not None
          bPD31 = re.match('UNION_PD_k31',buffer[id_idx]) is not None
          if bBC61: nTotBC61+=1
          if bBC31: nTotBC31+=1
          if bPD31: nTotPD31+=1


          ref = buffer[ref_idx]
          pos = int(buffer[pos_idx])
          chr = buffer[chr_idx]

          begin = pos-1
          end   = begin+len(ref)

          ref_in_dict = str(ref_dict[chr].seq[ begin:end])

          if ref!=ref_in_dict:
              print "hg19_ref=%s %s %s"%(str(ref_dict[chr].seq[ begin-10:begin])
                                        ,ref_in_dict
                                        ,str(ref_dict[chr].seq[ end:end+10]))
              print line
              nErr+=1
              assert int(bBC61)+int(bBC31)+(bPD31) == 1
              if bBC61: nErrBC61+=1
              if bBC31: nErrBC31+=1
              if bPD31: nErrPD31+=1

  if nTotPD31>0: sys.stdout.write("PD31:  %d/%d=%e\n"%(nErrPD31,nTotPD31,nErrPD31*1.0/nTotPD31))
  if nTotBC31>0: sys.stdout.write("BC31:  %d/%d=%e\n"%(nErrBC31,nTotBC31,nErrBC31*1.0/nTotBC31))
  if nTotBC61>0: sys.stdout.write("BC61:  %d/%d=%e\n"%(nErrBC61,nTotBC61,nErrBC61*1.0/nTotBC61))
  if nTot>0:     sys.stdout.write("Total: %d/%d=%e\n"%(nErr,nTot,nErr*1.0/nTot))
