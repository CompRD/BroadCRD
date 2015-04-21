
import sys
from Bio import SeqIO



if __name__ == '__main__':

  fName = sys.argv[1];

  data=dict()

  f=open(fName,'r');

  var_idx=0
  alt_idx=1

  alt_start_idx=2

  alt_end_idx=3

  alt_len_idx=4

  alt_rc_idx=5

  contig_idx=6

  nBlock_idx=10

  first_gap_idx=11
  first_mis_idx=13

  sBuffer=f.readline()

  print "#var\t#alt\t#alt_start\t#alt_end\t#alt_length\t#alt_rc\t#contig\t#contig_start\t#contig_end\t#contig_length\t#num_block"
  while len(sBuffer)>0:
    sStripped = sBuffer.strip('\n').split()

    var = sStripped[var_idx]

    alt = int(sStripped[alt_idx])
    if alt == 0:
      assert( var not in data.keys() )

      assert( int(sStripped[alt_start_idx]) == 0)
      assert( int(sStripped[alt_end_idx]) == int(sStripped[alt_len_idx]))

      assert( int(sStripped[nBlock_idx]) == 1)
      assert( int(sStripped[first_gap_idx]) == 0)
      assert( int(sStripped[first_mis_idx]) == 0)

      data[var] = [ sStripped ]
      print sBuffer.strip('\n')
    else:
      assert( var in data.keys() )

      if  sStripped[contig_idx] == data[var][0][contig_idx] :
        data[var].append(  sStripped )
        print sBuffer.strip('\n')

    sBuffer=f.readline()

  f.close()

