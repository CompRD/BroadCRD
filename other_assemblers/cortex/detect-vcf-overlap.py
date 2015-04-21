
import sys



def line_equal( l, r ):
    c_idx=0
    p_idx=1
    r_idx=3
    a_idx=4
    return l[c_idx]==r[c_idx] and l[p_idx]==r[p_idx] and l[r_idx]==r[r_idx] and l[a_idx]==r[a_idx]

def output_block( block ):
    f_idx=6
    if len(block)==1 :
#        line = block[0]
#        for ii in range( len(line) ):
#            if ii !=0: sys.stdout.write('\t')
#            sys.stdout.write( line[ii] )
#        sys.stdout.write('\n');
        return


    unique_lines = []
#    sys.stdout.write(  "block size: " + str(len(block)) + '\n') 

    if len(block) > 100 :
        for line in block:
            print line
        exit (1)





    for line in block:
        bTake = True
        if line[f_idx] != 'PASS' : continue
        for uLine in unique_lines:
            if line_equal(line,uLine) :
                bTake = False
                break;
        if bTake :
            unique_lines.append(line)
    if len(unique_lines) > 1:
        sys.stdout.write(  "overlaping entries of block size: " + str(len(block)) 
                         + ". Number of unique lines: " + str(len(unique_lines)) + '\n')
#        print unique_lines

        for line in unique_lines:
            for ii in range( len(line) ):
                if ii !=0: sys.stdout.write('\t')
                sys.stdout.write( line[ii] )
            sys.stdout.write('\n');
#        print 'from'
#
#        for line in block:
#            for ii in range( len(line) ):
#                if ii !=0: sys.stdout.write('\t')
#                sys.stdout.write( line[ii] )
#            sys.stdout.write('\n');


#    for line in block:
#        for ii in range( len(line) ):
#            if ii !=0: sys.stdout.write('\t')
#            sys.stdout.write( line[ii] )
#        sys.stdout.write('\n');

if __name__ == '__main__':

  nArg =  len(sys.argv) -1;

  if( nArg != 1):
      print nArg
      sys.stderr.write( "usage: " + sys.argv[0] + " [in_vcf]\n")
      sys.exit(1)


  fName = sys.argv[1];

  last_chrom=""
  last_start=-1;
  back=-1;

  entry_printed = 0;

  block = []

  c_idx=0
  p_idx=1
  r_idx=3
  a_idx=4
  f_idx=6


  with open( fName ) as f:
      for line in f:
          if len(line) > 0 and line[0] == '#':
              sys.stdout.write( line )
          else:
              buffer = line.split()
              if buffer[f_idx] != 'PASS': continue
              loc_start = int(buffer[p_idx])
              loc_back  = loc_start + len(buffer[r_idx]) - 1
              loc_chrom = buffer[c_idx]
              assert loc_chrom!=last_chrom or loc_start >= last_start
              if loc_chrom != last_chrom :
                  print "chrom= ", loc_chrom
                  back=-1
              if loc_chrom != last_chrom or loc_start > back:
                  output_block(block)
                  block = [ ]

              block.append( buffer )

              back= max( back,loc_back )
              last_chrom=loc_chrom
              last_start=loc_start
  #deal with block
  output_block(block)
