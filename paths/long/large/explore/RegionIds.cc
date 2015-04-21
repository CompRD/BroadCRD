///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Find the NA12878 newchem reads that go in a region, as specified by hg19 
// coordinates.

#include "FastIfstream.h"
#include "MainTools.h"
#include "paths/long/large/ReadNameLookup.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(X);
     EndCommandArguments;

     String iroot = "/wga/scr4/vendor/illumina/2014-07-04/Conversion";
     String abam = ".aligned.sorted.bam";
     vec<String> bam(2);
     bam[0] = iroot 
          + "/140613_HSQ1185_0757_BH9TEUADXX/1/L1_NA12878/H9TEUADXX.1" + abam;
     bam[1] = iroot
          + "/140613_HSQ1185_0757_BH9TEUADXX/2/L2_NA12878/H9TEUADXX.2" + abam;

     vec<String> rids;
     for ( int b = 0; b < 2; b++ )
     {    fast_pipe_ifstream in( "samtools view " + bam[b] + " " + X );    
          String line, readname;
          int flags;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               istringstream iline( line.c_str( ) );
               iline >> readname >> flags;
               Bool first = ( flags & 64 );
               if (first) readname += ".1";
               else readname += ".2";
               rids.push_back(readname);    }    }
     UniqueSort(rids);
     for ( int i = 0; i < rids.isize( ); i++ )
          rids[i][7] = '_', rids[i][17] = '_';

     readname_lookup look;
     String global = "/wga/scr4/wg_projects/H.sapien/NA12878_newchem";
     BinaryReader::readFile( global + "/frag_reads_orig.names.fixed.idx", &look );

     vec<int64_t> ids;
     for ( int i = 0; i < rids.isize( ); i++ )
          ids.push_back( look.GetReadId( rids[i] ) );
     Sort(ids);

     cout << printSeq(ids) << endl;    }
