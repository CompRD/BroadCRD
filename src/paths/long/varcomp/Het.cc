///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Het: estimate rate of heterozygous variants across the genome.  Uses GATK
// calls from 100 base reads from NA12878.

#include "Bitvector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "math/Functions.h"
#include "random/Random.h"

int main( )
{    RunTime( );
     vecbitvector amb( "/wga/scr4/bigrefs/human19/genome.lookup.fastamb" );
     fast_ifstream in( "/wga/scr4/human_data/CEPH/NA12878_S1.genome.vcf" );
     String line;
     vec<String> fields;
     const int bin = 5000000;
     vec<vec<int>> hets(23);
     for ( int c = 0; c < 23; c++ )
          hets[c].resize( ( (int) amb[c].size( ) + bin - 1 ) / bin );
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( !line.Contains( "PASS" ) ) continue;
          if ( !line.Contains( "/" ) ) continue;
          if ( line.Contains( "0/0" ) || line.Contains( "1/1" ) ) continue;
          Tokenize( line, {'\t'}, fields );
          int pos = fields[1].Int( ) - 1;
          if ( !line.Contains( "0/1" ) && !line.Contains( "1/2" ) ) PRINT(line);
          String chr = fields[0].After( "chr" );
          int ichr;
          if ( chr == "X" ) ichr = 22;
          else if ( chr == "Y" ) continue;
          else ichr = chr.Int( ) - 1;
          hets[ ichr ][ pos/bin ]++;    }
     for ( int c = 0; c < 23; c++ )
     for ( int p = 0; p < hets[c].isize( ); p++ )
     {    if ( hets[c][p] != 0 )
          {    int n = 0;
               for ( int j = p*bin; j < (p+1)*bin; j++ )
                    if ( j >= (int) amb[c].size( ) || amb[c][j] ) n++;
               // if ( n >= 20000 ) continue;
               cout << "chr" << ( c == 22 ? "X" : ToString(c+1) )
                    << "." << p*bin/1000000 << "-" << (p+1)*bin/1000000 
                    << " Mb: " << hets[c][p] / double( bin/1000000 )
                    << " (amb=" << n << ")" << endl;    }    }    }
