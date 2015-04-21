/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FastbToFastqs.  Split a fastb file up into a bunch of fastq files.  Assign
// constant quality scores.

#include "Basevector.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "math/Functions.h"
#include "random/Random.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FASTB);
     CommandArgument_String_Doc(FASTQ,
          "if FASTQ=x.fastq, output files are x.n.fastq, n = 1, 2, ...");
     CommandArgument_String(READ_NAME_HEAD);
     CommandArgument_Int_Doc(N, "number of fastq files to generate");
     CommandArgument_Int_OrDefault_Doc(Q, 30, "quality score to assign");
     CommandArgument_Bool_OrDefault_Doc(ADD_FAKE, False,
          "add a fake read to work around a bug in MAQ");
     EndCommandArguments;

     vecbasevector bases(FASTB);
     vecqualvector quals;
     Mimic( bases, quals );
     for ( vecqvec::size_type i = 0; i < quals.size( ); i++ )
     {    for ( qvec::size_type j = 0; j < quals[i].size( ); j++ )
               quals[i][j] = Q;    }

     basevector fake(80);
     for ( int i = 0; i < 80; i++ )
          fake.Set( i, randomx( ) % 4 );

     int per = ( bases.size( ) / N );
     if ( bases.size( ) % N > 0 ) ++per;

     int fid = 0;
     for ( size_t x = 0; x < bases.size( ); x += per )
     {    Ofstream( out, FASTQ.Before( "fastq" ) + ToString(fid++) + ".fastq" );
          if (ADD_FAKE)
          {    String fakeq( fake.size( ) );
               for ( int j = 0; j < fakeq.isize( ); j++ )
                    fakeq[j] = 30 + 33;
               out << "@" << "fake:" << fid << "\n"
                    << fake.ToString( ) << "\n+\n" << fakeq << "\n";    }
          for ( size_t i = x; i < Min( x + per, bases.size( ) ); i++ )
          {    String b = bases[i].ToString( );
               String q( b.size( ) );
               for ( qvec::size_type j = 0; j < quals[i].size( ); j++ )
                    q[j] = quals[i][j] + 33;
               out << "@" << READ_NAME_HEAD << ":" << i << "\n"
                    << b << "\n+\n" << q << "\n";    }    }    }

