///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SearchFastb.  Search a fastb file F for sequences containing a given sequence S 
// or its reverse complement.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"

void Search( const vecbasevector& B, String S, const Bool COUNT_ONLY,
     const Bool FW_ONLY, const Bool START0, const Bool QUIET, const int W )
{    for ( size_t i = 0; i < S.size( ); i++ )
          S[i] = toupper( S[i] );
     String Src;
     StringReverseComplement( S, Src );
     int count = 0;
     #pragma omp parallel for
     for ( size_t i = 0; i < B.size( ); i++ )
     {    String s;
          // if ( i % 5000000 == 0 && !COUNT_ONLY && !QUIET ) 
          // {
          //      #pragma omp critical
          //      {    PRINT2( i, B.size( ) );    }    }
          s = B[i].ToString( );
          if ( ( START0 && s.Contains(S,0) ) || ( !START0 && s.Contains(S) ) )
          {    
               #pragma omp critical
               {    if ( !COUNT_ONLY ) 
                         B[i].Print( cout, "seq_" + ToString(i) + "_fw" );
               ++count;    }    }
          if (FW_ONLY) continue;
          if ( ( START0 && s.Contains(Src,0) ) || ( !START0 && s.Contains(Src) ) )
          {    if (COUNT_ONLY)
               {
                    #pragma omp critical
                    {    ++count;    }    }
               else
               {    basevector b;
                    b.ReverseComplement( B[i] );
                    #pragma omp critical
                    {    b.PrintCol( cout, "seq_" + ToString(i) + "_rc", W );
                         ++count;    }    }    }    }
     cout << S << ": found " << count << " instances" << endl;    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(F);
     CommandArgument_String(S);
     CommandArgument_Bool_OrDefault(COUNT_ONLY, False);
     CommandArgument_Bool_OrDefault(FW_ONLY, False);
     CommandArgument_Bool_OrDefault(START0, False);
     CommandArgument_Int_OrDefault(SUBK, 0);
     CommandArgument_Bool_OrDefault(QUIET, False);
     CommandArgument_Int_OrDefault_Doc(W, 80, "number of bases to print per line");
     EndCommandArguments;

     vecbasevector B(F);
     if ( SUBK == 0 ) Search( B, S, COUNT_ONLY, FW_ONLY, START0, QUIET );
     else
     {    for ( int p = 0; p <= S.isize( ) - SUBK; p++ )
          {    Search( B, S.substr( p, SUBK ), 
                    COUNT_ONLY, FW_ONLY, START0, QUIET, W );    }    }    }
