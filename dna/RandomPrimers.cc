/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// RandomPrimers.  Given a sequence (the "original primer"), find other random
// primers that have almost the same base composition as the original primer and 
// no homopolymers (but see next).
//
// If INTERNAL_HOMOPOLYMERS = False, do not allow homopolymers, except at end.
// 
// If MAX_HOMOPOLYMER is positive, only return primers have homopolymers of length
// at most MAX_HOMOPOLYMER.
//
// If MIN_T or MAX_T are specified, require that the melting temperature of
// the primer satisfies the given bound, where the melting temperature is computed
// using Tm_NearestNeighbor (probably not completely right).

#include "dna/DNAHybridization.h"
#include "MainTools.h"
#include "random/Random.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(ORIG_PRIMER);
     CommandArgument_Int(NPRIMERS);
     CommandArgument_Bool(INTERNAL_HOMOPOLYMERS);
     CommandArgument_Int(MAX_HOMOPOLYMER);
     CommandArgument_Double_OrDefault(MIN_T, -1.0);
     CommandArgument_Double_OrDefault(MAX_T, -1.0);
     EndCommandArguments;

     int A = 0, C = 0, G = 0, T = 0;
     int n = ORIG_PRIMER.size( );
     for ( int i = 0; i < n; i++ )
     {    char c = ORIG_PRIMER[i];
          if ( c == 'A' ) ++A;
          if ( c == 'C' ) ++C;
          if ( c == 'G' ) ++G;
          if ( c == 'T' ) ++T;    }
     for ( int i = 0; i < NPRIMERS; i++ )
     {    String primer;
          int a = A, c = C, g = G, t = T;
          for ( int j = 0; j < n; j++ )
          {    char base;
               while(1)
               {    int r = randomx( ) % ( n - j );
                    if ( r < a ) base = 'A';
                    else if ( r < a + c ) base = 'C';
                    else if ( r < a + c + g ) base = 'G';
                    else base = 'T';
                    if (INTERNAL_HOMOPOLYMERS) break;
                    if ( j == 0 ) break;
                    if ( base != primer[j-1] ) break;
                    if ( c + g + t == 0 || a + g + t == 0 || a + c + t == 0
                         || a + c + g == 0 )
                    {    break;    }    }
               if ( base == 'A' ) --a;
               else if ( base == 'C' ) --c;
               else if ( base == 'G' ) --g;
               else --t;
               primer += base;    }
          Bool bad = False;
          for ( int j = 0; j < n; j++ )
          {    int k;
               for ( k = j + 1; k < n; k++ )
                    if ( primer[k] != primer[j] ) break;
               if ( k - j > MAX_HOMOPOLYMER ) bad = True;    }
          if ( MIN_T >= 0.0 || MAX_T >= 0.0 )
          {    double T = Tm_NearestNeighbor(primer);
               if ( MIN_T >= 0.0 && T < MIN_T ) bad = True;
               if ( MAX_T >= 0.0 && T > MAX_T ) bad = True;    }
          if (bad)
          {    --i;
               continue;    }
          cout << "[" << i+1 << "] " << primer << "\n";    }    }
