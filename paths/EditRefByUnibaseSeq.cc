///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Edit a reference sequence using a sequence of unibases from some assembly.
// The first and last unibases in the sequence must be uniquely located on the
// assembly.

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "paths/UnibaseUtils.h"

int main(int argc, char *argv[])
{    
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(UNIBASES, "all_reads");
     CommandArgument_String_Doc(GENOME, "fastb file for genome");
     CommandArgument_Int(K);
     CommandArgument_String_Doc(SEQ,
          "ordered list of integers in ParseIntSet format");
     CommandArgument_Bool_OrDefault_Doc(RC, False,
          "if the sequence is reverse complement to the genome");
     EndCommandArguments;

     // Define directories.
     
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Load genome.

     vecbasevector genome(GENOME);

     // Load unibases and generate ancillary data.

     vecbasevector unibases( run_dir + "/" + UNIBASES 
          + ".unibases.k" + ToString(K) );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc );

     // Get list of unibases.

     vec<int> u;
     ParseIntSet( SEQ, u, False );
     if (RC)
     {    u.ReverseMe( );
          for ( int j = 0; j < u.isize( ); j++ )
               u[j] = to_rc[ u[j] ];    }

     // Locate first and last unibases on reference.

     vec< pair<int,int> > locs1, locs2;
     String s1 = unibases[ u.front( ) ].ToString( );
     String s2 = unibases[ u.back( ) ].ToString( );
     for ( size_t g = 0; g < genome.size( ); g++ )
     {    String G = genome[g].ToString( );
          for ( int p = 0; p <= G.isize( ) - s1.isize( ); p++ )
               if ( G.Contains( s1, p ) ) locs1.push( g, p );
          for ( int p = 0; p <= G.isize( ) - s2.isize( ); p++ )
               if ( G.Contains( s2, p ) ) locs2.push( g, p );    }
     if ( !locs1.solo( ) ) FatalErr("Left end unibase not uniquely located.");
     if ( !locs2.solo( ) ) FatalErr("Right end unibase not uniquely located.");
     int g = locs1[0].first;
     if ( locs2[0].first != g )
     {    cout << "End unibases not located on same chromosome." << endl;
          exit(1);    }
     int start = locs1[0].second; 
     int stop = locs2[0].second + unibases[ u.back( ) ].isize( );
     if ( !( start <= stop ) )
     {    cout << "Start must be <= stop." << endl;
          exit(1);    }
     basevector g1( genome[g], 0, start );
     basevector g2 = unibases[ u[0] ];
     for ( int j = 1; j < u.isize( ); j++ )
     {    g2.resize( g2.isize( ) - (K-1) );
          g2 = Cat( g2, unibases[ u[j] ] );    }
     basevector g3( genome[g], stop, genome[g].isize( ) - stop );
     genome[g] = Cat( g1, g2, g3 );
     for ( size_t i = 0; i < genome.size( ); i++ )
          genome[i].Print( cout, i );    }
