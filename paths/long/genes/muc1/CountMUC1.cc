///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use,misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CountMUC1.  Given a LongProto assembly of human 1:155.15M-155.17M, which includes 
// the MUC1 tandem repeat, count the number of certain repeat units that might occur 
// in the original, initially corrected (mod0) and corrected (mod) reads.  Some of 
// these repeat units are known to be in the genome whereas others should not be 
// present.  The intent is that this code might be used as an assessment tool that 
// could aid in the improvement of the error correction algorithm.
//
// Some paths are hardcoded for the moment.

// Diploid reference sequence is believed to be:
//
// [maternal] C-X-D-E-C-F-X-X-A-B-D-E-C-X-X-W-A-A-B-X-X-G-A-B-X-X-X-X-X-X-V
//
// [paternal] C-X-D-E-C-X-X-X-A-A-B-X-X-X-X-X-A-A-A-B-X-X-X-X-X-A-A-X-X-X-X-
//            X-X-X-X-A-A-B-X-X-X-X-B-B-X-X-A-A-B-X-X-X-B-A-A-X-X-X-X-X-X-X-V.

#include "Basevector.h"
#include "MainTools.h"

int main( )
{
     vecbasevector bases( "tmp.xxx/frag_reads_orig.fastb" );
     vecbasevector bases0( "tmp.xxx/frag_reads_mod0.fastb" );
     vecbasevector bases1( "tmp.xxx/frag_reads_mod.fastb" );

     vec< pair<String,String> > units;

     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA","X");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGAgAgCAGGCCGGCCCCGGGCTCCACCGCgCCCgCA","A");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGAgAgCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA","B");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCaA","C");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCcGCCCCGGGCTCCACCGCCCCCCCA","D");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCcGCCCCGGGCTCCACCGCgCCCgCA","E");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCaCA","F");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCgCCCgCA","G");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCaCCCCCA","V");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCg","W");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCgCCCCA","H");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCgCCCCCA","I");
     units.push("GCCCACGGTGTCACCTCGGCaCCGGAgAgCAGGCCGGCCCCGGGCTCCACCGCgCCCgCA","J");
     units.push("GCCCACGGTGTCACCTCGGCCCCGGAgAgCAGGCCGGCCCtGGGCTCCACCGCCCCCCCA","K");

     int nu = units.size( );
     for ( int i = 0; i < nu; i++ )
     for ( int j = 0; j < units[i].first.isize( ); j++ )
          units[i].first[j] = toupper( units[i].first[j] );

     vec<String> status;
     for ( int i = 0; i < 10; i++ )
          status.push_back( "yes" );
     for ( int i = 0; i < 4; i++ )
          status.push_back( "no" );

     vec<int> abundance;
     abundance.push_back( 51, 17, 10, 5, 3, 3, 1, 1, 2, 1 );
     abundance.push_back( 0, 0, 0, 0 );

     vec<int> count( nu, 0 ), count0( nu, 0 ), count1( nu, 0 );

     for ( int i = 0; i < (int) bases.size( ); i++ )
     {    String b = bases[i].ToString( );
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) StringReverseComplement( b, b );
               for ( int j = 0; j <= b.isize( ) - 60; j++ )
               for ( int u = 0; u < nu; u++ )
                    if ( b.Contains( units[u].first, j ) ) count[u]++;    }    }

     for ( int i = 0; i < (int) bases0.size( ); i++ )
     {    String b = bases0[i].ToString( );
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) StringReverseComplement( b, b );
               for ( int j = 0; j <= b.isize( ) - 60; j++ )
               for ( int u = 0; u < nu; u++ )
                    if ( b.Contains( units[u].first, j ) ) count0[u]++;    }    }

     for ( int i = 0; i < (int) bases1.size( ); i++ )
     {    String b = bases1[i].ToString( );
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) StringReverseComplement( b, b );
               for ( int j = 0; j <= b.isize( ) - 60; j++ )
               for ( int u = 0; u < nu; u++ )
                    if ( b.Contains( units[u].first, j ) ) count1[u]++;    }    }

     vec< vec<String> > rows;
     vec<String> row;
     row.push_back( "unit", "orig", "mod0", "mod", "genome", "abundance" );
     rows.push_back(row);
     for ( int u = 0; u < nu; u++ )
     {    vec<String> row;
          row.push_back( units[u].second, ToString( count[u] ),
               ToString( count0[u] ), ToString( count1[u] ), status[u], 
               ToString( abundance[u] ) );
          rows.push_back(row);    }

     PrintTabular( cout, rows, 2, "lrrrlr" );    }
