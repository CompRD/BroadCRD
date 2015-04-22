///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"

/**
 * ReportUnmapped
 *
 * Report a list of unmapped nhoods (nhoods for which the seed was not
 * aligned to an assembly or reference).
 *
 * SEEDS: full path name of seeds.ids
 * ALIGNS: aligns of seeds onto a reference
 * BASES: if not empty, load fastb and report lentghs of unmapped seeds
 * BRIEF: only print brief summary
 */
int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String( SEEDS );
  CommandArgument_String( ALIGNS );
  CommandArgument_String_OrDefault( BASES, "" );
  CommandArgument_Bool_OrDefault( BRIEF, False );
  EndCommandArguments;

  // Load.
  READ( SEEDS, vec<int>, ids );
  ForceAssert( is_sorted( ids.begin( ), ids.end( ) ) );

  vec<look_align> hits;
  LoadLookAligns( ALIGNS, hits );

  // Count number of hit per seed.
  vec<int> counts( ids.size( ), 0 );
  for (int ii=0; ii<hits.isize( ); ii++) {
    int qid = hits[ii].query_id;
    vec<int>::iterator it = lower_bound( ids.begin( ), ids.end( ), qid );
    ForceAssert( it != ids.end( ) );
    counts[distance( ids.begin( ), it )] += 1;
  }
  
  // Print anomalies.
  int n_unmapped = 0;
  int n_multiplet = 0;
  vec<int> sel;
  for (int ii=0; ii<counts.isize( ); ii++) {
    if ( counts[ii] == 1 ) continue;
    if ( counts[ii] == 0 ) n_unmapped++;
    else n_multiplet++;
    sel.push_back( ids[ii] );
  }

  // Load bases (so can report lengths).
  vecbvec bases;
  if ( BASES != "" ) bases.SparseRead( BASES, sel, 0 );

  // Verbose report.
  if ( ! BRIEF ) {
    vec< vec<String> > table;
    vec<String> line;

    line.push_back( "nhood_id" );
    line.push_back( "unipath_id" );
    if ( BASES != "" ) line.push_back( "length (bp)" );
    line.push_back( "n_aligns" );
    table.push_back( line );

    for (int ii=0; ii<counts.isize( ); ii++) {
      if ( counts[ii] == 1 ) continue;
      line.clear( );
      line.push_back( ToString( ii ) );
      line.push_back( ToString( ids[ii] ) );
      if ( BASES != "" ) line.push_back( ToString( bases[ ids[ii] ].size( ) ) );
      line.push_back( ToString( counts[ii] ) );
      table.push_back( line );
    }

    String just;
    if ( BASES == "" ) just = "rrr";
    else just = "rrrr";

    PrintTabular( cout, table, 2, just );
    cout << endl;
  }

  // Print summary.
  double ratio_unmapped = 100.0 * SafeQuotient( n_unmapped, ids.isize( ) );
  double ratio_multiplet = 100.0 * SafeQuotient( n_multiplet, ids.isize( ) );
  cout << "SUMMARY\n\n"
       << "total nhoods: " << ids.size( ) << "\n"
       << "unmapped:     " << n_unmapped
       << " (" << ToString( ratio_unmapped, 1 ) << "%" << " of the total)\n"
       << "multiplets:   " << n_multiplet
       << " (" << ToString( ratio_multiplet, 1 ) << "%" << " of the total)\n"
       << "\n";

  // Done.
  cout << Date( ) << ": done" << endl;

}

