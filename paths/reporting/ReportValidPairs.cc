/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "PairsHandler.h"
#include "STLExtensions.h"
#include "lookup/DetectValidPairs.h"
#include "lookup/LookAlign.h"

/**
 * ReportValidPairs
 *
 * Find all possible valid pairs (reads may own multiple aligns, in
 * which case more than one valid pairs are reported). Output is sent
 * to cout.
 *
 * READS: full path name (it needs READS.fastb, and READS.pairto)
 * HITS: full path name for look aligns
 * PROPER_ONLY: only accept pairs for which both ends align properly
 * EXTEND_ONLY: only accept pairs if right read (rc) extend left read (fw)
 * MAX_ERROR_RATE: do not take into account hits with too many errors
 * MAX_STRETCH: defines a valid insert
 * DUMP_HITS: if not empty, save all the pairs (sequentially)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( HITS );
  CommandArgument_Bool_OrDefault( PROPER_ONLY, True );
  CommandArgument_Bool_OrDefault( EXTEND_ONLY, True );
  CommandArgument_Double_OrDefault( MAX_ERROR_RATE, 0.2 );
  CommandArgument_Double_OrDefault( MAX_STRETCH, 3.5 );
  CommandArgument_String_OrDefault( DUMP_HITS, "" );
  EndCommandArguments;

  // Dir and file names.
  String bases_file = READS + ".fastb";
  String pairs_file = READS + ".pairto";

  // Load.
  int n_reads = MastervecFileObjectCount( bases_file );

  cout << Date( ) << ": loading pairs" << endl;
  phandler pairs( n_reads, pairs_file );

  cout << Date( ) << ": loading hits" << endl;
  vec<look_align_plus> hits;
  LoadLookAlignPlus( HITS, hits );
  if ( ! is_sorted( hits.begin( ), hits.end( ) ) ) {
    cout << Date( ) << ": sorting hits" << endl;
    sort( hits.begin( ), hits.end( ) );
  }
  
  // Do not run max_mult test (argument of DetectValidPairs).
  double max_mult = 0.0;
  const vec<read_pairing> &r_pairs = pairs.Pairs( );
  const vec<int> to_pair = pairs.PairId( );

  // Detect valid pairs
  cout << Date( ) << ": detecting valid pairs" << endl;
  vec< pair<int,int> > valid;
  DetectValidPairs( hits, r_pairs, to_pair, valid, MAX_STRETCH, max_mult );
  
  // Dump hits as valid pairs.
  if ( DUMP_HITS != "" ) {
    cout << Date( ) << ": saving valid pairs" << endl;

    ofstream out( DUMP_HITS.c_str( ) );
    for (int valid_id=0; valid_id<valid.isize( ); valid_id++) {
      const look_align_plus &hit1 = hits[ valid[valid_id].first ];
      const look_align_plus &hit2 = hits[ valid[valid_id].second ];
      const look_align_plus &left = hit1.rc1 ? hit2 : hit1;
      const look_align_plus &right = hit1.rc1 ? hit1 : hit2;

      // Test error rates.
      if ( hit1.ErrorRate( ) > MAX_ERROR_RATE ) continue;
      if ( hit2.ErrorRate( ) > MAX_ERROR_RATE ) continue;

      // Test proper aligns.
      if ( PROPER_ONLY ) {
	if ( ! hit1.IsProper( ) ) continue;
	if ( ! hit2.IsProper( ) ) continue;
      }

      // Test right extends properly left.
      if ( EXTEND_ONLY ) {
	if ( right.pos2( ) < left.pos2( ) ) continue;
	if ( left.Pos2( ) > right.Pos2( ) ) continue;
      }

      // Ok, pair is good.
      int pair_id = pairs.GetPairId( hit1.query_id );
      const read_pairing &pair = pairs[ pair_id ];
      out << hit1.target_id << " ["
	  << left.pos2( ) << ","
	  << right.Pos2( ) << ") len="
	  << right.Pos2( ) - left.pos2( ) << " sep="
	  << ObservedSeparation( hit1, hit2 ) << " stretch="
	  << ToString( Stretch( hit1, hit2, pair ), 2 ) << "\n";

      left.PrintParseable( out );
      right.PrintParseable( out );
      out << "\n";
    }
    out.close( );
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;

}
