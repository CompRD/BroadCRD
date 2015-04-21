///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

#include "Basevector.h"
#include "SeqInterval.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "util/RunCommand.h"

void BuildRanges( const int MIN_FLEN,
		  const int MAX_FLEN,
		  const int npairs,
		  const vec<look_align> &als1,
		  const vec<look_align> &als2,
		  vec<seq_interval> &sis,
		  ostream *log = 0 )
{
  sis.clear( );

  // Mapping read ids to align ids.
  if ( log ) *log << Date( ) << ": mapping reads to aligns" << endl;
  vec<int> to_al1( npairs, -1 );
  for (int aid=0; aid<als1.isize( ); aid++) {
    int rid = als1[aid].query_id;
    if ( to_al1[rid] == -1 ) to_al1[rid] = aid;
    else to_al1[rid] = -2;
  }
  vec<int> to_al2( npairs, -1 );
  for (int aid=0; aid<als2.isize( ); aid++) {
    int rid = als2[aid].query_id;
    if ( to_al2[rid] == -1 ) to_al2[rid] = aid;
    else to_al2[rid] = -2;
  }

  // Loop over all pairs, create seq_intervals of fragment pairs.
  sis.reserve( npairs );
  if ( log ) *log << Date( ) << ": mapping fragment pairs to reference" << endl;
  for (int pid=0; pid<npairs; pid++) {
    
    // No aligns, or multiple aligns found.
    if ( to_al1[pid] < 0 || to_al2[pid] < 0 ) continue;

    // On different targets.
    const look_align &al1 = als1[ to_al1[pid] ];
    const look_align &al2 = als2[ to_al2[pid] ];
    if ( al1.target_id != al2.target_id ) continue;

    // Invalid orientation.
    if ( al1.Fw1( ) == al2.Fw1( ) ) continue;
    
    // Observed fragment length is out of boundary.
    const look_align &alF = al1.Fw1( ) ? al1 : al2;
    const look_align &alR = al1.Fw1( ) ? al2 : al1;
    int flen = alR.EndOnTarget( ) - alF.StartOnTarget( );
    if ( flen < MIN_FLEN || flen > MAX_FLEN ) continue;

    // A keeper.
    int tid = al1.target_id;
    int beg = alF.StartOnTarget( );
    int end = alR.EndOnTarget( );
    sis.push_back( seq_interval( pid, tid, beg, end ) );
  }
  const int nsis = sis.size( );
  const String str_nsis = ToStringAddCommas( nsis );

  // No ranges found?
  if ( nsis < 1 ) {
    if ( log ) *log << Date( ) << ": no ranges found" << endl;
    return;
  }

  // Sort seq_intervals.
  if ( log ) *log << Date( ) << ": sorting " << str_nsis << " ranges" << endl;
  sort( sis.begin( ), sis.end( ) );
  
  // Remove duplicates.
  if ( log ) *log << Date( ) << ": removing duplicates... " << flush;
  int ndups = 0;
  vec<seq_interval> dedupped;
  dedupped.reserve( nsis );
  dedupped.push_back( sis[0] );
  for (int ii=1; ii<nsis; ii++) {
    if ( sis[ii].SeqId( ) == dedupped.back( ).SeqId( ) &&
	 sis[ii].Begin( ) == dedupped.back( ).Begin( ) &&
	 sis[ii].End( ) == dedupped.back( ).End( ) ) {
      ndups++;
      continue;
    }
    dedupped.push_back( sis[ii] );
  }
  if ( log ) *log << ToStringAddCommas( ndups ) << " found and removed" << endl;
  swap( dedupped, sis );
  
}

/**
 * LocalizeIndexedReads
 *
 * Localize indexed reads, based on their alignments onto a reference.
 * Not robust for more than 2^32 reads.
 */
int main( int argc, char* argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( IN_DIR );
  CommandArgument_String( REF_FASTB );
  CommandArgument_String_OrDefault( HEAD, "reads" );

  // Min overlap size.
  CommandArgument_Int_OrDefault( K, 24 );
  
  // Read length (used to compute read coverage).
  CommandArgument_Int_OrDefault( RLEN, 150 );
  
  // Min and max expected fragment length.
  CommandArgument_Int_OrDefault( MIN_FLEN, 150 );
  CommandArgument_Int_OrDefault( MAX_FLEN, 1500 );

  // Min coverage to tag a region as "assemblable".
  CommandArgument_Int_OrDefault( MIN_COV, 10 );
  
  // Do not use cached data, if FORCE=True.
  CommandArgument_Bool_OrDefault( FORCE, True );

  EndCommandArguments;

  // Dir and file names.
  String head = IN_DIR + "/" + HEAD;
  String head1 = head + "R1";
  String head2 = head + "R2";

  String bases1_file = head1 + ".fastb";
  String aligns1_file = head1 + ".qlt";

  String bases2_file = head2 + ".fastb";
  String aligns2_file = head2 + ".qlt";

  String ranges_file = head + ".ranges";
  
  vec<String> needed;
  needed.push_back( REF_FASTB );
  needed.push_back( bases1_file );
  needed.push_back( aligns1_file );
  needed.push_back( bases2_file );
  needed.push_back( aligns2_file );

  // Load aligns.
  int npairs = MastervecFileObjectCount( bases1_file );
  int npairs2 = MastervecFileObjectCount( bases2_file );
  ForceAssertEq( npairs, npairs2 );

  const String str_npairs = ToStringAddCommas( npairs );
  cout << Date( ) << ": " << str_npairs << " pairs found" << endl;
  
  // Build ranges.
  vec<seq_interval> ranges;
  if ( FORCE || ! IsRegularFile( ranges_file ) ) {
    cout << Date( ) << ": loading R1 aligns... " << flush;
    vec<look_align> als1;
    LoadLookAligns( aligns1_file, als1 );
    cout << ToStringAddCommas( als1.size( ) ) << " found" << endl;
    
    cout << Date( ) << ": loading R2 aligns... " << flush;
    vec<look_align> als2;
    LoadLookAligns( aligns2_file, als2 );
    cout << ToStringAddCommas( als2.size( ) ) << " found" << endl;
    
    BuildRanges( MIN_FLEN, MAX_FLEN, npairs, als1, als2, ranges, &cout );
    
    cout << Date( ) << ": saving ranges" << endl;
    ofstream out( ranges_file.c_str( ) );
    out << ranges.size( ) << "\n";
    for (int ii=0; ii<ranges.isize( ); ii++)
      out << ranges[ii] << "\n";
    out.close( );
  }
  else {
    cout << Date( ) << ": loading cached ranges" << endl;
    ifstream in( ranges_file.c_str( ) );
    int nranges;
    in >> nranges;
    ranges.reserve( nranges );
    while ( in ) {
      seq_interval sis;
      in >> sis;
      if ( !in ) break;
      ranges.push_back( sis );
    }
  }
  
  // Load reference sizes.
  cout << Date( ) << ": loading reference sizes" << endl;
  vec<int> rlens;
  {
    vecbvec ref( REF_FASTB );
    rlens.reserve( ref.size( ) );
    for (int ii=0; ii<(int)ref.size( ); ii++)
      rlens.push_back( (int)ref[ii].size( ) );
  }
  
  // Create bags.
  vec<int> bag_starts( 1, 0 );
  
  for (int ii=0; ii<ranges.isize( ); ii++) {
    const seq_interval &curr = ranges[ii];
    int cid = curr.SeqId( );
    int cbeg = curr.Begin( );
    int cend = curr.End( );
    int gap = curr.Length( ) - 2 * RLEN;
    if ( ii < ranges.isize( )-1 ) {
      const seq_interval &next = ranges[ii+1];
      int nid = next.SeqId( );
      int nbeg = next.Begin( );
      if ( nid != cid || nbeg - cend > -K ) bag_starts.push_back( ii+1 );
    }
  }

  // Done loading.
  cout << endl;
  
  // Counters.
  const int nbags = bag_starts.size( );
  vec<int> ntigs_in_bag;
  vec<int> pc_lens;
  ntigs_in_bag.resize( nbags, 0 );

  // Loop over all bags.
  for (int bid=0; bid<nbags; bid++) {
    const int rbeg = bag_starts[bid];
    const int rend = ( bid < nbags-1 ) ? bag_starts[bid+1] : ranges.isize( );
    const int bag_size = rend - rbeg;
    
    // Turn ranges back into pseudo reads (preads).
    vec<ho_interval> preads;
    preads.reserve( 2 * bag_size );
    for (int ii=rbeg; ii<rend; ii++) {
      const seq_interval &range = ranges[ii];
      preads.push_back( ho_interval( range.Begin( ), range.Begin( ) + RLEN ) );
      preads.push_back( ho_interval( range.End( ) - RLEN, range.End( ) ) );
    }
    sort( preads.begin( ), preads.end( ) );
    preads.erase( unique( preads.begin( ), preads.end( ) ), preads.end( ) );
    ForceAssertGt( preads.isize( ), 0 );

    // Merge preads into pseudo contigs (pcontigs).
    vec<ho_interval> pcontigs;
    pcontigs.reserve( preads.size( ) );
    pcontigs.push_back( preads[0] );
    for (int ii=1; ii<preads.isize( ); ii++) {
      if ( Overlap( pcontigs.back( ), preads[ii] ) >= K )
	pcontigs[ pcontigs.size( )-1 ].Merge( preads[ii] );
      else
	pcontigs.push_back( preads[ii] );
    }
    ntigs_in_bag[bid] += pcontigs.isize( );

    for (int ii=0; ii<pcontigs.isize( ); ii++)
      pc_lens.push_back( pcontigs[ii].Length( ) );
  }

  // Sort counters.
  sort( ntigs_in_bag.begin( ), ntigs_in_bag.end( ) );
  sort( pc_lens.begin( ), pc_lens.end( ) );
  
  // Prepare table (report).
  vec< vec<String> > table;

  table.push_back( MkVec( ToStringAddCommas( nbags ),
			  String( "number of bags" ) ) );
  
  table.push_back( MkVec( ToStringAddCommas( Median( ntigs_in_bag ) ),
			  String( "median number of contigs in each bag" ) ) );
  
  table.push_back( MkVec( ToString( Mean( ntigs_in_bag ), 1 ),
			  String( "mean number of contigs in each bag" ) ) );
  
  table.push_back( MkVec( ToStringAddCommas( N50( ntigs_in_bag ) ),
			  String( "N50 of number of contigs in each bag" ) ) );

  table.push_back( MkVec( ToStringAddCommas( pc_lens.size( ) ),
			  String( "number of contigs" ) ) );
  
  table.push_back( MkVec( ToStringAddCommas( BigSum( pc_lens ) ),
			  String( "total contig length" ) ) );

  table.push_back( MkVec( ToStringAddCommas( N50( pc_lens ) ),
			  String( "N50 length of contigs" ) ) );

  PrintTabular( cout, table, 3, "rl" );
  cout << endl;
  
  // Done.
  cout << Date( ) << ": LocalizeIndexedReads done" << endl;

}
