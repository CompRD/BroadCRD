/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Fastavector.h"
#include "Superb.h"
#include "lookup/LookAlign.h"
#include "paths/SaveScaffoldGraph.h"
#include "util/ScaffoldContigsOnRef.h"

/**
 * ScaffoldContigsOnRef
 */
void ScaffoldContigsOnRef( const vec<superb> &initial_supers,
			   const vec<fastavector> &initial_contigs,
			   const vec<look_align> &aligns,
			   const String &ASSEMBLY_OUT,
			   const int MAX_GAP,
			   const int MAX_OVERLAP,
			   const int MIN_GAP_DEV,
			   ostream *log )
{
  // Make a local copy of the initial contigs (these may be rc-ed).
  vec<fastavector> contigs = initial_contigs;

  // Flag which contigs are actually in the inital supers
  vec<bool> initial( contigs.size( ), false );
  for (int s = 0; s < initial_supers.isize(); ++s) {
    for (int t = 0; t < initial_supers[s].Ntigs(); ++t) {
      initial[initial_supers[s].Tig(t)] = True;
    }      
  }

  // Place contigs on reference (rc if needed).
  vec<bool> placed( contigs.size( ), false );
  vec<seq_interval> placs;
  for (size_t ii=0; ii<aligns.size( ); ii++) {
    const look_align_plus &al = aligns[ii];
    int cid = al.query_id;
    int tid = al.target_id;
    int clen = al.query_length;
    int tlen = al.target_length;

    // Skip align.
    if ( ! initial[cid] ) continue;   // skip contig
    if ( placed[cid] ) continue;   // already placed
    if ( ! al.IsProper( ) ) continue;   // improper
    
    // Place contig (eventually rc it).
    placed[cid] = true;
    placs.push_back( seq_interval( cid, tid, al.a.pos2( ), al.a.Pos2( ) ) );
    if ( al.Rc1( ) ) contigs[cid].ReverseComplement( );
  }

  sort( placs.begin( ), placs.end( ) );
  
  // Create scaffolds (placed contigs).
  if ( log ) *log << Date( ) << ": creating scaffolds" << endl;
  vec<superb> supers;
  for (size_t ii=0; ii<placs.size( ); ii++) {
    int cid = placs[ii].IntervalId( );
    int clen = contigs[cid].size( );

    // New super (first contig, or contigs in different supers).
    if ( ii == 0 || placs[ii].SeqId( ) != placs[ii-1].SeqId( )) {
      superb sup;
      sup.PlaceFirstTig( cid, clen );
      supers.push_back( sup );
      continue;
    }

    // Append contig.
    int gap = placs[ii].Begin( ) - placs[ii-1].End( );
    int dev = Max( MIN_GAP_DEV, gap / 4 );
    if ( gap >= -MAX_OVERLAP && gap <= MAX_GAP ) {
      supers[supers.size( )-1].AppendTig( cid, clen, gap, dev );
      continue;
    }
    
    // New super (gap out of allowed sizes).
    superb sup;
    sup.PlaceFirstTig( cid, clen );
    supers.push_back( sup );
  }

  // Add unplaced contigs.
  if ( log ) *log << Date( ) << ": adding unplaced contigs" << endl;
  for (int cid=0; cid<(int)contigs.size( ); cid++) {
    if ( ! initial[cid] ) continue;
    if ( placed[cid]  ) continue;

    int clen = contigs[cid].size( );
    superb sup;
    sup.PlaceFirstTig( cid, clen );
    supers.push_back( sup );    
  }

  // Save.
  if ( log ) *log << Date( ) << ": saving assembly" << endl;
  SaveScaffoldAssembly( ASSEMBLY_OUT, supers, contigs, NULL, true );
  
}

