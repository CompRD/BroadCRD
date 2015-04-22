// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "ReadPairing.h"
#include "STLExtensions.h"
#include "lookup/DetectValidPairs.h"
#include "lookup/LookAlign.h"

/**
 * DetectValidPairs
 */
void DetectValidPairs( const vec<look_align_plus> &hits,
		       const vec<read_pairing> &pairs,
		       const vec<int> &to_pair,
		       vec< pair<int,int> > &valid,
		       double max_stretch,
		       double max_mult,
		       int *cap )
{
  valid.clear( );

  // Check hits is sorted.
  order_lookalign_Query sorterQ;
  ForceAssert( is_sorted( hits.begin( ), hits.end( ), sorterQ ) );
  
  // Map (first hits).
  vec<int> fhits( to_pair.size( ), -1 );
  for (int ii=hits.size( )-1; ii>=0; ii--)
    fhits[ hits[ii].query_id ] = ii;
  
  // If cap is given, generate a map to store the number of valid pairs.
  vec<int> count;
  if ( cap )
    count.resize( pairs.size( ), 0 );

  // Loop over all hits.
  for (int id1=0; id1<(int)hits.size( ); id1++) {
    
    // Ids of read and partner.
    int rid = hits[id1].query_id;
    int pair_id = to_pair[rid];
    if ( pair_id < 0 )
      continue;
    const read_pairing &pair = pairs[pair_id];
    int pid = pair.id1 == rid ? pair.id2 : pair.id1;
    
    // Need to cap?
    if ( cap && count[pair_id] > *cap )
      continue;

    // No hits contains partner.
    int fhit = fhits[pid];
    if ( fhit < 0 )
      continue;

    // Loop over all hits containing partner.
    for (int id2=fhit; id2<(int)hits.size( ); id2++) {
      if ( hits[id2].query_id != pid )
	break;

      // Avoid duplication.
      if ( id2 < id1 )
	continue;
      
      // Invalid pair (check which test we want to run).
      const look_align_plus &hit1 = hits[id1];
      const look_align_plus &hit2 = hits[id2];

      bool is_valid = false;
      if ( max_stretch < 0 && max_mult < 0 ) {
	if ( IsLogicalPair( hit1, hit2 ) )
	  is_valid = true;
      }
      else {
	const read_pairing &pair = pairs[pair_id];
	if ( IsValidPair( hit1, hit2, pair, max_stretch, max_mult ) )
	  is_valid = true;
      }
      if ( ! is_valid )
	continue;
      
      // Is valid.
      if ( cap )
	count[pair_id] += 1;
      
      valid.push_back( make_pair( id1, id2 ) );
    }
  }

}

