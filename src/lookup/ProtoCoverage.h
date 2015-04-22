// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Simple minded tools to deal with coverage implied by a set of aligns
// onto a reference genome.

#ifndef PROTO_COVERAGE_H
#define PROTO_COVERAGE_H

#include "lookup/LookAlign.h"

// Build vectors of coverages for targets (use only specified hits) (output
// is saved in covers).
void GenerateCovers( int n_targets,
		     const vec<look_align_plus> &hits,
		     const vec<Bool> keeper,
		     vec< vec<int> > &covers )
{
  covers.clear( );
  covers.resize( n_targets );

  // Resize covers.
  int checked = 0;
  for (int ii=0; ii<(int)hits.size( ); ii++) {
    if ( checked == n_targets )
      break;
    int tid = hits[ii].target_id;
    if ( covers[tid].size( ) < 1 ) {
      covers[tid].resize( hits[ii].target_length, 0 );
      checked++;
    }
  }

  // Eval actual coverage.
  for (int ii=0; ii<(int)hits.size( ); ii++) {
    if ( ! keeper[ii] )
      continue;

    int tid = hits[ii].target_id;
    int tlen = hits[ii].target_length;
    int begin = Max( 0, hits[ii].a.pos2( ) - hits[ii].a.pos1( ) );
    int end = Min( begin + (int)hits[ii].query_length, tlen );

    for (int jj=begin; jj<end; jj++)
      covers[tid][jj] += 1;
  }
}

// What is (roughly) the max coverage of a region an hit aligns to?
int MaxCoverage( const vec< vec<int> > &covers,
		 const look_align_plus &hit )
{
  const vec<int> &cov = covers[hit.target_id];
  int begin = Max( 0, hit.a.pos2( ) - hit.a.pos1( ) );
  int end = Min( begin + (int)hit.query_length, (int)cov.size( ) );
  
  return *max_element( cov.begin( ) + begin, cov.begin( ) + end );
}

#endif
