///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__DEDUP_WINNER_H
#define BTL__DEDUP_WINNER_H

#include "Equiv.h"
#include "Qualvector.h"

/**
 * DedupWinner
 *
 * Pick a winner from a set of duplicates. It returns -1 on error.
 */
int DedupWinner( const vecqvec &quals, const vec<int> &orbit )
{
  if ( orbit.size( ) < 1 ) return -1;
  
  vec< triple<int,double,int> > len2qual;
  len2qual.reserve( orbit.size( ) );
  for (int ii=0; ii<orbit.isize( ); ii++) {
    int id = orbit[ii];
    int len = quals[id].size( );
    double qual = Mean( quals[id] );
    len2qual.push_back( triple<int,double,int>( len, qual, id ) );
  }
  sort( len2qual.rbegin( ), len2qual.rend( ) );
  
  return len2qual[0].third;
}

#endif
