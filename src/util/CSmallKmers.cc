///////////////////////////////////////////////////////////////////////////////
//                   Software COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Intvector.h"
#include "util/CSmallKmers.h"

/**
 * class CSmallKmers
 * Constructor
 */
CSmallKmers::CSmallKmers( const int k, const vecbvec *target ) :
  k_ ( k )
{
  ForceAssertLt( k_, 13 );
  this->Setup( );
  if ( target ) this->SetTargetBases( *target );
}

/**
 * class CSmallKmers
 * SetTargetBases
 */
void CSmallKmers::SetTargetBases( const vecbvec &target ) {
  target_ = target;
  target_ids_.resize( fourPow_[k_] );
  target_pos_.resize( fourPow_[k_] );
  for (size_t tid=0; tid<target.size( ); tid++) {
    for (size_t pos=0; pos<target[tid].size( )-k_; pos++) {
      int id = this->ToInt( tid, pos );
      target_ids_[id].push_back( tid );
      target_pos_[id].push_back( pos );
    }
  }
}

/**
 * class CSmallKmers
 * ToInt
 */
int CSmallKmers::ToInt( const bvec &bases, const int pos ) const
{
  int id = 0;
  for (int ii=pos; ii<pos+k_; ii++)
    id += bases[ii] * fourPow_[ii-pos];
  return id;
}

/**
 * class CSmallKmers
 * ToInt
 */
int CSmallKmers::ToInt( const int tid, const int pos ) const
{
  return this->ToInt( target_[tid], pos );
}

/**
 * class CSmallKmers
 * ToInt
 */
void CSmallKmers::ToBases( int id, bvec &bases ) const
{
  bases.resize( k_ );
  for (int ii=0; ii<k_; ii++) {
    bases.Set(ii,id&3);
    id >>= 2;
  }
}

/**
 * class CSmallKmers
 * NPlacements
 */
size_t CSmallKmers::NPlacements( const int id ) const
{
  return target_ids_[id].size( );
}

/**
 * class CSmallKmers
 * TargetId
 */
int64_t CSmallKmers::TargetId( const int id, const int ii ) const
{
  return target_ids_[id][ii];
}

/**
 * class CSmallKmers
 * TargetPos
 */
int CSmallKmers::TargetPos( const int id, const int ii ) const
{
  return target_pos_[id][ii];
}

/**
 * class CSmallKmers
 * Setup
 * private
 */
void CSmallKmers::Setup( )
{
  fourPow_.resize( k_+1, 1 );
  for (int ii=1; ii<k_+1; ii++)
    fourPow_[ii] = fourPow_[ii-1] * 4;
}

