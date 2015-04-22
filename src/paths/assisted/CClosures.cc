///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "paths/assisted/CClosures.h"

/**
 * CClosures
 * Clear
 */
void CClosures::Clear( )
{
  bases_.clear( );
  overlaps_.clear( );
  ctypes_.clear( );
}

/**
 * CClosures
 * Reserve
 */
void CClosures::Reserve( size_t nclosures )
{
  bases_.reserve( nclosures );
  overlaps_.reserve( nclosures );
  ctypes_.reserve( nclosures );
}

/**
 * CClosures
 * Resize
 */
void CClosures::Resize( size_t nclosures )
{
  bases_.resize( nclosures );
  overlaps_.resize( nclosures );
  ctypes_.resize( nclosures );
}

/**
 * CClosures
 * AddOverlap
 */
void CClosures::AddOverlap( const int over )
{
  this->AddClosure( bvec( ), over, OVERLAP );
}

/**
 * CClosures
 * AddTrueClosure
 */
void CClosures::AddTrueClosure( const bvec &bases, const closure_t ct )
{
  this->AddClosure( bases, -1, ct );
}

/**
 * CClosures
 * AddClosure
 */
void CClosures::AddClosure( const bvec &bases,
			    const int over,
			    const closure_t ct )
{
  if ( over == -1 ) {
    ForceAssert( ct != OVERLAP );
  }
  else {
    ForceAssert( over > 0 );
    ForceAssert( ct == OVERLAP );
  }

  bases_.push_back( bases );
  overlaps_.push_back( over );
  ctypes_.push_back( ct );
}
  
/**
 * CClosures
 * AddClosure
 */
void CClosures::AddClosure( const CClosures &other, const int ii )
{
  this->AddClosure( other.Bases( ii ), other.Overlap( ii ), other.Ctype( ii ) );
}

/**
 * CClosures
 * ChompAdjacencies
 *
 * Remove K-1 bases on the left and/or right of closures. This may
 * cause a true closure to be turned into an overlap.
 */
void CClosures::ChompAdjacencies( const int K )
{
  CClosures closures_out;
  closures_out.Reserve( bases_.size( ) );

  for (size_t ii=0; ii<bases_.size( ); ii++) {
    if ( ctypes_[ii] == OVERLAP ) {
      closures_out.AddClosure( bases_[ii], overlaps_[ii], ctypes_[ii] );
      continue;
    }
    
    int start = 0;
    int len = bases_[ii].size( );
    if ( ctypes_[ii] == LEFT_CLOSURE || ctypes_[ii] == FULL_CLOSURE )
      start = K-1; 
    if ( ctypes_[ii] == RIGHT_CLOSURE || ctypes_[ii] == FULL_CLOSURE )
      len = len - ( K-1 ) - start;
    
    // A true closure (of size >= 0).
    if ( len >=0 ) {
      bvec newb = bvec( bases_[ii], start, len );
      closures_out.AddClosure( newb, overlaps_[ii], ctypes_[ii] );
      continue;
    }
    
    // Turn closure into overlap.
    closures_out.AddOverlap( -len );
  }

  swap( closures_out, *this );
  
}

/**
 * CClosures
 * SortClosures
 *
 * NB: you must run ChompAdjacencies first.
 *
 * WARNING! As closures are sorted by length (among other things), the
 * final output is not completely deterministic! This should be fixed.
 */
void CClosures::GapBasedSort( int gap_size )
{
  vec<int> sorted_ids;
  sorted_ids.reserve( bases_.size( ) );

  // First, overlaps and full closures.
  vec< pair<int,int> > dist2ids;
  dist2ids.reserve( bases_.size( ) );
  for (int ii=0; ii<(int)bases_.size( ); ii++) {
    if ( ctypes_[ii] == LEFT_CLOSURE || ctypes_[ii] == RIGHT_CLOSURE )
      continue;
    
    int implied_gap = 0;
    if ( ctypes_[ii] == OVERLAP ) implied_gap = - overlaps_[ii];
    else implied_gap = bases_[ii].size( );
    int dist = Abs( implied_gap - gap_size );
    dist2ids.push_back( make_pair( dist, ii ) );
  }
  sort( dist2ids.begin( ), dist2ids.end( ) );

  for (int ii=0; ii<dist2ids.isize( ); ii++)
    sorted_ids.push_back( dist2ids[ii].second );

  // Left closures.
  vec< pair<int,int> > lens2ids;
  for (int ii=0; ii<(int)bases_.size( ); ii++)
    if ( ctypes_[ii] == LEFT_CLOSURE )
      lens2ids.push_back( make_pair( (int)bases_[ii].size( ), ii ) );
  sort( lens2ids.begin( ), lens2ids.end( ) );
  for (int ii=0; ii<lens2ids.isize( ); ii++)
    sorted_ids.push_back( lens2ids[ii].second );

  // Right closures.
  lens2ids.clear( );
  for (int ii=0; ii<(int)bases_.size( ); ii++)
    if ( ctypes_[ii] == RIGHT_CLOSURE )
      lens2ids.push_back( make_pair( (int)bases_[ii].size( ), ii ) );
  sort( lens2ids.begin( ), lens2ids.end( ) );
  for (int ii=0; ii<lens2ids.isize( ); ii++)
    sorted_ids.push_back( lens2ids[ii].second );

  // This should not be needed...
  ForceAssert( sorted_ids.size( ) == bases_.size( ) );

  // Build new data, swap and done.
  vecbvec new_bases;
  vec<int> new_overlaps;
  vec<closure_t> new_ctypes;
  for (int ii=0; ii<sorted_ids.isize( ); ii++) {
    int new_id = sorted_ids[ii];
    new_bases.push_back( bases_[new_id] );
    new_overlaps.push_back( overlaps_[new_id] );
    new_ctypes.push_back( ctypes_[new_id] );
  }
  swap( bases_, new_bases );
  swap( overlaps_, new_overlaps );
  swap( ctypes_, new_ctypes );
  
}

/**
 * CClosures
 * PrintInfo
 */
void CClosures::PrintInfo( int ii,
			   int super_id,
			   int gap_id, 
			   ostream &out ) const
{
  String str_over = "0";
  String str_type = "";
  if ( ctypes_[ii] == OVERLAP ) {
    str_type = "O";
    str_over = ToString( overlaps_[ii] );
  }
  else if ( ctypes_[ii] == FULL_CLOSURE ) str_type = "F";
  else if ( ctypes_[ii] == LEFT_CLOSURE ) str_type = "L";
  else if ( ctypes_[ii] == RIGHT_CLOSURE ) str_type = "R";
  else ForceAssert( 1 == 0 );
  
  out << super_id << "\t"
      << gap_id << "\t"
      << str_over << "\t" 
      << str_type << "\n";
}

/**
 * CClosures
 * PrintReport
 */
void CClosures::PrintReport( ostream &out ) const
{
  int n_overlaps = 0;
  int n_full = 0;
  int n_left = 0;
  int n_right = 0;
  for (size_t ii=0; ii<bases_.size( ); ii++) {
    if ( ctypes_[ii] == OVERLAP ) n_overlaps++;
    else if ( ctypes_[ii] == FULL_CLOSURE ) n_full++;
    else if ( ctypes_[ii] == LEFT_CLOSURE ) n_left++;
    else if ( ctypes_[ii] == RIGHT_CLOSURE ) n_right++;
    else ( ForceAssert( 1 == 0 ) );
  }

  if ( bases_.size( ) < 1 ) {
    out << "no closures found\n";
    return;
  }

  out << "found " << bases_.size( ) << " closures:   ";
  if ( n_overlaps > 0 ) out << n_overlaps << " overlaps   ";
  if ( n_full > 0 ) out << n_full << " full   ";
  if ( n_left > 0 ) out << n_left << " left   ";
  if ( n_right > 0 ) out << n_right << " right   ";
  out << "\n";
  
}
