///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_LIB_LG_H
#define C_LIB_LG_H

#include "PairsManager.h"
#include "VecString.h"
#include "Vec.h"
#include "paths/Alignlet.h"

/**
 * class CLibLG
 *
 * It deals with stats for an ALLPATH's library.
 */
class CLibLG {

public:

  CLibLG( ) :
    max_stretch_ ( 3.5 ), pairs_ ( 0 ), aligns_ ( 0 ), index_ ( 0 ) { }
  
  CLibLG( const int &id, const String &name ) :
    max_stretch_ ( 3.5 ), pairs_ ( 0 ), aligns_ ( 0 ), index_ ( 0 ) {
    this->SetLib( id, name ); }
  
  // Set core info.
  void SetLib( const int &id, const String &name ) {
    lib_id_ = id;
    lib_name_ = name;
    this->ResetCounters( );
  }

  // Set pointers.
  void SetPointers( const PairsManager *pairs,
		    const vec<alignlet> *aligns,
		    const vec<int> *index ) {
    pairs_ = pairs;
    aligns_ = aligns;
    index_ = index;
  }
  
  // Update counters for this pair.
  void AddPairInfo( int pair_id ) {
    ForceAssert( pairs_ );
    ForceAssertEq( pairs_->libraryID( pair_id ), lib_id_ );
    
    // Total.
    n_tot_ += 1;

    // End reads align uniquely.
    int id1 = pairs_->ID1( pair_id );
    int id2 = pairs_->ID2( pair_id );
    if ( (*index_)[id1] < 0 || (*index_)[id2] < 0 )
      return;
    n_11_ += 1;

    // In different supers.
    const alignlet &al1 = (*aligns_)[ (*index_)[id1] ];
    const alignlet &al2 = (*aligns_)[ (*index_)[id2] ];
    if ( al1.TargetId( ) != al2.TargetId( ) ) {
      n_sep_ += 1;
      return;
    }

    // In the same super. Valid?
    if ( al1.Fw1( ) == al2.Fw1( ) ) return;
    int fw_end = al1.Fw1( ) ? al1.Pos2( ) : al2.Pos2( );
    int rc_beg = al1.Fw1( ) ? al2.pos2( ) : al1.pos2( );
    int observed_sep = rc_beg - fw_end;
    int given_sep = pairs_->sep( pair_id );
    int given_sd = pairs_->sd( pair_id );
    if ( given_sd > 0 ) {
      double stretch = SafeQuotient( observed_sep - given_sep, given_sd );
      if ( - max_stretch_ <= stretch && stretch <= max_stretch_ )
	n_val_ += 1;
    }
  }

  // Build a line for a table entry (or a legend line).
  vec<String> TableLine( bool legend = false ) const {
    vec<String> line;

    if ( legend ) line.push_back( "library name" );
    else line.push_back( lib_name_ );

    if ( legend ) line.push_back( "total" );
    else line.push_back( ToString( n_tot_ ) );
    
    if ( legend ) line.push_back( "n11" );
    else line.push_back( ToString( n_11_ ) );
    
    double ratio = -1.0;
    if ( n_tot_ > 0 ) ratio = 100.0 * SafeQuotient( n_11_, n_tot_ );
    if ( legend ) line.push_back( "% n11/tot" );
    else ( line.push_back( ratio < 0 ? "na" : ToString( ratio, 2 ) ) );

    if ( legend ) line.push_back( "valid" );
    else line.push_back( ToString( n_val_ ) );

    ratio = -1.0;
    if ( n_11_ > 0 ) ratio = 100.0 * SafeQuotient( n_val_, n_11_ );
    if ( legend ) line.push_back( "% val/n11" );
    else ( line.push_back( ratio < 0 ? "na" : ToString( ratio, 2 ) ) );
    
    int n_inv = n_11_ - n_val_ - n_sep_;
    if ( legend ) line.push_back( "invalid" );
    else line.push_back( ToString( n_inv ) );

    ratio = -1.0;
    if ( n_11_ > 0 ) ratio = 100.0 * SafeQuotient( n_inv, n_11_ );
    if ( legend ) line.push_back( "% inv/n11" );
    else ( line.push_back( ratio < 0 ? "na" : ToString( ratio, 2 ) ) );
    
    if ( legend ) line.push_back( "separated" );
    else line.push_back( ToString( n_sep_ ) );

    ratio = -1.0;
    if ( n_11_ > 0 ) ratio = 100.0 * SafeQuotient( n_sep_, n_11_ );
    if ( legend ) line.push_back( "% sep/n11" );
    else ( line.push_back( ratio < 0 ? "na" : ToString( ratio, 2 ) ) );

    return line;
  }
  
  
private:

  void ResetCounters( ) {
    n_tot_ = 0;
    n_11_ = 0;
    n_val_ = 0;
    n_sep_ = 0;
  }
  
  
private:

  double max_stretch_; // defines valid inserts

  const PairsManager *pairs_;    // pairs manager
  const vec<alignlet> *aligns_;  // alignments of reads onto assembly
  const vec<int> *index_;        // -2: no align; -1: multiple; >=0: index

  int lib_id_;       // id in the pairs manager
  String lib_name_;  // library name

  longlong n_tot_;   // total number of pairs
  longlong n_11_;    // both end reads align uniquely
  longlong n_val_;   // reads align uniquely and validly on the same super
  longlong n_sep_;   // reads align uniquely on different supers
  
};

#endif
