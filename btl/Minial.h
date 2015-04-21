///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__MINIAL__H
#define BTL__MINIAL__H

#include "lookup/LookAlign.h"

/**
 * struct minial
 *
 * It encapsulates the alignment of a read on a vector in a very
 * compact way. Sign of read_id_ captures orientation. Note the
 * unusual setup: reads (long) are stored as queries (1), and vectors
 * (short) as targets (2).
 */
struct minial {

  minial( ) : read_id_ ( -1 ), vector_id_ ( -1 ), offset_ ( 0 ) { }

  minial( const look_align_plus &hit ) {
    this->SetFromLookAlignPlus( hit );
  }

  void SetFromLookAlignPlus( const look_align_plus &hit ) {
    read_id_ = hit.target_id;   // note that target_id is just a plain int
    vector_id_ = hit.query_id;
    offset_ = hit.a.pos2( ) - hit.a.pos1( );
    if ( hit.rc1 ) read_id_ = - 1 - read_id_;
  }
    
  bool Rc( ) const { return read_id_ < 0; }
  int64_t ReadId( ) const { return read_id_ >= 0 ? read_id_ : -read_id_-1; }
  int Offset( ) const { return offset_;}

  void PrintInfo( ostream &out ) const {
    out << vector_id_ << " aligns"
	<< ( this->Rc( ) ? " rc" : "" ) << " on read "
	<< this->ReadId( ) << " @"
	<< offset_ << "\n";
  }

  friend bool operator== ( const minial &left, const minial &right ) {
    if ( left.read_id_ != right.read_id_ ) return false;
    if ( left.vector_id_ != right.vector_id_ ) return false;
    return left.offset_ == right.offset_;
  }

  friend bool operator< ( const minial &left, const minial &right ) {
    if ( Abs( left.ReadId( ) ) < Abs( right.ReadId( ) ) ) return true;
    if ( Abs( left.ReadId( ) ) > Abs( right.ReadId( ) ) ) return false;
    if ( left.read_id_ > right.read_id_ ) return true;   // fw first...
    if ( left.read_id_ < right.read_id_ ) return false;  // ...then rc
    if ( left.vector_id_ < right.vector_id_ ) return true;
    if ( left.vector_id_ > right.vector_id_ ) return false;
    return left.offset_ < right.offset_;
  }

  friend ostream &operator<< ( ostream &out, const minial &al ) {
    out << al.read_id_ << "\t" << al.vector_id_ << "\t" << al.offset_ << "\n";
    return out;
  }

  friend istream &operator>> ( istream &in, minial &al ) {
    in >> al.read_id_ >> al.vector_id_ >> al.offset_;
    return in;
  }

  int64_t read_id_;  // signed, to capture orientation
  int vector_id_;
  int offset_;

};

#endif
