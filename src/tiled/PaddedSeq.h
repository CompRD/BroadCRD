// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef PADDED_SEQUENCE_H
#define PADDED_SEQUENCE_H

#include "PackAlign.h"
#include "String.h"
#include "Vec.h"
#include "tiled/CharBaseUtil.h"

/**
 * class padded_seq
 *
 * It organizes pads in a given sequence (the type of sequence does not
 * matter). Think of pads as of events happening in a Galilean reference
 * system, so that you may move the sequence back and forth, keeping the
 * informations about gaps.
 *
 * Be extra careful when you read the code in distinguishing local from
 * global (Galilean) coordinates. Pads are in local coordinates: pads_[2]=5
 * means that the fifth base in the sequence is a pad (actually the third
 * pad). However AddPad( pos ) means that if x is such that begin_+x=pos,
 * then the the x-th base in the sequence is a pad.
 *
 * begin_ denotes the start of the the read on the global system; al_begin_
 * (which is >=0) the beginning of the alignment (ie of the significative
 * part).
 *
 * len_ is the length of the aligned portion of sequence.
 */
class padded_seq {
  
public:
  
  padded_seq( );
  
  padded_seq( int begin, int al_begin, int len, const vec<int> *pads = 0 );
  
  void Set( int begin, int al_begin, int len, const vec<int> *pads = 0 );
  
  void ResetBegin( int begin ) { begin_ = begin; }
  
  void ResetAlBegin( int al_begin ) { al_begin_ = al_begin; }

  void ResetUnpaddedLength( int len ) { len_ = len; }
  
  void Clear( );

  void AddPad( int pos );
  
  void RemovePad( int pos );

  void ShiftPad( int pos, int shift );

  void ShiftPadsAfter( int pos, int shift );
  
  void FromAlign( const align &al );
  
  void ToAlign( align &al ) const;

  int AlBegin( ) const { return al_begin_; }

  int Begin( ) const { return begin_; }
  
  int PadsCount( ) const { return pads_.size( ); }

  int UnpaddedLength( ) const { return len_; }
  
  int PaddedLength( ) const { return len_ + pads_.size( ); }
  
  bool IsPad( int pos ) const { return ( this->PadIndex(pos-begin_) > -1 ); }
  
  int Pad( int ii ) const { return pads_[ii]; }

  int ToPaddedPos( int pos ) const;

  int ToUnpaddedPos( int pos ) const;

  char Base( const basevector *fastb, int pos ) const;
  
  int Qual( const qualvector *qualb, int pos ) const;

  void ToVec( const basevector *fastb,
	      const qualvector *qualb,
	      vec<char> &bases,
	      vec<int> &quals,
	      int from = -1,
	      int to = -1 ) const;
  
  void Unpack( const basevector *fastb,
	       const qualvector *qualb,
	       vec<char> &bases,
	       vec<int> &quals,
	       bool filter_ends = true ) const;
  
  void PrintBrief( ostream &out, const bool &new_line = true ) const;

  int &operator[] ( int pos ) { return pads_[pos]; }

  const int &operator[] ( int pos ) const { return pads_[pos]; }

  friend istream &operator>> ( istream &in, padded_seq &pad );

  friend ostream &operator<< ( ostream &out, const padded_seq pad );


private:
  
  void TrimEndingPads( );

  int PadIndex( int pos ) const;
  
  
private:
  
  int begin_;        // begin of sequence on a Galilean reference system
  int al_begin_;     // >=0, it is the start of the alignment
  int len_;          // length of aligned portion (without pads)
  vec<int> pads_;    // sorted pads positions in the sequence
  
};

#endif
