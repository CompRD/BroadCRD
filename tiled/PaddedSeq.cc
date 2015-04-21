// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "PackAlign.h"
#include "Vec.h"
#include "tiled/CharBaseUtil.h"
#include "tiled/PaddedSeq.h"



/*
 * padded_seq
 * Constructor
 */
padded_seq::padded_seq( ) :
  begin_ ( 0 ),
  al_begin_ ( 0 ),
  len_ ( 0 )
{ }



/*
 * padded_seq
 * Constructor
 */
padded_seq::padded_seq( int begin,
			int al_begin,
			int len,
			const vec<int> *pads ) :
  begin_ ( begin ),
  al_begin_ ( al_begin ),
  len_ ( len )
{
  if ( pads )
    pads_ = *pads;
}



/*
 * padded_seq
 * Set
 */
void padded_seq::Set( int begin, int al_begin, int len, const vec<int> *pads )
{
  begin_ = begin;
  al_begin_ = al_begin;
  len_ = len;

  if ( pads )
    pads_ = *pads;
}



/*
 * padded_seq
 * Clear
 */
void padded_seq::Clear( )
{
  begin_ = 0;
  al_begin_ = 0;
  len_ = 0;
  pads_.resize( 0 );
}



/*
 * padded_seq
 * AddPad
 */
void padded_seq::AddPad( int pos )
{
  // Reserve memory.
  if ( (int)pads_.capacity( ) < len_ )
    pads_.reserve( len_ );

  // Translate pos in local coordinates.
  pos += -begin_;

  if ( pos > len_ + al_begin_ + (int)pads_.size( ) - 1 )
    return;
  
  if ( pos < 1 + al_begin_ ) {
    begin_++;
    return;
  }
  
  // Insert pad.
  vec<int>::iterator iter = lower_bound( pads_.begin( ), pads_.end( ), pos );
  for (int ii=iter-pads_.begin( ); ii<(int)pads_.size( ); ii++)
    pads_[ii] +=1;
  pads_.insert( iter, pos );

  this->TrimEndingPads( );
}



/*
 * padded_seq
 * RemovePad
 */
void padded_seq::RemovePad( int pos )
{
  // Translate pos in local coordinates.
  pos += -begin_;

  if ( pos > len_ + al_begin_ + (int)pads_.size( ) - 1 )
    return;
  
  if ( pos < 1 + al_begin_ ) {
    begin_--;
    return;
  }
  
  // Find and remove pad.
  int index = this->PadIndex( pos );
  ForceAssert( index > -1 );
  for (int ii=index+1; ii<(int)pads_.size( ); ii++)
    pads_[ii] -= 1;
  pads_.erase( pads_.begin() + index );

  this->TrimEndingPads( );
}



/*
 * padded_seq
 * ShiftPad
 */
void padded_seq::ShiftPad( int pos, int shift )
{
  // Translate pos in local coordinates.
  pos += -begin_;

  // Find and shift pad.
  int index = this->PadIndex( pos );
  ForceAssert( index > -1 );
  pads_[index] += shift;

  this->TrimEndingPads( );
}



/*
 * padded_seq
 * ShiftPadsAfter
 */
void padded_seq::ShiftPadsAfter( int pos, int shift )
{
  // Translate pos in local coordinates.
  pos += -begin_;
  
  for (int ii=0; ii<(int)pads_.size( ); ii++)
    if ( pads_[ii] > pos )
      pads_[ii] += shift;

  this->TrimEndingPads( );
}



/*
 * padded_seq
 * FromAlign
 *
 * Inverse to ToAlign.
 */
void padded_seq::FromAlign( const align &al )
{
  ForceAssert( al.Nblocks( ) > 0 );
  
  // begin_, al_begin_, and len_.
  al_begin_ = al.pos1( );
  begin_ = al.pos2( );
  len_ = 0;
  for (int ii=0; ii<(int)al.Nblocks( ); ii++)
    len_ += al.Lengths(ii);
  
  // pads_.
  int n_pads = 0;
  for (int ii=0; ii<(int)al.Nblocks( ); ii++)
    n_pads += al.Gaps(ii);
  pads_.resize( 0 );
  pads_.reserve( n_pads );
  int pos = al_begin_;
  for (int ii=0; ii<(int)al.Nblocks( ); ii++) {
    for (int jj=0; jj<(int)al.Gaps(ii); jj++) {
      pads_.push_back( pos );
      pos++;
    }
    pos += al.Lengths(ii);
  }

}



/*
 * padded_seq
 * ToAlign
 *
 * Generates an alignments with pos1=al_begin_; pos2=begin_, and gaps(0)=0.
 */
void padded_seq::ToAlign( align &al ) const
{
  int pos1 = al_begin_;
  int pos2 = begin_;
  avector<int> gaps( 0 );
  avector<int> lens( 0 );
  
  // There are no pads.
  if ( pads_.size( ) == 0 ) {
    gaps.Append( 0 );
    lens.Append( len_ );
    al.Set( pos1, pos2, gaps, lens );
    return;
  }
  
  // There are pads.
  gaps.Append( 0 );
  lens.Append( pads_[0] - al_begin_ );
  for (int ii=0; ii<(int)pads_.size( ); ii++) {
    
    int gap = 1;
    while ( ii < (int)pads_.size( ) - 1 && pads_[ii+1] == pads_[ii] + 1 ) {
      gap++;
      ii++;
    }
    
    int len
      = ( ii < (int)pads_.size( ) - 1 )
      ? pads_[ii+1] - pads_[ii] - 1
      : al_begin_ + this->PaddedLength( ) - pads_[ii] - 1;

    gaps.Append( gap );
    lens.Append( len );
  }
  
  // Set al.
  al.Set( pos1, pos2, gaps, lens );
  
}



/*
 * padded_seq
 * ToPaddedPos
 *
 * Translate (unpadded) given pos to its padded value.
 */
int padded_seq::ToPaddedPos( int pos ) const
{
  int current_pos = 0;
  for (int jj=0; jj<this->PaddedLength( ); jj++) {
    if ( this->IsPad( jj ) )
      continue;
    if ( current_pos == pos )
      return jj;
    current_pos++;
  }
  if ( current_pos == pos )
    return this->PaddedLength( );

  // Should never get here.
  FatalErr("ERROR in PaddedSeq::ToPaddedPos(int pos).\nThe input value of pos (you gave " << pos << ") must not be greater than this PaddedSeq's unpadded length (which is " << len_ << ").\n");
  return -666;
}



/*
 * padded_seq
 * ToUnpaddedPos
 *
 * Translate (padded) given pos to its unpadded value.
 */
int padded_seq::ToUnpaddedPos( int pos ) const
{
  int npads = 0;
  for (int ii=0; ii<(int)pads_.size( ); ii++)
    if ( pads_[ii] >= pos )
      break;
  else
    npads++;

  return pos - npads;
}



/*
 * padded_seq
 * Base
 *
 * It returns the base (as a char) at pos. For example, if fastb is AGA**CC,
 * then Base(3)=*, Base(4)=*, and Base(5)=C.
 */
char padded_seq::Base( const basevector *fastb, int pos ) const
{
  if ( pos < 0 || pos >= (int)fastb->size( ) + this->PadsCount( ) )
    return empty_base;
  if ( this->PadIndex( pos )  >= 0 )
    return gap_base;
  int index = pos;
  for (int ii=0; ii<(int)pads_.size( ); ii++) {
    if ( pads_[ii] > pos )
      break;
    index--;
  }
  return as_base( (*fastb)[ index ] );
}



/*
 * padded_seq
 * Qual
 *
 * The same as Base, only it returns the quality score (as an int). If qualb
 * is empty, then it returns the default finished base quality score.
 */
int padded_seq::Qual( const qualvector *qualb, int pos ) const
{
  if ( pos < 0 || pos >= (int)qualb->size( ) + this->PadsCount( ) )
    return empty_qual;
  if ( this->PadIndex( pos )  >= 0 )
    return gap_qual;
  if ( !qualb )
    return fin_qual;
  int index = pos;
  for (int ii=0; ii<(int)pads_.size( ); ii++) {
    if ( pads_[ii] > pos )
      break;
    index--;
  }
  return int( (*qualb)[index] );
}



/*
 * padded_seq
 * ToVec
 *
 * qualb may be null, in which case quality scores will be defaulted to
 * "finished base" quality score.
 */
void padded_seq::ToVec( const basevector *fastb,
			const qualvector *qualb,
			vec<char> &bases,
			vec<int> &quals,
			int from,
			int to ) const
{
  bases.resize( 0 );
  quals.resize( 0 );

  // Reserve memory.
  int seq_end = fastb->size( ) + this->PadsCount( );
  from = ( from < 0 ) ? 0 : from;
  to = ( to < 0 ) ? seq_end : to;
  
  bases.reserve( to - from );
  quals.reserve( to - from );

  // Fill vectors.
  int index = 0;
  int pos = ( from <= begin_ ) ? 0 : from - begin_;
  for (int ii=from; ii<to; ii++) {
    int loc_ii = ii - begin_;
    if ( loc_ii < 0 ) {
      bases.push_back( empty_base );
      quals.push_back( empty_qual );
    }
    else if ( loc_ii < al_begin_ ) {
      if ( pos >= (int)fastb->size( ) ) {
	// We should not get in here (temporary hack to deal
	//  with an upstream bug).
	bases.push_back( 'x' );
	quals.push_back( 0 );
	pos++;
      }
      else { 
	// Switch to lower case.
	bases.push_back( char( 32 + as_base( (*fastb)[ pos ] ) ) );
	qualb ? quals.push_back((*qualb)[pos]) : quals.push_back(fin_qual);
	pos++;
      }
    }
    else if ( loc_ii > al_begin_ + this->PaddedLength( ) - 1 ) {
      if ( loc_ii > (int)fastb->size( ) + this->PadsCount( ) - 1 )
	break;
      else {
	if ( pos >= (int)fastb->size( ) ) {
	  // We should not get in here (temporary hack to deal
	  //  with an upstream bug).
	  bases.push_back( 'x' );
	  quals.push_back( 0 );
	  pos++;
	}
	else { 
	  // Switch to lower case.
	  bases.push_back( char( 32 + as_base( (*fastb)[ pos ] ) ) );
	  qualb ? quals.push_back((*qualb)[pos]) : quals.push_back(fin_qual);
	  pos++;
	}
      }
    }
    else {
      if ( index < (int)pads_.size( ) ) {
	if ( loc_ii < pads_[index] ) {
	  if ( pos >= (int)fastb->size( ) ) {
	    // We should not get in here (temporary hack to deal
	    //  with an upstream bug).
	    bases.push_back( 'x' );
	    quals.push_back( 0 );
	    pos++;
	  }
	  else { 
	    bases.push_back( as_base( (*fastb)[ pos ] ) );
	    qualb ? quals.push_back((*qualb)[pos]) : quals.push_back(fin_qual);
	    pos++;
	  }
	}
	else {
	  bases.push_back( gap_base );
	  quals.push_back( gap_qual );
	  index++;
	}
      }
      else {
	if ( pos >= (int)fastb->size( ) ) {
	  // We should not get in here (temporary hack to deal
	  //  with an upstream bug).
	  bases.push_back( 'x' );
	  quals.push_back( 0 );
	  pos++;
	}
	else {
	  bases.push_back( as_base( (*fastb)[ pos ] ) );
	  qualb ? quals.push_back( (*qualb)[pos] ) : quals.push_back(fin_qual);
	  pos++;
	}
      }
    }
  }
  
}



/*
 * padded_seq
 * Unpack
 *
 * This is similar to ToVec. The difference is that bases and quals (the
 * output) are the full thing, i.e. all the bases in fastb are placed in
 * bases, and all the scores in qualb are placed in quals.
 * 
 * If filter_ends = True, the fastb will be checked and tagged for low
 * quality end bases (printed in lower case).
 */
void padded_seq::Unpack( const basevector *fastb,
			 const qualvector *qualb,
			 vec<char> &bases,
			 vec<int> &quals,
			 bool filter_ends ) const
{
  bases.clear( );
  quals.clear( );

  if ( qualb )
    ForceAssert( (int)fastb->size( ) == (int)qualb->size( ) );
  
  int full_length = this->PadsCount( ) + (int)fastb->size( ) ;
  bases.resize( full_length, empty_base );
  quals.resize( full_length, empty_qual );

  for (int ii=0; ii<(int)pads_.size( ); ii++) {
    bases[ pads_[ii] ] = gap_base;
    quals[ pads_[ii] ] = gap_qual;
  }
  
  int bpos = 0;
  for (int ii=0; ii<(int)fastb->size( ); ii++) {
    while ( bases[bpos] == gap_base )
      bpos++;
    bases[bpos] = as_base( (*fastb)[ii] );
    quals[bpos] = qualb ? (int)(*qualb)[ii] : fin_qual;
    bpos++;
  }

  // Detect and lower-case low quality end bases.
  if ( ! filter_ends )
    return;

  // No quality scores.
  if ( !qualb )
    return;
  
  // Good portion starts where there are good_len bases of qual>=min_good.
  const int min_good = 20;
  const int good_len = 12;
  
  // Detect boundaries of "good portion" of read.
  int wbeg = 0;
  while ( wbeg < (int)qualb->size( ) - good_len ) {
    bool is_good = true;
    for (int ii=wbeg; ii<wbeg+good_len; ii++) {
      if ( (*qualb)[ii] < min_good ) {
	wbeg = ii+1;
	is_good = false;
	break;
      }
    }
    if ( is_good )
      break;
  }

  int wend = (int)qualb->size( )-1;
  while ( wend >= good_len ) {
    bool is_good = true;
    for (int ii=wend; ii>=wend-good_len; ii--) {
      if ( (*qualb)[ii] < min_good ) {
	wend = ii-1;
	is_good = false;
	break;
      }
    }
    if ( is_good ) {
      wend++;
      break;
    }
  }

  // Replace upper case CGTA with lower case cgta outside good interval.
  bpos = 0;
  for (int ii=0; ii<(int)fastb->size( ); ii++) {
    while ( bases[bpos] == gap_base ) {
      if ( ii < wbeg || ii >= wend )
	ToLowerCase( bases[bpos] );
      bpos++;
    }
    if ( ii < wbeg || ii >= wend )
      ToLowerCase( bases[bpos] );
    bpos++;
  }
  
}



/*
 * padded_seq
 * PrintBrief
 *
 * Send to out a brief description of the padded sequence (if new_line is set
 * to true, terminate with a new line).
 */
void padded_seq::PrintBrief( ostream &out, const bool &new_line ) const
{
  out << "begin: " << begin_ << ", "
      << "al_begin: " << al_begin_ << ", "
      << "length: " << len_ << ", "
      << "#pads: " << pads_.size( );

  if ( pads_.size( ) > 0 ) {
    out << " (";
    if ( (int)pads_.size( ) < 6 ) {
      for (int ii=0; ii<(int)pads_.size( ); ii++) {
	out << pads_[ii];
	if ( ii < (int)pads_.size( ) - 1 )
	  out << " ";
      }
    }
    else {
      for (int ii=0; ii<4; ii++)
	out << pads_[ii] << " ";
      out << "... " << pads_.back( );
    }
    out << ")";
  }
  
  if ( new_line )
    out << "\n";
}



/*
 * padded_seq
 * istream>>
 */
istream &operator>> ( istream &in, padded_seq &pad )
{
  int pos1;
  int pos2;
  int n_blocks;

  in >> pos1
     >> pos2
     >> n_blocks;

  if ( !in )
    return in;

  avector<int> gaps( n_blocks );
  avector<int> lens( n_blocks );

  for (int ii=0; ii<n_blocks; ii++)
    in >> gaps(ii)
       >> lens(ii);

  static align al;
  al.Set( pos1, pos2, gaps, lens );

  pad.FromAlign( al );
  
  return in;
}



/*
 * padded_seq
 * ostream<<
 */
ostream &operator<< ( ostream &out, const padded_seq pad )
{
  static align al;
  pad.ToAlign( al );
  
  out << al.pos1( ) << "  "
      << al.pos2( ) << "  "
      << al.Nblocks( ) << "  ";

  for (int ii=0; ii<(int)al.Nblocks( ); ii++)
    out << al.Gaps(ii) << "  "
	<< al.Lengths(ii) << "  ";
  
  return out;
}



/*
 * padded_seq
 * TrimEndingPads
 *
 * Remove ending pads. In other wors, first and last objects of a padded_seq
 * cannot be pads.
 */
void padded_seq::TrimEndingPads( )
{
  while ( pads_.size( ) > 0 &&
	  pads_[pads_.size( ) - 1] == al_begin_ + this->PaddedLength( ) - 1 )
    pads_.resize( pads_.size( ) - 1 );
  while ( pads_.size( ) > 0 && pads_[0] == 0 )
    pads_.erase( pads_.begin( ) );
}



/*
 * padded_seq
 * PadIndex
 *
 * Returns -1 if pos is not a pad, or else the index in pads_ of pos.
 * Warning: pos is in local coordinates.
 */
int padded_seq::PadIndex( int pos ) const
{
  vec<int>::const_iterator iter;
  iter = lower_bound( pads_.begin( ), pads_.end( ), pos );
  if ( iter == pads_.end( ) || *iter != pos )
    return -1;
  return ( iter - pads_.begin( ) );
}

