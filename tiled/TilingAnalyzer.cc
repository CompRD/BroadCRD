// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "STLExtensions.h"
#include "tiled/CharBaseUtil.h"
#include "tiled/PaddedSeq.h"
#include "tiled/Tiling.h"
#include "tiled/TilingAnalyzer.h"



/*
 * tiling_analyzer
 * Constructor
 */
tiling_analyzer::tiling_analyzer( ) :
  k_bases_ ( 0 ),
  c_bases_ ( 0 ),
  r_bases_ ( 0 ),
  k_quals_ ( 0 ),
  c_quals_ ( 0 ),
  r_quals_ ( 0 ),
  the_tiling_ ( 0 )
{ }



/*
 * tiling_analyzer
 * SetPointers
 *
 * Although both contig and known sequence pointers may be null, they
 * cannot be null at the same time. Notice that reads quals may be
 * zero as well. Finally, if k_quals_ is null, then known sequence
 * quality scores will be set to a given finished-grade quality.
 */
void tiling_analyzer::SetPointers( const vecbasevector *r_bases,
				   const vecqualvector *r_quals,
				   const vecbasevector *c_bases,
				   const vecqualvector *c_quals,
				   const vecbasevector *k_bases,
				   const vecqualvector *k_quals )
{
  r_bases_ = r_bases;
  r_quals_ = r_quals;
  c_bases_ = c_bases;
  c_quals_ = c_quals;
  k_bases_ = k_bases;
  k_quals_ = k_quals;
}



/*
 * tiling_analyzer
 * SetTiling
 *
 * Clean ups the previous tiling data.
 */
void tiling_analyzer::SetTiling( const tiling *the_tiling )
{
  // Consistency check, and set tiling.
  ForceAssert( r_bases_ );
  
  if ( the_tiling->ContigId( ) > -1 ) {
    ForceAssert( c_bases_ );
    if ( r_quals_ )
      ForceAssert( c_quals_ );
  }
  
  if ( the_tiling->KnownId( ) > -1 )
    ForceAssert( k_bases_ );
  
  the_tiling_ = the_tiling;
  
  // Reset ready-to-print bases and quals.
  str_k_bases_.clear( );
  str_k_quals_.clear( );
  str_c_bases_.clear( );
  str_c_quals_.clear( );
  str_r_bases_.clear( );
  str_r_quals_.clear( );
  str_r_bases_.resize( the_tiling_->ReadsCount( ) );
  str_r_quals_.resize( the_tiling_->ReadsCount( ) );

  // Reset windows of alignments on master.
  if ( the_tiling_->KnownId( ) == -1 || the_tiling_->ContigId( ) == -1 )
    c_win_ = make_pair( -1, -1 );
  else {
    const padded_seq &pads = the_tiling_->Contig( ).PaddedSeq( );
    int cg_len = (*c_bases_)[the_tiling_->ContigId( )].size( );
    int beg = pads.Begin( );
    int end = beg + pads.PadsCount( ) + cg_len;
    c_win_ = make_pair( beg, end );
  }

  r_win_.clear( );
  r_win_.reserve( the_tiling_->ReadsCount( ) );
  for (int ii=0; ii<the_tiling_->ReadsCount( ); ii++ ) {
    const padded_seq &pads = the_tiling_->Read( ii ).PaddedSeq( );
    int read_id = the_tiling_->Read( ii ).Id( );
    int read_len = (*r_bases_)[read_id].size( );
    int beg = pads.Begin( );
    int end = beg + pads.PadsCount( ) + read_len;
    r_win_.push_back( make_pair( beg, end ) );
  }

}



/*
 * tiling_analyzer
 * PrintKnown
 *
 * For any ii, it saves in base and qual base and quality score of known
 * sequence.
 */
void tiling_analyzer::PrintKnown( int ii, char &base, int &qual ) const
{
  base = empty_base;
  qual = empty_qual;

  int known_id = the_tiling_->KnownId( );

  if ( !k_bases_ || known_id < 0 )
    return;

  if ( str_k_bases_.size( ) < 1 )
    this->BuildStringsKnown( );

  if ( ii < 0 || ii > (int)str_k_bases_.size( ) - 1 )
    return;

  base = str_k_bases_[ii];
  if ( base != empty_base )
    ToUpperCase( base );
  qual = str_k_quals_[ii];
}



/*
 * tiling_analyzer
 * PrintContig
 *
 * Same, for the contig.
 */
void tiling_analyzer::PrintContig( int ii, char &base, int &qual ) const
{
  base = empty_base;
  qual = empty_qual;

  int contig_id = the_tiling_->ContigId( );

  if ( !c_bases_ || contig_id < 0 )
    return;

  if ( str_c_bases_.size( ) < 1 )
    this->BuildStringsContig( );
  
  ii -= ( the_tiling_->KnownId( ) > -1 ) ? c_win_.first : 0;

  if ( ii < 0 || ii > (int)str_c_bases_.size( ) - 1 )
    return;
  
  base = str_c_bases_[ii];
  if ( base != empty_base )
    ToUpperCase( base );
  qual = str_c_quals_[ii];
}



/*
 * tiling_analyzer
 * PrintRead
 *
 * Same, for the the given read.
 */
void tiling_analyzer::PrintRead( int read, int ii, char &base, int &qual ) const
{
  base = empty_base;
  qual = empty_qual;

  if ( !r_bases_ )
    return;

  if ( str_r_bases_[read].size( ) < 1 )
    this->BuildStringsRead( read );
  
  ii -= r_win_[read].first;
  
  if ( ii < 0 || ii > (int)str_r_bases_[read].size( ) - 1 )
    return;
  
  base = str_r_bases_[read][ii];
  qual = str_r_quals_[read][ii];
}



/*
 * tiling_analyzer
 * QualKnown
 *
 * If the base at ii is a regular (non-pad) base it returns its quality
 * score. If it is a pad then it returns the minimum of the quality
 * scores of the two closest (left and right) non-pad bases. In all
 * the other cases (or on error) it returns empty_qual.
 */
int tiling_analyzer::QualKnown( int ii ) const
{
  int known_id = the_tiling_->KnownId( );

  if ( !k_bases_ || known_id < 0 )
    return empty_qual;

  if ( str_k_bases_.size( ) < 1 )
    this->BuildStringsKnown( );

  if ( ii < 0 || ii > (int)str_k_bases_.size( ) - 1 )
    return empty_qual;
  
  char base = str_k_bases_[ii];
  if ( base != gap_base )
    return str_k_quals_[ii];
  
  int left_qual = empty_qual;
  for (int jj=ii-1; jj>=0; jj--) {
    if ( str_k_bases_[jj] != gap_base ) {
      left_qual = str_k_quals_[jj];
      break;
    }
  }
  int right_qual = empty_qual;
  for (int jj=ii+1; jj<(int)str_k_bases_.size( ); jj++) {
    if ( str_k_bases_[jj] != gap_base ) {
      right_qual = str_k_quals_[jj];
      break;
    }
  }

  return Min( left_qual, right_qual );
}



/*
 * tiling_analyzer
 * QualContig
 *
 * Same, for the contig.
 */
int tiling_analyzer::QualContig( int ii ) const
{
  int contig_id = the_tiling_->ContigId( );

  if ( !c_bases_ || contig_id < 0 )
    return empty_qual;

  if ( str_c_bases_.size( ) < 1 )
    this->BuildStringsContig( );
  
  ii -= ( the_tiling_->KnownId( ) > -1 ) ? c_win_.first : 0;

  if ( ii < 0 || ii > (int)str_c_bases_.size( ) - 1 )
    return empty_qual;
  
  char base = str_c_bases_[ii];
  if ( base != gap_base )
    return str_c_quals_[ii];
  
  int left_qual = empty_qual;
  for (int jj=ii-1; jj>=0; jj--) {
    if ( str_c_bases_[jj] != gap_base ) {
      left_qual = str_c_quals_[jj];
      break;
    }
  }
  int right_qual = empty_qual;
  for (int jj=ii+1; jj<(int)str_c_bases_.size( ); jj++) {
    if ( str_c_bases_[jj] != gap_base ) {
      right_qual = str_c_quals_[jj];
      break;
    }
  }

  return Min( left_qual, right_qual );
}



/*
 * tiling_analyzer
 * PrintRead
 *
 * Same, for the the given read.
 */
int tiling_analyzer::QualRead( int read, int ii ) const
{
  if ( !r_bases_ )
    return empty_qual;

  if ( str_r_bases_[read].size( ) < 1 )
    this->BuildStringsRead( read );
  
  ii -= r_win_[read].first;
  
  if ( ii < 0 || ii > (int)str_r_bases_[read].size( ) - 1 )
    return empty_qual;
  
  char base = str_r_bases_[read][ii];
  if ( base != gap_base )
    return str_r_quals_[read][ii];
  
  int left_qual = empty_qual;
  for (int jj=ii-1; jj>=0; jj--) {
    if ( str_r_bases_[read][jj] != gap_base ) {
      left_qual = str_r_quals_[read][jj];
      break;
    }
  }
  int right_qual = empty_qual;
  for (int jj=ii+1; jj<(int)str_r_bases_[read].size( ); jj++) {
    if ( str_r_bases_[read][jj] != gap_base ) {
      right_qual = str_r_quals_[read][jj];
      break;
    }
  }

  return Min( left_qual, right_qual );
}



/*
 * tiling_analyzer
 * QualToChar
 *
 * It returns the first digit of the given quality score.
 */
char tiling_analyzer::QualToChar( int n_qual ) const
{
  if ( gap_qual == n_qual )
    return gap_base;

  if ( empty_qual == n_qual )
    return empty_base;
  
  if ( n_qual < 10 )
    return '0';
  
  return ToString( n_qual )[0];
}



/*
 * tiling_analyzer
 * StartOnContig
 *
 * It returns the start on contig in the read_location sense (in other
 * words after "removing" pads).
 */
int tiling_analyzer::StartOnContig( int read ) const
{
  return this->OnContig( r_win_[read].first );
}



/*
 * tiling_analyzer
 * StopOnContig
 *
 * It returns the stop on contig in the read_location sense (in other
 * words after "removing" pads). Notice this is the stop, not the end.
 */
int tiling_analyzer::StopOnContig( int read ) const
{
  return this->OnContig( r_win_[read].second - 1 );
}



/*
 * tiling_analyzer
 * WinBegin
 *
 * It returns the leftmost value for all windows of aligned objects on master.
 */
int tiling_analyzer::WinBegin( ) const
{
  int begin = 0;

  if ( the_tiling_->KnownId( ) != -1 && the_tiling_->ContigId( ) != -1 )
    begin = c_win_.first;
  else
    if ( (int)r_win_.size( ) > -1 )
      begin = r_win_[0].first;
  
  for (int ii=0; ii<(int)r_win_.size( ); ii++)
    begin = Min( begin, r_win_[ii].first );

  return begin;
}



/*
 * tiling_analyzer
 * WinEnd
 *
 * It returns the rightmost value for all windows of aligned objects on master.
 */
int tiling_analyzer::WinEnd( ) const
{
  int end = 0;
  
  if ( the_tiling_->KnownId( ) != -1 && the_tiling_->ContigId( ) != -1 )
    end = c_win_.second;
  else
    if ( (int)r_win_.size( ) > -1 )
      end = r_win_[0].second;
  
  for (int ii=0; ii<(int)r_win_.size( ); ii++)
    end = Max( end, r_win_[ii].second );

  return end;
}



/*
 * tiling_analyzer
 * CWin
 *
 * Window on master of alignment of contig. Notice that a default pair
 * for "empty" contig (when no contig alignment is given) is built in
 * SetTiling.
 */
pair<int, int> tiling_analyzer::CWin( ) const
{
  return c_win_;
}



/*
 * tiling_analyzer
 * RWin
 *
 * Window on master of alignment of pos-th read
 */
pair<int, int> tiling_analyzer::RWin( int pos ) const
{
  return r_win_[pos];
}



/*
 * tiling_analyzer
 * Tiling
 */
const tiling &tiling_analyzer::Tiling( ) const
{
  return *the_tiling_;
}



/*
 * tiling_analyzer
 * PaddedLengthKnown
 */
int tiling_analyzer::PaddedLengthKnown( ) const
{
  int known_id = the_tiling_->KnownId( );
  if ( known_id < 0 || !k_bases_ )
    return 0;
  const padded_seq &pads = the_tiling_->Known( ).PaddedSeq( );
  return (*k_bases_)[known_id].size( ) + pads.PadsCount( );
}



/*
 * tiling_analyzer
 * PaddedLengthContig
 */
int tiling_analyzer::PaddedLengthContig( ) const
{
  int contig_id = the_tiling_->ContigId( );
  if ( contig_id < 0 || !c_bases_ )
    return 0;
  const padded_seq &pads = the_tiling_->Contig( ).PaddedSeq( );
  return (*c_bases_)[contig_id].size( ) + pads.PadsCount( );
}



/*
 * tiling_analyzer
 * PaddedLengthRead
 */
int tiling_analyzer::PaddedLengthRead( int pos ) const
{
  if ( !r_bases_ )
    return 0;
  int read_id = the_tiling_->Read( pos ).Id( );
  const padded_seq &pads = the_tiling_->Read( pos ).PaddedSeq( );
  return (*r_bases_)[read_id].size( ) + pads.PadsCount( );
}



/*
 * tiling_analyzer
 * TamedRectangles
 *
 * The difference with Rectangles down below is that in here read_pos
 * is given as an input, not as an output. You can specify which reads
 * (and/or known, contig) to print. For example, if read_pos={-2,-1,7},
 * then known (for the -2), contig (for the -1), and reads_[7] will
 * be printed in rbases and rquals. Notice read_pos must be sorted.
 * 
 *  read_pos[ii] = -2 means that ii is known sequence;
 *  read_pos[ii] = -1 means that ii is contig.
 */
void tiling_analyzer::TamedRectangles( int post,
				       int width,
				       const vec<int> &read_pos,
				       vec< vec<char> > &rbases,
				       vec< vec<int> > &rquals ) const
{
  // Clean up.
  rbases.clear( );
  rquals.clear( );

  // Check sortedness.
  ForceAssert( is_sorted( read_pos.begin( ), read_pos.end( ) ) );
  
  // Known.
  vec<char> row_bases;
  vec<int> row_quals;
  
  if ( binary_search( read_pos.begin( ), read_pos.end( ), -2 ) ) {
    row_bases.reserve( width );
    row_quals.reserve( width );
    
    if ( the_tiling_->KnownId( ) > -1 ) {
      for (int ii=post; ii<post+width; ii++) {
	char the_b;
	int the_q;
	this->PrintKnown( ii, the_b, the_q );
	row_bases.push_back( the_b );
	row_quals.push_back( the_q );
      }
      
      rbases.push_back( row_bases );
      rquals.push_back( row_quals );
    }
  }
  
  // Contig, when contig is master sequence (i.e. there is no known).
  row_bases.clear( );
  row_quals.clear( );
  
  if ( binary_search( read_pos.begin( ), read_pos.end( ), -1 ) ) {
    if ( the_tiling_->ContigId( ) > -1 && the_tiling_->KnownId( ) == -1 ) {
      for (int ii=post; ii<post+width; ii++) {
	char the_b;
	int the_q;
	this->PrintContig( ii, the_b, the_q );
	row_bases.push_back( the_b );
	row_quals.push_back( the_q );
      }
      
      rbases.push_back( row_bases );
      rquals.push_back( row_quals );
    }
  }
  
  // Contig, when contig is not master sequence (i.e. known is present).
  row_bases.clear( );
  row_quals.clear( );

  if ( binary_search( read_pos.begin( ), read_pos.end( ), -1 ) ) {
    if ( the_tiling_->ContigId( ) > -1 && the_tiling_->KnownId( ) > -1 ) {
      pair<int, int> c_win = this->CWin( );
      
      int aa = c_win.first;
      int bb = c_win.second;
      int cc = post;
      int dd = post + width;
      if ( IntervalOverlap( aa, bb, cc, dd ) > 0 ) {
	for (int ii=post; ii<post+width; ii++) {
	  char the_b;
	  int the_q;
	  this->PrintContig( ii, the_b, the_q );
	  row_bases.push_back( the_b );
	  row_quals.push_back( the_q );
	}
	
	rbases.push_back( row_bases );
	rquals.push_back( row_quals );
      }
    }
  }

  // Reads.
  for (int pos_id=0; pos_id<(int)read_pos.size( ); pos_id++) {
    int read_id = read_pos[pos_id];
    if ( read_id < 0 || read_id >= the_tiling_->ReadsCount( ) )
      continue;

    row_bases.clear( );
    row_quals.clear( );
    
    pair<int, int> r_win = this->RWin( read_id );
    
    int aa = r_win.first;
    int bb = r_win.second;
    int cc = post;
    int dd = post + width;
    if ( IntervalOverlap( aa, bb, cc, dd ) ) {
      for (int ii=post; ii<post+width; ii++) {
	char the_b;
	int the_q;
	this->PrintRead( read_id, ii, the_b, the_q );
	row_bases.push_back( the_b );
	row_quals.push_back( the_q );
      }
      
      rbases.push_back( row_bases );
      rquals.push_back( row_quals );
    }
  }

}



/*
 * tiling_analyzer
 * Rectangles
 *
 * Generates rectangles of bases and quals starting at post with given width.
 * Will also fill read_pos with the positions of the reads in the tiling
 * sets (so rbases[ii] corresponds to reads_[pos[ii]] in the tiling). Notice
 * that known and contig (if exist) are both assigned read_pos < 0:
 *
 *  read_pos[ii] = -2 means that ii is known sequence;
 *  read_pos[ii] = -1 means that ii is contig.
 */
void tiling_analyzer::Rectangles( int post,
				  int width,
				  vec<int> &read_pos,
				  vec< vec<char> > &rbases,
				  vec< vec<int> > &rquals ) const
{
  // Clean up.
  read_pos.clear( );
  rbases.clear( );
  rquals.clear( );

  // Known.
  vec<char> row_bases;
  vec<int> row_quals;
  row_bases.reserve( width );
  row_quals.reserve( width );
  
  if ( the_tiling_->KnownId( ) > -1 ) {
    for (int ii=post; ii<post+width; ii++) {
      char the_b;
      int the_q;
      this->PrintKnown( ii, the_b, the_q );
      row_bases.push_back( the_b );
      row_quals.push_back( the_q );
    }
    
    rbases.push_back( row_bases );
    rquals.push_back( row_quals );
    read_pos.push_back( -2 );
  }
  
  // Contig, when contig is master sequence (i.e. there is no known).
  row_bases.clear( );
  row_quals.clear( );
  
  if ( the_tiling_->ContigId( ) > -1 && the_tiling_->KnownId( ) == -1 ) {
    for (int ii=post; ii<post+width; ii++) {
      char the_b;
      int the_q;
      this->PrintContig( ii, the_b, the_q );
      row_bases.push_back( the_b );
      row_quals.push_back( the_q );
    }
    
    rbases.push_back( row_bases );
    rquals.push_back( row_quals );
    read_pos.push_back( -1 );
  }
  
  // Contig, when contig is not master sequence (i.e. known is present).
  row_bases.clear( );
  row_quals.clear( );
  
  if ( the_tiling_->ContigId( ) > -1 && the_tiling_->KnownId( ) > -1 ) {
    pair<int, int> c_win = this->CWin( );
    
    int aa = c_win.first;
    int bb = c_win.second;
    int cc = post;
    int dd = post + width;
    if ( IntervalOverlap( aa, bb, cc, dd ) > 0 ) {
      for (int ii=post; ii<post+width; ii++) {
	char the_b;
	int the_q;
	this->PrintContig( ii, the_b, the_q );
	row_bases.push_back( the_b );
	row_quals.push_back( the_q );
      }
      
      rbases.push_back( row_bases );
      rquals.push_back( row_quals );
      read_pos.push_back( -1 );
    }
  }

  // Reads.
  for (int read_id=0; read_id<the_tiling_->ReadsCount( ); read_id++) {
    row_bases.clear( );
    row_quals.clear( );
    
    pair<int, int> r_win = this->RWin( read_id );
    
    int aa = r_win.first;
    int bb = r_win.second;
    int cc = post;
    int dd = post + width;
    if ( IntervalOverlap( aa, bb, cc, dd ) ) {
      for (int ii=post; ii<post+width; ii++) {
	char the_b;
	int the_q;
	this->PrintRead( read_id, ii, the_b, the_q );
	row_bases.push_back( the_b );
	row_quals.push_back( the_q );
      }
      
      rbases.push_back( row_bases );
      rquals.push_back( row_quals );
      read_pos.push_back( read_id );
    }
  }

}



/**
 * tiling_analyzer
 * ContigToReadPos
 *
 * Given pos (position) on the contig, find all the reads that cover
 * it, and return pairs of read_id, pos_on_read for all of them.
 *
 * WARNING: both contig pos and pos on reads are NOT in padded
 * coordinates, but in their natural basevector coordinates.
 *
 * WARNING: if a read is rc, then you first need to rc the basevector,
 * and then fish the base at the given pos on the rc-ed vector.
 */
void tiling_analyzer::ContigToReadPos( int pos, vec<CPos> &rpos ) const
{
  rpos.clear( );

  // Cannot run if there is no contig.
  int contig_id = the_tiling_->ContigId( );
  ForceAssert( contig_id > -1 );
  
  // Translate pos on contig to padded coordinate.
  int padpos = the_tiling_->Contig( ).PaddedSeq( ).ToPaddedPos( pos );
  
  // Loop over all reads covering the given spot on contig.
  for (int read_id=0; read_id<the_tiling_->ReadsCount( ); read_id++) {
    pair<int,int> rwin = this->RWin( read_id );
    if ( IntervalOverlap( rwin.first, rwin.second, padpos, padpos+1 ) ) {
      const tile &rtile = the_tiling_->Read( read_id );
      const padded_seq &padseq = rtile.PaddedSeq( );
      int the_id = rtile.Id( );
      int the_padpos = padpos - r_win_[read_id].first;
      bool is_pad = padseq.IsPad( the_padpos );
      int the_pos = is_pad ? -1 : padseq.ToUnpaddedPos( the_padpos );
      Bool rc = rtile.RC( );
      
      CPos newrpos( the_id, the_pos, rc );
      rpos.push_back( newrpos );
    }
  }
  
}



/*
 * tiling_analyzer
 * BuildStringsKnown
 */
void tiling_analyzer::BuildStringsKnown( ) const
{
  int known_id = the_tiling_->KnownId( );
  if ( known_id < 0 )
    return;
  const padded_seq &pads = the_tiling_->Known( ).PaddedSeq( );
  const basevector *bases = &((*k_bases_)[known_id]);
  const qualvector *quals = k_quals_ ? &((*k_quals_)[known_id]) : 0;
  
  pads.Unpack( bases, quals, str_k_bases_, str_k_quals_ );
}



/*
 * tiling_analyzer
 * BuildStringsContig
 */
void tiling_analyzer::BuildStringsContig( ) const
{
  bool RC = the_tiling_->Contig( ).RC( );
  int contig_id = the_tiling_->ContigId( );
  if ( contig_id < 0 )
    return;
  const padded_seq &pads = the_tiling_->Contig( ).PaddedSeq( );
  
  basevector bases_rc;
  qualvector quals_rc;
  if ( RC ) {
    bases_rc = (*c_bases_)[contig_id];
    bases_rc.ReverseComplement( );
    if ( c_quals_ ) {
      quals_rc = (*c_quals_)[contig_id];
      quals_rc.ReverseMe( );
    }
  }
  const basevector *p_bases_rc = &bases_rc;
  const basevector *bases = RC ? p_bases_rc : &((*c_bases_)[contig_id]);

  const qualvector *quals = 0;
  if ( c_quals_ ) {
    const qualvector *p_quals_rc = &quals_rc;
    quals = RC ? p_quals_rc : &((*c_quals_)[contig_id]);
  }
  
  pads.Unpack( bases, quals, str_c_bases_, str_c_quals_ );
}



/*
 * tiling_analyzer
 * BuildStringsRead
 */
void tiling_analyzer::BuildStringsRead( int pos ) const
{
  bool RC = the_tiling_->Read( pos ).RC( );
  int read_id = the_tiling_->Read( pos ).Id( );
  const padded_seq &pads = the_tiling_->Read( pos ).PaddedSeq( );
  
  basevector bases_rc;
  qualvector quals_rc;
  if ( RC ) {
    bases_rc = (*r_bases_)[read_id];
    bases_rc.ReverseComplement( );
    if ( r_quals_ ) {
      quals_rc = (*r_quals_)[read_id];
      quals_rc.ReverseMe( );
    }
  }
  const basevector *p_bases_rc = &bases_rc;
  const basevector *bases = RC ? p_bases_rc : &((*r_bases_)[read_id]);

  const qualvector *quals = 0;
  if ( r_quals_ ) {
    const qualvector *p_quals_rc = &quals_rc;
    quals = RC ? p_quals_rc : &((*r_quals_)[read_id]);
  }

  pads.Unpack( bases, quals, str_r_bases_[pos], str_r_quals_[pos] );
}



/*
 * tiling_analyzer
 * OnContig
 *
 * Translates padded into unpadded coordinate on the contig.
 */
int tiling_analyzer::OnContig( int coordinate ) const
{
  if ( str_c_bases_.size( ) < 1 )
    this->BuildStringsContig( );
  
  bool have_k = ( the_tiling_->Known( ).Id( ) != -1 );
  int result = coordinate;
  int c_end = ( have_k ) ? c_win_.second : (int)str_c_bases_.size( );
  int beg_loop = ( have_k ) ? c_win_.first : 0;
  int end_loop = Min( c_end, coordinate );

  for (int ii=beg_loop; ii<end_loop; ii++) {
    char base = ' ';
    int qual = 0;
    this->PrintContig( ii, base, qual );
    if ( base == gap_base )
      result--;
  }
  
  return result;
}



