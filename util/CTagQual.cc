/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoverageAnalyzer.h"
#include "Qualvector.h"
#include "SeqInterval.h"
#include "SeqIntervalUtil.h"
#include "math/Functions.h"
#include "util/CTagQual.h"

/**
 * CTagQual
 * Constructor
 */
CTagQual::CTagQual( const vecqvec *cg_quals ) :
  cg_quals_ ( cg_quals ),
  r_quals_ ( 0 ),
  locs_ ( 0 )
{
  this->SetDefaults( );
}

/**
 * CTagQual
 * SetCgQuals
 */
void CTagQual::SetCgQuals( const vecqvec *cg_quals )
{
  cg_quals_ = cg_quals;
}

/**
 * CTagQual
 * SetRPointers
 */
void CTagQual::SetRPointers( const vecqvec *r_quals, const lhandler *locs )
{
  r_quals_ = r_quals;
  locs_ = locs;
}

/**
 * CTagQual
 * SetPit
 */
void CTagQual::SetPit( int min_cov, int min_qual, int min_size )
{
  pit_min_cov_ = min_cov;
  pit_min_qual_ = min_qual;
  pit_min_size_ = min_size;
}

/**
 * CTagQual
 * SelectOverlapping
 */
void CTagQual::SelectOverlapping( const vec<seq_interval> &annot )
{
  vec<seq_interval> select_tags;
  SelectIntersecting( cg_quals_->size(), tags_, annot, select_tags );
  swap( tags_, select_tags );
  this->SortAndRenum( );
}

/**
 * CTagQual
 * LowQual
 *
 * It will merge tag windows closer than fudge_.
 */
const vec<seq_interval> &CTagQual::LowQual( bool append )
{
  if ( ! append ) this->CleanUp( );

  // Tag windows.
  for (vecqvec::size_type contig_id=0; contig_id<cg_quals_->size( ); contig_id++) {
    const qvec &qual = (*cg_quals_)[contig_id];
    if ( qual.size( ) < 1 ) continue;

    int pos = 0;
    int begin = 0;
    bool inside_bad = (int)qual[0] < min_good_;

    while ( pos < (int)qual.size( ) ) {
      if ( inside_bad ) {
	while ( pos < (int)qual.size( ) && qual[pos] < min_good_ ) pos++;
	// TODO: Potentially dangerous truncation of quals size
	seq_interval interv( (int)tags_.size( ), contig_id, begin, pos );
	tags_.push_back( interv );
	inside_bad = false;
      }
      else {
	while ( pos < (int)qual.size( ) && qual[pos] >= min_good_ ) pos++;
	inside_bad = true;
	begin = pos;
      }
    }
  }

  // Merge close windows.
  this->FudgeMerge( );

  // Finalize windows.
  this->SortAndRenum( );

  // Return.
  return tags_;
}

/**
 * CTagQual
 * LowPitCoverage
 *
 * It will merge tag windows closer than fudge_. Only the good portions from a
 * read will be used to compute coverage, where the "good portion" starts where
 * there are pit_min_size bases of qual >= pit_min_qual_.
 */
const vec<seq_interval> &CTagQual::LowPitCoverage( bool append )
{
  if ( ! append ) this->CleanUp( );

  // Lengths of contigs.
  vec<int> cglens;
  cglens.reserve( cg_quals_->size( ) );
  for (vecqvec::size_type ii=0; ii<cg_quals_->size( ); ii++)
    cglens.push_back( (*cg_quals_)[ii].size( ) );

  // Intervals of reads on contigs (pitted, only good chunks are used).
  vec<seq_interval> intervals;
  intervals.reserve( r_quals_->size( ) );

  // Loop over locs.
  for (int loc_id=0; loc_id<locs_->Size( ); loc_id++) {
    const read_location &loc = (*locs_)[loc_id];
    int contig_id = loc.Contig( );
    int read_id = loc.ReadId( );
    const qvec *qualb = &(*r_quals_)[read_id];

    // Detect boundaries of "good portion" of read.
    int wbeg = 0;
    while ( wbeg < (int)qualb->size( ) - pit_min_size_ ) {
      bool is_good = true;
      for (int ii=wbeg; ii<wbeg+pit_min_size_; ii++) {
	if ( (*qualb)[ii] < pit_min_qual_ ) {
	  wbeg = ii+1;
	  is_good = false;
	  break;
	}
      }
      if ( is_good )
	break;
    }

    int wend = (int)qualb->size( )-1;
    while ( wend >= pit_min_size_ ) {
      bool is_good = true;
      for (int ii=wend; ii>=wend-pit_min_size_; ii--) {
	if ( (*qualb)[ii] < pit_min_qual_ ) {
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

    // This should not happen often.
    if ( wend <= wbeg )
      continue;

    // The seq_interval for this read (use loc_id for interval_id).
    int begin = loc.StartOnContig( ) + wbeg;
    int end = begin + ( wend - wbeg );
    seq_interval newseq( loc_id, contig_id, begin, end );
    intervals.push_back( newseq );
  }

  // Find low coverage windows.
  CoverageAnalyzer analyzer( intervals, &cglens );
  analyzer.GetCoveragesAtMost( pit_min_cov_ - 1, tags_ );

  // Merge close windows.
  this->FudgeMerge( );

  // Finalize windows.
  this->SortAndRenum( );

  // Return.
  return tags_;
}

/**
 * CTagQual
 * PrintInfo
 *
 * Warning! ftags_ must exist.
 */
void CTagQual::PrintInfo( ostream &out ) const
{
  out << "Tagging as bad:\n";
  for (size_t contig_id=0; contig_id<cg_quals_->size( ); contig_id++) {
    String cg_line = "c_" + ToString( contig_id );
    int ftag = ftags_[contig_id];
    bool empty = true;
    if ( ftag < 0 ) continue;
    for (int ii=ftag; ii<(int)tags_.size( ); ii++) {
      if ( static_cast<size_t>(tags_[ii].SeqId()) != contig_id ) break;
      String str_beg = ToString( tags_[ii].Begin( ) );
      String str_end = ToString( tags_[ii].End( ) );
      cg_line += " [" + str_beg + "," + str_end + ")";
      empty = false;
    }
    cg_line += " len=" + ToString( (*cg_quals_)[contig_id].size( ) );
    if ( ! empty ) out << cg_line << "\n";
  }
  out << endl;
}

/**
 * CTagQual
 * SetDefaults
 */
void CTagQual::SetDefaults( )
{
  min_good_ = 45;
  fudge_ = 12;
  tag_edges_ = true;

  pit_min_cov_ = 12;
  pit_min_qual_ = 25;
  pit_min_size_ = 12;
}

/**
 * CTagQual
 * CleanUp
 */
void CTagQual::CleanUp( )
{
  tags_.clear( );
  ftags_.clear( );
}

/**
 * CTagQual
 * GenerateFtags
 *
 * Will also renum interval ids (setting tags_[ii].IntervalId( ) = ii).
 */
void CTagQual::GenerateFtags( )
{
  sort( tags_.begin( ), tags_.end( ) );

  ftags_.clear( );
  ftags_.resize( cg_quals_->size( ), -1 );
  for (int ii=tags_.size( )-1; ii>=0; ii--) {
    tags_[ii].SetIntervalId( ii );
    ftags_[ tags_[ii].SeqId( ) ] = ii;
  }
}

/**
 * CTagQual
 * FudgeMerge
 */
void CTagQual::FudgeMerge( )
{
  // Do ftags.
  this->GenerateFtags( );

  // Merge intervals closer together than fudge_.
  for (size_t contig_id=0; contig_id<cg_quals_->size( ); contig_id++) {
    int fbeg = ftags_[contig_id];
    if ( fbeg < 0 )
      continue;
    size_t fend = fbeg;
    while ( fend < tags_.size() && static_cast<size_t>(tags_[fend].SeqId()) == contig_id )
      fend++;
    for (size_t ii=fend-1; ii>static_cast<size_t>(fbeg); ii--) {
      if ( Abs( tags_[ii].Begin( ) - tags_[ii-1].End( ) ) < fudge_ ) {
	tags_[ii-1].SetEnd( tags_[ii].End( ) );
	tags_[ii].SetBegin( -1 );
      }
    }
  }

  // Clean up and regenerate ftags.
  vec<seq_interval> new_tags;
  for (int ii=0; ii<(int)tags_.size( ); ii++)
    if ( tags_[ii].Begin( ) > -1 )
      new_tags.push_back( tags_[ii] );
  swap( tags_, new_tags );

  // Redo ftags.
  this->GenerateFtags( );
}

/**
 * CTagQual
 * SortAndRenum
 */
void CTagQual::SortAndRenum( )
{
  // This will also sort the tags.
  this->GenerateFtags( );

  // No discards, can return.
  if ( tag_edges_ )
    return;

  // Discard windows at ends of contigs.
  vec<seq_interval> new_tags;
  new_tags.reserve( tags_.size( ) );
  for (int ii=0; ii<(int)tags_.size( ); ii++) {
    int contig_id = tags_[ii].SeqId( );
    int contig_len = (int)(*cg_quals_)[contig_id].size( );
    if ( tags_[ii].Begin( ) == 0 ) continue;
    if ( tags_[ii].End( ) == contig_len ) continue;
    new_tags.push_back( tags_[ii] );
  }
  swap( tags_, new_tags );

  // Renum (and sort).
  this->GenerateFtags( );
}

