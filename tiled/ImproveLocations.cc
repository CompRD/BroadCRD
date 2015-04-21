// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include <map>

#include "math/Functions.h"
#include "STLExtensions.h"
#include "VecAlignmentPlus.h"
#include "tiled/BaggedReads.h"
#include "tiled/ImproveLocations.h"



/*
 * improve_locations
 * Constructor
 */
improve_locations::improve_locations( ) :
  aligns_ ( 0 ),
  locs_ ( 0 ),
  log_ ( 0 ),
  exclude_log_ ( 0 ),
  devnull_ ( "/dev/null" )
{ }



/*
 * improve_locations
 * Constructor
 */
improve_locations::improve_locations( const vec_alignment_plus *aligns,
				      const vec<read_location> *locs,
				      ostream *log,
				      ostream *exclude_log ) :
  aligns_ ( aligns ),
  locs_ ( locs ),
  log_ ( log ),
  exclude_log_ ( exclude_log ),
  devnull_ ( "/dev/null" )
{
  this->Setup( );
}



/*
 * improve_locations
 * SetPointers
 */
void improve_locations::SetPointers( const vec_alignment_plus *aligns,
				     const vec<read_location> *locs,
				     ostream *log,
				     ostream *exclude_log )
{
  aligns_ = aligns;
  locs_ = locs;
  log_ = log;
  exclude_log_ = exclude_log;

  this->Setup( );
}



/*
 * improve_locations
 * Improve
 *
 * Improve locations for all contigs.
 */
void improve_locations::Improve( )
{
  for (int ii=0; ii<(int)first_locs_.size( ); ii++)
    this->Improve( ii );
}



/*
 * improve_locations
 * Improve
 *
 * Improve only locations for the given contig.
 */
void improve_locations::Improve( int contig_id )
{
  ostream &log = ( log_ ) ? *log_ : devnull_;

  int n_locs = locs_->size( );
  int n_contigs = first_locs_.size( );
  int begin = first_locs_[contig_id];
  int end = ( contig_id+1 < n_contigs ) ? first_locs_[contig_id+1] : n_locs;
  
  log << "\ncontig_" << contig_id
      << " (" << end - begin
      << " reads in contig)\n";
  
  // One read contig.
  if ( end - begin < 2 ) {
    for (int ii=begin; ii<end; ii++)
      ilocs_.push_back( (*locs_)[ii] );
    return;
  }
       
  // Minilocs: a small vector of read_locations for contig_id.
  minilocs_.clear( );
  minilocs_.reserve( end - begin );
  for (int ii=begin; ii<end; ii++)
    minilocs_.push_back( (*locs_)[ii] );
  
  // Eliminate trouble locs.
  this->EliminateBadLocs( );

  // Locate chimp only regions and shift locs accordingly.
  this->AdjustLocs( );
  
  // Append contig to ilocs_.
  copy( minilocs_.begin(), minilocs_.end(), back_inserter( ilocs_ ) );

}



/*
 * improve_locations
 * SortAndFill
 */
void improve_locations::SortAndFill( vec<read_location> &improved_locs )
{
  sort( ilocs_.begin( ), ilocs_.end( ) );
  improved_locs = ilocs_;
}



/*
 * improve_locations
 * SortAndSave
 */
void improve_locations::SortAndSave( const String &out_file )
{
  Sort(ilocs_);
  WriteLocs( out_file, ilocs_ );
}



/*
 * improve_locations
 * Setup
 */
void improve_locations::Setup( )
{
  // locs_ must be sorted.
  ForceAssert( is_sorted( locs_->begin( ), locs_->end( ) ) );

  // Reserve memory for improved locs.
  ilocs_.reserve( locs_->size( ) );
  ilocs_.clear( );

  // Generate first locs for original locations.
  first_locs_.clear( );
  first_locs_.resize( 1 + locs_->back( ).Contig( ) );
  for (int ii=(int)locs_->size( ) - 1; ii>=0; ii--)
    first_locs_[ (*locs_)[ii].Contig( ) ] = ii;
}



/*
 * improve_locations
 * EliminateBadLocs
 *
 * For each read r in minilocs_ look at all adjacent reads r' in minilocs_
 * such that r and r' loc-overlap. If none of the pairs as above has
 * an align-overlap, then discard the read.
 */
void improve_locations::EliminateBadLocs( )
{
  ostream &ex_log = ( exclude_log_ ) ? *exclude_log_ : devnull_;
  ostream &log = ( log_ ) ? *log_ : devnull_;
  
  for (int ii=0; ii<(int)minilocs_.size( ); ii++) {
    vec<int> aligners;
    this->OverlapByLoc( ii, aligners );
    
    // There are no aligners.
    if ( aligners.size( ) < 1 ) {
      log << " read " << minilocs_[ii].ReadId( )
	  << " has no loc_aligners. Removing it\n";
      ex_log << "-1\t" << minilocs_[ii].ReadId( )
	     << "\tno_loc_aligners\n";
      minilocs_.erase( minilocs_.begin( ) + ii );
      ii--;
      continue;
    }
    
    // None of the aligners does actually have an alignment.
    int n_found = 0;
    vec<int> align_overlap( (int)aligners.size( ), 0 );
    for (int jj=0; jj<(int)aligners.size( ); jj++) {
      align_overlap[jj] = this->AlignOverlap( ii, aligners[jj] );
      if ( align_overlap[jj] > 0 )
	n_found++;
    }
    if ( 0 == n_found  ) {
      log << " read " << minilocs_[ii].ReadId( )
	  << " has " << aligners.size( ) 
	  << " loc_aligners, but no aligns where found. Removing it\n";
      ex_log << "-1\t" << minilocs_[ii].ReadId( )
	     << "\tno_aligns\n";
      minilocs_.erase( minilocs_.begin( ) + ii );
      ii--;
      continue;
    }
  }
  
}



/*
 * improve_locations
 * AdjustLocs
 */
void improve_locations::AdjustLocs( )
{
  ostream &log = ( log_ ) ? *log_ : devnull_;

  const int radius = 12;
  const int max_discrep = 20;

  // First entry in minilocs_.
  if ( minilocs_.size( ) < 1 )
    return;
  minilocs_[0].SetStartOnContig( 0 );
  
  // Loop over locs at pos>=1 in minilocs_.
  for (int ii=1; ii<(int)minilocs_.size( ); ii++) {
    
    // Starts implied by already placed locs.
    vec< pair<int, int> > mass_start;
    for (int jj=Max( 0, ii-radius ); jj<ii; jj++) {
      int mass = this->AlignOverlap( jj, ii );
      if ( mass > 0 ) {
	int offset = this->AlignOffset( jj, ii );
	int start = minilocs_[jj].StartOnContig( ) + offset;
	mass_start.push_back( make_pair( mass, start ) );
      }
    }
    sort( mass_start.begin( ), mass_start.end( ) );
    
    // Miniloc ii does not align any already placed locs (no starts found).
    if ( mass_start.size( ) < 1 ) {
      log << " " << ii
	  << " (read " << minilocs_[ii].ReadId( )
	  << ") does not align any locs on the left\n";
      int prev_start = minilocs_[ii-1].StartOnContig( );
      int prev_len = minilocs_[ii-1].LengthOfContig( );
      int start = max_discrep + prev_start + prev_len;
      minilocs_[ii].SetStartOnContig( start );
      continue;
    }
    
    // Gauge the window of possible starts.
    vec<int> all_starts;
    for (int jj=0; jj<(int)mass_start.size( ); jj++)
      all_starts.push_back( mass_start[jj].second );

    int low_start = Min( all_starts );
    int high_start = Max( all_starts );
    int delta = high_start - low_start;
    
    // Solid agreement found.
    if ( delta <= max_discrep ) {
      int the_start = (int)Mean( all_starts );
      log << " " << ii
	  << " (read " << minilocs_[ii].ReadId( )
	  << ") placed at " << the_start
	  << " by using " << all_starts.size( )
	   << " align-overlaps\n";
      minilocs_[ii].SetStartOnContig( the_start );
      continue;
    }
    
    // Starts do not completely agree, just take the best mass starts.
    int loc_low = all_starts[0];
    int loc_high = all_starts[0];
    vec<int> acceptable_starts;
    for (int jj=0; jj<(int)all_starts.size( ); jj++) {
      loc_low = Min( loc_low, all_starts[jj] );
      loc_high = Max( loc_high, all_starts[jj] );
      if ( loc_high - loc_low > max_discrep )
	break;
      else
	acceptable_starts.push_back( all_starts[jj] );
    }
    
    int the_start = (int)Mean( acceptable_starts );
    log << " " << ii
	<< " (read " << minilocs_[ii].ReadId( )
	<< ") placed at " << the_start
	<< " by using " << acceptable_starts.size( )
	<< " align-overlaps, out of a total of  " << all_starts.size( )
	<< " found align-overlaps: \n";
    for (int jj=Max( 0, ii-radius ); jj<ii; jj++) {
      int mass = this->AlignOverlap( jj, ii );
      if ( mass > 0 ) {
	int offset = this->AlignOffset( jj, ii );
	int start = minilocs_[jj].StartOnContig( ) + offset;
	log << "   [" << jj
	    << ", " << ii
	    << "]   (" << minilocs_[jj].ReadId( )
	    << ", " << minilocs_[ii].ReadId( )
	    << ")   implied start: " << start
	    << "\n";
      }
    }
    minilocs_[ii].SetStartOnContig( the_start );
  }
  
  // Clean up (first read starts at 0, adjust contig length, and sort).
  ForceAssert( minilocs_.size( ) > 0 );

  int leftmost_base = minilocs_[0].StartOnContig( );
  for (int ii=1; ii<(int)minilocs_.size( ); ii++)
    leftmost_base = Min( leftmost_base, minilocs_[ii].StartOnContig( ) );
  for (int ii=0; ii<(int)minilocs_.size( ); ii++) {
    int old_start = minilocs_[ii].StartOnContig( );
    int new_start = old_start - leftmost_base;
    minilocs_[ii].SetStartOnContig( new_start );
  }

  int rightmost_base = minilocs_[0].StopOnContig( );
  for (int ii=1; ii<(int)minilocs_.size( ); ii++)
    rightmost_base = Max( rightmost_base, minilocs_[ii].StopOnContig( ) );
  int new_contig_length = 1 + rightmost_base;
  for (int ii=0; ii<(int)minilocs_.size( ); ii++)
    minilocs_[ii].SetLengthOfContig( new_contig_length );
  
  sort( minilocs_.begin( ), minilocs_.end( ) );
  
}  



/*
 * improve_locations
 * OverlapByLoc
 *
 * Merge the mini_ids from OverlapLeftByLoc and OverlapRightByLoc.
 */
void improve_locations::OverlapByLoc( int id, vec<int> &mini_ids ) const
{
  mini_ids.clear( );
  vec<int> right_ids;
  
  this->OverlapLeftByLoc( id, mini_ids );
  this->OverlapRightByLoc( id, right_ids );

  copy( right_ids.begin( ), right_ids.end( ), back_inserter( mini_ids ) );
}



/*
 * improve_locations
 * OverlapLeftByLoc
 *
 * Find locations which loc-overlap with miniloc_[id] from the left.
 */
void improve_locations::OverlapLeftByLoc( int id, vec<int> &mini_ids ) const
{
  mini_ids.clear( );
  
  int max_skip = 3;
  int n_skipped = 0;
  for (int ii=id-1; ii>=0; ii--) {
    if ( n_skipped >= max_skip )
      break;
    if ( this->LocOverlap( ii, id ) > 0 )
      mini_ids.push_back( ii );
    else
      n_skipped++;
  }
}



/*
 * improve_locations
 * OverlapRightByLoc
 *
 * Find locations which loc-overlap with miniloc_[id] from the right.
 */
void improve_locations::OverlapRightByLoc( int id, vec<int> &mini_ids ) const
{
  mini_ids.clear( );
  
  int max_skip = 3;
  int n_skipped = 0;
  for (int ii=id+1; ii<(int)minilocs_.size( ); ii++) {
    if ( n_skipped >= max_skip )
      break;
    if ( this->LocOverlap( id, ii ) > 0 )
      mini_ids.push_back( ii );
    else
      n_skipped++;
  }     
}



/*
 * improve_locations
 * GetAlignIndex
 *
 * Returns the index of the alignment between reads at mini1 and mini2
 * in minilocs_.
 */
int improve_locations::GetAlignIndex( int mini1, int mini2 ) const
{
  int read1 = minilocs_[mini1].ReadId( );
  int read2 = minilocs_[mini2].ReadId( );
  int index = aligns_->GetAlignsIndex( read1 );
  if ( index < 0 )
    return -1;
  for (int ii=index; ii<(int)aligns_->GetNumberAlignments( ); ii++) {
    if ( aligns_->GetAlignmentId1( ii ) != read1 )
      break;
    if ( aligns_->GetAlignmentId2( ii ) == read2 )
      return ii;
  }
  return -1;
}



/*
 * improve_locations
 * LocOverlap
 *
 * It returns the amount of overlap betwen minilocs_[mini1] and
 * minilocs_[mini2] based only on locs information.
 */
int improve_locations::LocOverlap( int mini1, int mini2 ) const
{
  const read_location &loc1 = minilocs_[mini1];
  const read_location &loc2 = minilocs_[mini2];

  int left = Max( loc1.StartOnContig( ), loc2.StartOnContig( ) );
  int right = 1 + Min( loc1.StopOnContig( ), loc2.StopOnContig( ) );
  
  return Max( 0 , right - left );
}



/*
 * improve_locations
 * LocOffset
 *
 * It returns the amount of offset betwen minilocs_[mini1] and
 * minilocs_[mini2] based only on locs information.
 */
int improve_locations::LocOffset( int mini1, int mini2 ) const
{
  const read_location &loc1 = minilocs_[mini1];
  const read_location &loc2 = minilocs_[mini2];

  return loc2.StartOnContig( ) - loc1.StartOnContig( );
}



/*
 * improve_locations
 * LocOverlap
 *
 * It returns the amount of overlap betwen minilocs_[mini1] and
 * minilocs_[mini2] based only on the alignment between the reads
 * (if it exists).
 */
int improve_locations::AlignOverlap( int mini1, int mini2 ) const
{
  int align_id = this->GetAlignIndex( mini1, mini2 );

  if ( align_id < 0 )
    return 0;

  alignment_plus ap;
  aligns_->GetAlignment( ap, align_id );
  
  return ap.a.Pos1( ) - ap.a.pos1( );
}



/*
 * improve_locations
 * AlignOffset
 *
 * It returns the amount of align-offset betwen minilocs_[mini1] and
 * minilocs_[mini2]. Warning: if there is no align between the to
 * then it will ForceAssert.
 */
int improve_locations::AlignOffset( int mini1, int mini2 ) const
{
  int align_id = this->GetAlignIndex( mini1, mini2 );

  ForceAssert( align_id >= 0 );
  
  alignment_plus the_ap;
  aligns_->GetAlignment( the_ap, align_id );
  align al = the_ap.a;
  if ( minilocs_[mini1].OrientationOnContig( ) == ReverseOr ) {
    int len1 = minilocs_[mini1].LengthOfRead( );
    int len2 = minilocs_[mini2].LengthOfRead( );
    al.ReverseThis( len1, len2 );
  }

  return al.pos1( ) - al.pos2( );
}



/*
 * improve_locations
 * EvalShift
 *
 * Look at minilocs_ mini1 and mini2. If they do no align_overlap return
 * false, otherwise return the unique shift such that loc_offset + shift
 * = align_offset.
 */
bool improve_locations::EvalShift( int mini1, int mini2, int &shift ) const
{
  shift = 0;

  if ( this->AlignOverlap( mini1, mini2 ) == 0 )
    return false;
  
  int align_offset = this->AlignOffset( mini1, mini2 );
  int loc_offset = this->LocOffset( mini1, mini2 );
  shift = align_offset - loc_offset;

  return true;
}



