// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "ScoreAlignment.h"
#include "tiled/AlignsPath.h"



/*
 * aligns_path
 * Constructor
 */
aligns_path::aligns_path( ) :
  aligns_ ( 0 ),
  aligns_index_ ( 0 ),
  hits_ ( 0 ),
  read_bags_ ( 0 ),
  bases_ ( 0 ),
  quals_ ( 0 ),
  log_ ( 0 )
{ }



/*
 * aligns_path
 * Constructor
 */
aligns_path::aligns_path( const vec<alignment_plus> *aligns,
			  const vec<int> *aligns_index,
			  const vec<look_align_plus> *hits,
			  const bagged_reads *read_bags,
			  const vecbasevector *bases,
			  const vecqualvector *quals,
			  ostream *log ) :
  aligns_ ( aligns ),
  aligns_index_ ( aligns_index ),
  hits_ ( hits ),
  read_bags_ ( read_bags ),
  bases_ ( bases ),
  quals_ ( quals ),
  log_ ( log )
{
  this->Setup( );
}



/*
 * aligns_path
 * SetPointers
 */
void aligns_path::SetPointers( const vec<alignment_plus> *aligns,
			       const vec<int> *aligns_index,
			       const vec<look_align_plus> *hits,
			       const bagged_reads *read_bags,
			       const vecbasevector *bases,
			       const vecqualvector *quals,
			       ostream *log )
{
  aligns_ = aligns;
  aligns_index_ = aligns_index;
  hits_ = hits;
  read_bags_ = read_bags;
  bases_ = bases;
  quals_ = quals;
  log_ = log;

  this->Setup( );
}



/*
 * aligns_path
 * SetMaxAlignScore
 */
void aligns_path::SetMaxAlignScore( float max_align_score )
{
  max_align_score_ = max_align_score;
}



/*
 * aligns_path
 * SetMaxPolyScore
 */
void aligns_path::SetMaxPolyScore( int max_poly_score )
{
  max_poly_score_ = max_poly_score;
}



/*
 * aligns_path
 * SetUsePolyScore
 */
void aligns_path::SetUsePolyScore( bool use_poly_score )
{
  use_poly_score_ = use_poly_score;
}



/*
 * aligns_path
 * AreConnected
 *
 * Recursive algorithm to decide if the pos1-th and pos2-th locations in
 * loc_indices can be joined by a walk on good aligns only.
 */
bool aligns_path::AreConnected( int pos1,
				int pos2,
				const vec<int> &loc_indices,
				vec<int> &dead ) const
{
  if ( pos2 < pos1 )
    swap( pos1, pos2 );

  int loc1 = loc_indices[pos1];
  bool is_alive_loc1 = false;

  for (int ii=pos1+1; ii<pos2+1; ii++) {
    int loc2 = loc_indices[ii];
    if ( find( dead.begin( ), dead.end( ), loc2 ) != dead.end( ) )
      continue;
    if ( !this->LocOverlap( loc1, loc2 ) )
      continue;
    if ( !this->AlignOverlap( loc1, loc2 ) )
      continue;
    is_alive_loc1 = true;
    if ( log_ ) 
      *log_ << loc1 << " (id="
	    << (*read_bags_)[loc1].ReadId( ) << ") is well connected with "
	    << loc2 << " (id="
	    << (*read_bags_)[loc2].ReadId( ) << ")\n";
    if ( loc2 == loc_indices[pos2] )
      return true;
    if ( this->AreConnected( ii, pos2, loc_indices, dead ) )
      return true;
  }
  if ( !is_alive_loc1 ) {
    if ( log_ )
      *log_ << loc1 << " (id="
	    << (*read_bags_)[loc1].ReadId( ) << ") is a dead branch\n";
    dead.push_back( loc1 );
  }
  
  return false;
}



/*
 * aligns_path
 * IsBadLoc
 */
bool aligns_path::IsBadLoc( int pos1, const vec<int> &loc_indices ) const
{
  return ( IsDeadForward( pos1, loc_indices ) &&
	   IsDeadBackward( pos1, loc_indices ) );
}



/*
 * aligns_path
 * IsDeadForward
 */
bool aligns_path::IsDeadForward( int pos1,
				 const vec<int> &loc_indices ) const
{
  if ( pos1 > (int)loc_indices.size( ) - 2 )
    return true;

  int max_skip_align = 3;
  int skipped_count = 0;
  int loc1 = loc_indices[pos1];
  for (int ii=pos1+1; ii<(int)loc_indices.size( ); ii++) {
    int loc2 = loc_indices[ii];
    if ( skipped_count > max_skip_align )
      break;
    if ( !this->LocOverlap( loc1, loc2 ) ) {
      skipped_count++;
      continue;
    }
    if ( !this->AlignOverlap( loc1, loc2 ) )
      continue;
    return false;
  }
  
  return true;
}



/*
 * aligns_path
 * IsDeadBackward
 */
bool aligns_path::IsDeadBackward( int pos1,
				  const vec<int> &loc_indices ) const
{
  if ( pos1 < 1 )
    return true;
  
  int max_skip_align = 3;
  int skipped_count = 0;
  int loc1 = loc_indices[pos1];
  for (int ii=pos1-1; ii>-1; ii--) {
    int loc2 = loc_indices[ii];
    if ( skipped_count > max_skip_align )
      break;
    if ( !this->LocOverlap( loc1, loc2 ) ) {
      skipped_count++;
      continue;
    }
    if ( !this->AlignOverlap( loc1, loc2 ) )
      continue;
    return false;
  }
  
  return true;
}



/*
 * aligns_path
 * CountAligners
 *
 * It counts the number of reads aligning pos1 (as per locations file). If
 * aligner_ids is not null, it is filled with the read ids of all the reads
 * pos1 should align to.
 */
int aligns_path::CountAligners( int pos1,
				const vec<int> &loc_indices,
				vec<int> *aligner_ids ) const
{
  int count = 0;
  int max_skip_align = 3;
  int loc1 = loc_indices[pos1];
  if ( aligner_ids )
    aligner_ids->clear( );
  
  // On the left.
  int skipped_count = 0;
  for (int ii=pos1-1; ii>-1; ii--) {
    int loc2 = loc_indices[ii];
    if ( skipped_count > max_skip_align )
      break;
    if ( !this->LocOverlap( loc1, loc2 ) ) {
      skipped_count++;
      continue;
    }
    count++;
    if ( aligner_ids )
      aligner_ids->push_back( (*read_bags_)[loc2].ReadId( ) );
  }
  
  // On the right.
  skipped_count = 0;
  for (int ii=pos1+1; ii<(int)loc_indices.size( ); ii++) {
    int loc2 = loc_indices[ii];
    if ( skipped_count > max_skip_align )
      break;
    if ( !this->LocOverlap( loc1, loc2 ) ) {
      skipped_count++;
      continue;
    }
    count++;
    if ( aligner_ids )
      aligner_ids->push_back( (*read_bags_)[loc2].ReadId( ) );
  }

  return count;
}



/*
 * aligns_path
 * Setup
 */
void aligns_path::Setup( )
{
  max_align_score_ = 100.0;
  max_poly_score_ = 2;
  use_poly_score_ = ( bases_ && quals_ ) ? true : false;
}



/*
 * aligns_path
 * LocOverlap
 *
 * If the two location overlap (by read_bags_ only).
 */
bool aligns_path::LocOverlap( int loc1, int loc2 ) const
{
  const read_location &left = (*read_bags_)[loc1];
  const read_location &right = (*read_bags_)[loc2];
  
  if ( left.Contig( ) != right.Contig( ) )
    return false;
  if ( right.StartOnContig( ) < left.StartOnContig( ) )
    return ( left.StartOnContig( ) <= right.StopOnContig( ) );
  return ( right.StartOnContig( ) <= left.StopOnContig( ) );
}



/*
 * aligns_path
 * AlignOverlap
 *
 * If there is a good alignmnent between the two locations.
 */
bool aligns_path::AlignOverlap( int loc1, int loc2 ) const
{
  int al_index = FindAlign( loc1, loc2 );
  
  if ( al_index < 0 )
    return false;
  
  if ( (*aligns_)[al_index].score <= max_align_score_ )
    return true;
  
  if ( use_poly_score_ ) {
    ForceAssert( bases_ && quals_ );
    const alignment_plus &alplus = (*aligns_)[al_index];
    const align &al = align( alplus.a );
    bool isRC = alplus.Rc2( );
    int id1 = (*read_bags_)[loc1].ReadId( );
    int id2 = (*read_bags_)[loc2].ReadId( );
    const basevector &b1 = (*bases_)[id1];
    const basevector &b2 = (*bases_)[id2];
    const qualvector &q1 = (*quals_)[id1];
    const qualvector &q2 = (*quals_)[id2];
    
    int poly = ScoreAlignmentPoly( isRC, al, b1, q1, b2, q2 );
    
    if ( 0 < poly && poly <= max_poly_score_ )
      return true;
  }
  
  return false;
}



/*
 * aligns_path
 * FindAlign
 *
 * Returns either the index for the alignment between loc1 and loc2
 * or -1.
 */
int aligns_path::FindAlign( int loc1, int loc2 ) const
{
  int read1 = (*read_bags_)[loc1].ReadId( );
  int read2 = (*read_bags_)[loc2].ReadId( );
  
  int al_index = (*aligns_index_)[read1];
  if ( al_index < 0 )
    return -1;

  for (int ii=al_index; ii<(int)aligns_->size( ); ii++) {
    if ( (*aligns_)[ii].Id1( ) != read1 )
      break;
    if ( (*aligns_)[ii].Id2( ) == read2 )
      return ii;
  }

  return -1;
}



