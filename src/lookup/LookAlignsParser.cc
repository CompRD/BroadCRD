// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "Basevector.h"
#include "Qualvector.h"
#include "ReadLocation.h"
#include "ScoreAlignment.h"
#include "STLExtensions.h"
#include "Vec.h"
#include "VecAlignmentPlus.h"
#include "lookup/LookAlign.h"
#include "lookup/LookAlignsParser.h"



/*
 * look_aligns_parser
 * Constructor
 */
look_aligns_parser::look_aligns_parser( ) :
  hits_ ( 0 ),
  aligns_ ( 0 ),
  bases_ ( 0 ),
  quals_ ( 0 )
{ }



/*
 * look_aligns_parser
 * Constructor
 */
look_aligns_parser::look_aligns_parser( const vec<look_align_plus> *hits,
					const vec_alignment_plus *aligns,
					const vecbasevector *bases,
					const vecqualvector *quals ) :
  hits_ ( hits ),
  aligns_ ( aligns ),
  bases_ ( bases ),
  quals_ ( quals ) 
{ }



/*
 * look_aligns_parser
 * SetPointers
 */
void look_aligns_parser::SetPointers( const vec<look_align_plus> *hits,
				      const vec_alignment_plus *aligns,
				      const vecbasevector *bases,
				      const vecqualvector *quals )
{
  hits_ = hits;
  aligns_ = aligns;
  bases_ = bases;
  quals_ = quals;
}



/*
 * look_aligns_parser
 * PartitionHits
 *
 * It partitions the hits in bags (proto contigs), saving in bag_begins the
 * position of the first alignment of each bag.
 */
void look_aligns_parser::PartitionHits( vec<int> &bag_begins ) const
{
  if ( hits_->empty() )
      return;

  // Look aligns must be sorted.
  order_lookalign_TargetLookAlignOffset order;
  ForceAssert( is_sorted( hits_->begin( ), hits_->end( ), order ) );

  // Initialize bag_begins.
  bag_begins.clear( );
  bag_begins.push_back( 0 );

  // Loop over the hits.
  for (int al_pos=1; al_pos<(int)(*hits_).size( ); al_pos++) {
    int this_bag_begin = bag_begins[ bag_begins.size( ) - 1 ];
    bool overlap_found = false;
    
    for (int jj=al_pos-1; jj>=this_bag_begin; jj--) {
      // Overlap on master.
      if ( LookAlignOffsetOverlap( (*hits_)[jj], (*hits_)[al_pos] ) > 0 ) {
	overlap_found = true;
	break;
      }
      
      // Overlap by align.
      if ( this->AlignId( jj, al_pos ) > -1 ) {
	overlap_found = true;
	break;
      }
    }
    
    if ( !overlap_found )
      bag_begins.push_back( al_pos );
  }
  
}



/*
 * look_aligns_parser
 * AppendLocs
 *
 * After having partitioned the hits in bags, it builds an initial
 * approximation for the vector of read locations.
 */
void look_aligns_parser::AppendLocs( vec<read_location> &locs ) const
{
  locs.reserve( hits_->size( ) );

  vec<int> bag_begins;
  this->PartitionHits( bag_begins );
  
  for (int bag_id=0; bag_id<(int)bag_begins.size( ); bag_id++) {
    int contig_id = locs.size( ) == 0 ? 0 : locs[locs.size()-1].Contig( ) + 1;
    int begin = bag_begins[bag_id];
    int end = ( bag_id < (int)bag_begins.size( ) - 1 )
      ? bag_begins[bag_id+1]
      : (int)hits_->size( );

    if ( end - begin < 1 )
      continue;

    // Build a first approximation for temp_locs.
    vec<read_location> temp_locs;
    temp_locs.reserve( end - begin );
    for (int al_pos=begin; al_pos<end; al_pos++) {
      const look_align &the_al = (*hits_)[al_pos];
      int id = the_al.query_id;
      int len = the_al.query_length;
      int contig_len = 0;      // temp value, filled below.

      // The first read in the contig starts at 0.
      int start = 0;
      bool start_found = ( al_pos == begin ) ? true : false;

      // For the successive reads, start is determined by the offset
      //  with the closest already placed read, where the offset is
      //  calculated either with LookAlignOffset (in LookAlign.h), or by using
      //  the align between the two reads.
      int al_pos_pre = al_pos-1;
      while ( al_pos_pre >= begin && start_found == false ) {
	int offset;
	if( LookAlignOffset( (*hits_)[al_pos_pre], (*hits_)[al_pos], offset ) ||
	    this->AlignOffset( al_pos_pre, al_pos, offset ) ) {
	  int loc_id = -1;
	  for (int ii=temp_locs.size( )-1; ii>=0; ii--) {
	    if ( temp_locs[ii].ReadId( ) == (*hits_)[al_pos_pre].query_id ) {
	      loc_id = ii;
	      break;
	    }
	  }
	  ForceAssert( loc_id > -1 );
	  start = temp_locs[loc_id].StartOnContig( ) + offset;
	  start_found = true;
	  break;
	}
	al_pos_pre--;
      }
      ForceAssert( start_found );

      orientation orient = the_al.rc1 ? ReverseOr : ForwardOr;
      read_location new_loc(id, len, contig_id, start, orient, contig_len);
      temp_locs.push_back( new_loc );
    }

    // Adjust temp_locs (fix contig lengths);
    ForceAssert( temp_locs.size( ) > 0 );

    sort( temp_locs.begin( ), temp_locs.end( ) );
  
    int cg_len = temp_locs[0].LengthOfRead( );
    for (int ii=1; ii<(int)temp_locs.size( ); ii++) {
      int start = temp_locs[ii].StartOnContig( );
      int len = temp_locs[ii].LengthOfRead( );
      int read_end = start + len;
      if ( read_end > cg_len )
	cg_len = read_end;
    }
    for (int ii=0; ii<(int)temp_locs.size( ); ii++)
      temp_locs[ii].SetLengthOfContig( cg_len );
    
    // Append temp_locs to locs.
    copy( temp_locs.begin( ), temp_locs.end( ), back_inserter( locs ) );
  }
  
}



/*
 * look_aligns_parser
 * GetBracket
 *
 * Given a set of read ids, it returns the smallest bracket [begin, end) on
 * the master sequence containing all read ids.
 */
pair<int, int> look_aligns_parser::GetBracket( const vec<int> &read_ids ) const
{
  if ( id_to_hit_.size( ) < 1 )
    this->GenerateIdToHit( );
  
  if ( read_ids.size( ) < 1 )
    return make_pair( -1, -1 );
  
  int leftmost = -1;
  int rightmost = -1;
  map<int, int>::const_iterator iter;
  for (int ii=0; ii<(int)read_ids.size( ); ii++) {
    iter = id_to_hit_.find( read_ids[ii] );
    if ( iter == id_to_hit_.end( ) )
      continue;
    const look_align &hit = (*hits_)[iter->second];
    int hit_begin = hit.a.pos2( ) - hit.a.pos1( );
    int hit_end = hit_begin + hit.query_length;
    leftmost = ( leftmost < 0 ) ? hit_begin : Min( leftmost, hit_begin );
    rightmost = Max( rightmost, hit_end );
  }
  
  if ( leftmost < 0 || rightmost < 0 )
    return make_pair( -1, -1 );
  
  return make_pair( leftmost, rightmost );
}



/*
 * look_aligns_parser
 * GenerateIdToHit
 */
void look_aligns_parser::GenerateIdToHit( ) const
{
  for (int ii=0; ii<(int)hits_->size( ); ii++)
    id_to_hit_[ (*hits_)[ii].query_id ] = ii;
}



/*
 * look_aligns_parser
 * AlignId
 *
 * If aligns_ is not null, check if the two reads share a good align. It
 * returns either the good align id or -1.
 */
int look_aligns_parser::AlignId( int hit1, int hit2 ) const
{
  if ( !aligns_ )
    return -1;
  
  int read_id1 = (*hits_)[hit1].query_id;
  int read_id2 = (*hits_)[hit2].query_id;
  
  int al_idx = aligns_->GetAlignsIndex( read_id1 );
  if ( al_idx < 0 )
    return -1;
  for (int ii=al_idx; ii<(int)aligns_->GetNumberAlignments( ); ii++) {
    if ( aligns_->GetAlignmentId1( ii ) != read_id1 )
      break;
    if ( aligns_->GetAlignmentId2( ii ) == read_id2 )
      return ( this->IsGoodAlign( ii ) ) ? ii : -1;
  }
  
  return -1;
}



/*
 * look_aligns_parser
 * AlignOffset
 *
 * If aligns_ is not null, fill offset with the offset induced by the
 * alignment between the two reads. Returns false if the reads do not
 * align.
 */
bool look_aligns_parser::AlignOffset( int hit1, int hit2, int &offset ) const
{
  offset = 0;
  
  if ( !aligns_ )
    return false;
  
  int align_index = this->AlignId( hit1, hit2 );
  if ( align_index < 0 )
    return false;
  
  alignment_plus the_ap;
  aligns_->GetAlignment( the_ap, align_index );
  const alignment &al = the_ap.a;

  offset = al.pos1( ) - al.pos2( );
  return true;
}



/*
 * look_aligns_parser
 * IsGoodAlign
 */
bool look_aligns_parser::IsGoodAlign( int align_index ) const
{
  if ( !aligns_ )
    return false;

  // Some heuristic.
  float very_good = 100.0;
  float barely_good = 1000.0;
  int max_poly_score = 4;
  
  float align_score = aligns_->GetAlignmentScore( align_index );
  
  // Very good align.
  if ( align_score <= very_good )
    return true;

  // Align is not acceptable.
  if ( align_score > barely_good )
    return false;

  // Align may be salvaged, provided poly_score is low.
  if ( !bases_ || !quals_ )
    return false;
  static alignment_plus alplus;
  aligns_->GetAlignment( alplus, align_index );
  const align &al = align( alplus.a );
  bool isRC = alplus.Rc2( );
  int id1 = aligns_->GetAlignmentId1( align_index );
  int id2 = aligns_->GetAlignmentId2( align_index );
  const basevector &b1 = (*bases_)[id1];
  const basevector &b2 = (*bases_)[id2];
  const qualvector &q1 = (*quals_)[id1];
  const qualvector &q2 = (*quals_)[id2];
  
  int poly = ScoreAlignmentPoly( isRC, al, b1, q1, b2, q2 );
  
  return ( 0 < poly && poly <= max_poly_score );
    
}



