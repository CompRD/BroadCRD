// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "math/Functions.h"
#include "String.h"
#include "tiled/PaddedSeq.h"
#include "tiled/TAlign.h"
#include "tiled/Tile.h"
#include "tiled/Tiling.h"



/*
 * tiling
 * Constructor
 */
tiling::tiling( ) :
  known_id_ ( -1 ),
  contig_id_ ( -1 ),
  n_tiles_ ( 0 )
{ }



/*
 * tiling
 * Constructor
 */
tiling::tiling( int known_id, int contig_id ) :
  known_id_ ( known_id ),
  contig_id_ ( contig_id ),
  n_tiles_ ( 0 )
{
  ForceAssert( known_id_ > -1 || contig_id_ > -1 );
}



/*
 * tiling
 * SetIds
 */
void tiling::SetIds( int known_id, int contig_id )
{
  known_id_ = known_id;
  contig_id_ = contig_id;
  n_tiles_ = 0;

  ForceAssert( known_id_ > -1 || contig_id_ > -1 );
}



/*
 * tiling
 * ResetContigId
 */
void tiling::ResetContigId( int contig_id )
{
  if ( contig_id_ == contig_id )
    return;

  contig_id_ = contig_id;
  contig_.SetContig( contig_id, contig_.RC( ), contig_.PaddedSeq( ) );
}



/*
 * tiling
 * SetFromAligns
 *
 * Return < 0 on error, or 1 if all right. It assumes that read_al (and
 * eventually cg_al) are the alignments of reads (eventually contig) onto
 * the known sequence.
 */
int tiling::SetFromAligns( const vec<t_align> &read_al, const t_align *cg_al )
{
  // No alignments in input.
  if ( read_al.size( ) < 1 && cg_al == 0 )
    return -1;

  // If both known_id_ and contig_id_ are given, then cg_al must be given too.
  if ( ( known_id_ > -1 && contig_id_ > -1 ) && !cg_al )
    return -2;

  // If cg_al != 0 then contig id from t_align must match contig_id_.
  if ( cg_al  && cg_al->id_ != contig_id_ )
    return -3;

  // Master sequence provides the coordinate system: its base 0 is the origin.
  n_tiles_ = 1 + ( cg_al ? 1 : 0 ) + (int)read_al.size( );
  int tiles_count = n_tiles_ - 1;
  
  vec<int> ends;
  ends.reserve( tiles_count );

  if ( cg_al ) {
    ends.push_back( cg_al->al_.Pos2( ) );
    tiles_count--;
  }
  for (int ii=0; ii<tiles_count; ii++)
    ends.push_back( read_al[ii].al_.Pos2( ) );

  int master_len = Max( ends );
  
  // Create tiles. At this time all padded_seq's are empty.
  if ( known_id_ > -1 ) {
    int begin = 0;
    int al_begin = 0;
    int len = master_len;
    padded_seq pads( begin, al_begin, len );

    known_.SetKnown( known_id_, pads );
  }

  if ( contig_id_ > -1 ) {
    if ( known_id_ > -1 ) {
      bool isRC = cg_al->isRC_;
      int begin = cg_al->al_.pos2( ) - cg_al->al_.pos1( );
      int al_begin = cg_al->al_.pos1( );
      int len = cg_al->al_.Pos1( ) - cg_al->al_.pos1( );
      padded_seq pads( begin, al_begin, len );

      contig_.SetContig( cg_al->id_, isRC, pads );
    }
    else {
      bool isRC = false;
      int begin = 0;
      int al_begin = 0;
      int len = master_len;
      padded_seq pads( begin, al_begin, len );
      
      contig_.SetContig( contig_id_, isRC, pads );
    }
  }

  reads_.resize( read_al.size( ) );
  for (int ii=0; ii<(int)read_al.size( ); ii++) {
    int id = read_al[ii].id_;
    bool isRC = read_al[ii].isRC_;
    const align &al = read_al[ii].al_;
    int begin = al.pos2( ) - al.pos1( );
    int al_begin = al.pos1( );
    int len = al.Pos1( ) - al.pos1( );
    padded_seq pads( begin, al_begin, len );

    reads_[ii].SetRead( id, isRC, pads );
  }
  
  // Adjust for gaps.
  if ( known_id_ > -1 && contig_id_ > -1 )
    this->Adjust( *cg_al, -1 );    
  
  for (int ii=0; ii<(int)reads_.size( ); ii++)
    this->Adjust( read_al[ii], ii );

  // Sort read tiles.
  sort( reads_.begin( ), reads_.end( ) );

  // Return ok.
  return 1;
}



/*
 * tiling
 * SetFromAlignsAlt
 *
 * Return < 0 on error, or 1 if all right. It assumes that the read_al's
 * represent the aligmnents of reads onto the contig, and cg_al the
 * alignment of the contig onto the known sequence.
 *
 * Warning! If contig aligns known as RC, then the alignments in read_al
 * are assumed to be those of the reads onto the rc of the contig.
 */
int tiling::SetFromAlignsAlt( const vec<t_align> &read_al,
			      const t_align &cg_al )
{
  // Swap known with contig
  t_align swap_al = cg_al;
  swap_al.al_.Flip( );
  swap_al.id_ = known_id_;
  swap_al.isRC_ = false;
  
  swap( known_id_, contig_id_ );
  swap( known_, contig_ );

  // Call SetFromAligns.
  this->SetFromAligns( read_al, &swap_al );

  // Swap back contig with known.
  swap( known_id_, contig_id_ );

  int old_len_k = known_.PaddedSeq( ).UnpaddedLength( );
  vec<int> vpads_k;
  vpads_k.reserve( known_.PaddedSeq( ).PadsCount( ) );
  for (int ii=0; ii<known_.PaddedSeq( ).PadsCount( ); ii++)
    vpads_k.push_back( known_.PaddedSeq( ).Pad( ii ) );

  int old_begin_c = contig_.PaddedSeq( ).Begin( );
  int old_al_begin_c = contig_.PaddedSeq( ).AlBegin( );
  int old_len_c = contig_.PaddedSeq( ).UnpaddedLength( );
  vec<int> vpads_c;
  vpads_c.reserve( contig_.PaddedSeq( ).PadsCount( ) );
  for (int ii=0; ii<contig_.PaddedSeq( ).PadsCount( ); ii++)
    vpads_c.push_back( contig_.PaddedSeq( ).Pad( ii ) );

  int new_len_k = old_al_begin_c + old_len_c;
  padded_seq new_pads_k( 0, 0, new_len_k, &vpads_c );
  known_.SetKnown( known_id_, new_pads_k );
  
  int new_begin_c = - old_begin_c;
  int new_al_begin_c = old_begin_c + old_al_begin_c;
  int new_len_c = old_len_k - ( old_begin_c + old_al_begin_c );
  padded_seq new_pads_c( new_begin_c, new_al_begin_c, new_len_c, &vpads_k );
  contig_.SetContig( contig_id_, cg_al.isRC_, new_pads_c );

  for (int ii=0; ii<(int)reads_.size( ); ii++) {
    int old_begin_r = reads_[ii].PaddedSeq( ).Begin( );
    reads_[ii].PaddedSeq( ).ResetBegin( new_begin_c + old_begin_r );
  }

  // Return ok.
  return 1;
}



/*
 * tiling
 * Shift
 *
 * Shifts origin of global coordinate system by given amount.
 */
void tiling::Shift( int amount )
{
  if ( known_id_ > -1 ) {
    known_.PaddedSeq( ).ShiftPadsAfter( 0, amount );
  }

  if ( contig_id_ > -1 ) {
    if ( known_id_ == -1 )
      contig_.PaddedSeq( ).ShiftPadsAfter( 0, amount );
    else {
      int new_begin = contig_.PaddedSeq( ).Begin( ) + amount;
      contig_.PaddedSeq( ).ResetBegin( new_begin );
    }
  }

  for (int ii=0; ii<(int)reads_.size( ); ii++) {
    int new_begin = reads_[ii].PaddedSeq( ).Begin( ) + amount;
    reads_[ii].PaddedSeq( ).ResetBegin( new_begin );
  }
}



/*
 * tiling
 * operator>>
 */
istream &operator>> ( istream &in, tiling &tiles )
{
  // known_id_, contig_id_, and n_tiles_.
  int known_id = 0;
  int contig_id = 0;
  int n_tiles = 0;
  
  in >> known_id;
  if ( !in )
    return in;
  in >> contig_id;
  if ( !in )
    return in;
  in >> n_tiles;
  if ( !in || n_tiles < 1 )
    return in;
  
  tiles.known_id_ = known_id;
  tiles.contig_id_ = contig_id;
  tiles.n_tiles_ = n_tiles;

  // The very first tile: either known or contig.
  tile master_tile;
  in >> master_tile;
  if ( master_tile.IsKnown( ) ) {
    ForceAssert( known_id == master_tile.Id( ) );
    tiles.known_ = master_tile;
  }
  else if ( master_tile.IsContig( ) ) {
    ForceAssert( contig_id == master_tile.Id( ) );
    tiles.contig_ = master_tile;
  }
  else
    ForceAssert( 1 == 0 );
  n_tiles--;

  // The next tile (first non-master tile): either contig or a read.
  tile first_tile;
  in >> first_tile;
  if ( first_tile.IsContig( ) ) {
    tiles.contig_ = first_tile;
    tiles.reads_.reserve( n_tiles - 1 );
  }
  else {
    tiles.reads_.reserve( n_tiles );
    tiles.reads_.push_back( first_tile );
  }
  n_tiles--;

  // The remaining reads.
  for (int ii=0; ii<n_tiles; ii++) {
    tile a_read_tile;
    in >> a_read_tile;
    tiles.reads_.push_back( a_read_tile );
  }

  return in;
}



/*
 * tiling
 * operator<<
 */
ostream &operator<< ( ostream &out, const tiling &tiles )
{
  out << tiles.known_id_ << "\t"
      << tiles.contig_id_ << "\t"
      << tiles.n_tiles_ << "\n";
  
  if ( !tiles.known_.IsEmpty( ) )
    out << tiles.known_;

  if ( !tiles.contig_.IsEmpty( ) )
    out << tiles.contig_;

  for (int ii=0; ii<(int)tiles.reads_.size( ); ii++)
    out << tiles.reads_[ii];

  return out;
}



/*
 * tiling
 * Adjust
 *
 * tile_id = -1 refers to contig_ aligning to master (in this case master
 * is known_). tile_id = ii, for ii >= 0, refers to reads_[ii] aligning to
 * master (in this case master can either be known_ or contig_).
 */
void tiling::Adjust( const t_align &tile_align, int tile_id )
{
  const align &the_al = tile_align.al_;

  // Consistency check.
  if ( tile_id == -1 && contig_id_ < 0 )
    return;
  
  tile &the_tile = ( tile_id == -1 ) ? contig_ : reads_[tile_id];
  ForceAssert( tile_align.id_ == the_tile.Id( ) );

  // Go at the beginning of alignment.
  padded_seq &the_pads = the_tile.PaddedSeq( );
  int n_blocks = (int)the_al.Nblocks( );
  int beg = the_pads.Begin( ) + the_al.pos1( );
  int pos = beg;

  // Walk along the alignment.
  for (int ii=0; ii<n_blocks; ii++) {

    // Gap on master sequence.
    if ( the_al.Gaps( ii ) < 0 )
      for (int jj=0; jj<-(int)the_al.Gaps(ii); jj++)
	this->PadClone( pos, tile_id );
    
    // Gap on aligned sequence.
    if ( the_al.Gaps( ii ) > 0 )
      for (int jj=0; jj<(int)the_al.Gaps(ii); jj++)
	this->PadRead( pos, tile_id );
    
    // Length.
    int len = the_al.Lengths( ii );
    
    if ( the_al.Gaps( ii ) < 0 )
      len += -the_al.Gaps( ii );

    while ( len > 0 ) {
      if ( pos >= beg + the_pads.PaddedLength( ) )
	break;
      if ( !the_pads.IsPad( pos ) )
	len--;
      pos++;
    }
  }
}



/*
 * tiling
 * PadClone
 */
void tiling::PadClone( int &pos, int tile_id )
{
  tile &the_tile = ( tile_id == -1 ) ? contig_ : reads_[tile_id];
  padded_seq &the_pads = the_tile.PaddedSeq( );

  if ( tile_id == -1 && contig_id_ < 0 )
    return;

  if ( the_pads.IsPad( pos ) ) {
    // read_index has bases extending into a master-padded region.
    int loc_pos = pos;
    while ( the_pads.IsPad( loc_pos ) )
      loc_pos++;
    loc_pos--;
    the_pads.RemovePad( loc_pos );
    the_pads.ShiftPadsAfter( loc_pos, +1 );
    
    // If last base is pad, remove it.
    int last_pos
      = the_pads.Begin( ) 
      + the_pads.AlBegin()
      + the_pads.PaddedLength()
      - 1;
    while ( the_pads.IsPad( last_pos ) ) {
      the_pads.RemovePad( last_pos );
      last_pos--;
    }
  }
  else {
    // Padding master.
    tile &master_tile = ( known_id_ > -1 ) ? known_ : contig_;
    padded_seq &master_pads = master_tile.PaddedSeq( );
    master_pads.AddPad( pos );
    
    // Push pads on tile_id mirroring pads on master after pos.
    the_pads.ShiftPadsAfter( pos, +1 );

    // Padding non-master sequences (with the exclusion of tile_id).
    for (int jj=-1; jj<(int)reads_.size( ); jj++) {
      if ( jj == -1 && known_id_ == -1 )
	continue;
      if ( jj == tile_id )
	continue;
      
      tile &other_tile = ( jj == -1 ) ? contig_ : reads_[jj];
      padded_seq &other_pads = other_tile.PaddedSeq( );

      // Reads other than read_index must be padded too.
      other_pads.AddPad( pos );
      
      // Adding a pad may push the read into a master padded region.
      int loc_pos
	= other_pads.Begin( )
	+ other_pads.AlBegin( )
	+ other_pads.PaddedLength( )
	- 1;
      while ( master_pads.IsPad( loc_pos ) ) {
	other_pads.AddPad( loc_pos );
	loc_pos++;
      }
    }
  }
}



/*
 * tiling
 * PadRead
 */
void tiling::PadRead( int &pos, int tile_id )
{
  tile &the_tile = ( tile_id == -1 ) ? contig_ : reads_[tile_id];
  tile &master_tile = ( known_id_ > -1 ) ? known_ : contig_;
  padded_seq &the_pads = the_tile.PaddedSeq( );
  padded_seq &master_pads = master_tile.PaddedSeq( );
  
  if ( tile_id == -1 && contig_id_ < 0 )
    return;

  // Move to the first non-padded master base.
  while ( master_pads.IsPad( pos ) )
    pos++;
  
  // Add the pad.
  the_pads.AddPad( pos );

  // Shift all pads after pos.
  the_pads.ShiftPadsAfter( pos, -1 );
  
  // Adding a pad may push the read into a master padded region.
  int loc_pos
    = the_pads.Begin( )
    + the_pads.AlBegin( )
    + the_pads.PaddedLength( )
    - 1;
  while ( master_pads.IsPad( loc_pos ) ) {
    the_pads.AddPad( loc_pos );
    loc_pos++;
  }

  // Move on.
  pos++;
}



/*
 * tiling
 * Compactify
 *
 * Remove pads-only columns.
 */
void tiling::Compactify( )
{
  tile &master_tile = ( known_id_ > -1 ) ? known_ : contig_;
  padded_seq &master_pads = master_tile.PaddedSeq( );

  // Loop over all master pads.
  for (int ii=0; ii<master_pads.PadsCount( ); ii++) {
    int pos = master_pads.Pad( ii );
    
    // Check if this is an all-pads column.
    bool all_pads = true;

    for (int jj=-1; jj<(int)reads_.size( ); jj++) {
      if ( jj == -1 && contig_id_ < 0 )
	continue;

      const tile &the_tile = ( jj == -1 ) ? contig_ : reads_[jj];
      const padded_seq &the_pads = the_tile.PaddedSeq( );
      
      int actual_begin = the_pads.Begin( ) + the_pads.AlBegin( );
      if ( pos < actual_begin )
	continue;
      if ( pos >= actual_begin + the_pads.PaddedLength( ) )
	continue;
      if ( !the_pads.IsPad( pos ) ) {
	all_pads = false;
	break;
      }
    }
  
    // Remove pads.
    if ( all_pads ) {
      master_pads.RemovePad( pos );
      for (int jj=-1; jj<(int)reads_.size( ); jj++) {
	if ( jj == -1 && contig_id_ < 0 )
	  continue;

	tile &the_tile = ( jj == -1 ) ? contig_ : reads_[jj];
	padded_seq &the_pads = the_tile.PaddedSeq( );
	the_pads.RemovePad( pos );
      }
    }
  }
}



