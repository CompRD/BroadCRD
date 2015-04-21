// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "Basevector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "lookup/LookAlign.h"
#include "tiled/PrimateContig.h"
#include "tiled/ReadsTiling.h"



/*
 * primate_contig
 * Constructor
 */
primate_contig::primate_contig( ) :
  clone_bases_ ( 0 ),
  primate_bases_ ( 0 ),
  primate_quals_ ( 0 ),
  log_ ( 0 )
{
  aligns_.resize( 0 );
}



/*
 * primate_contig
 * Constructor
 */
primate_contig::primate_contig( const basevector *clone_bases,
				const vecbasevector *primate_bases,
				const vecqualvector *primate_quals,
				ostream *log ) :
  clone_bases_ ( clone_bases ),
  primate_bases_ ( primate_bases ),
  primate_quals_ ( primate_quals ),
  log_ ( log )
{
  aligns_.resize( 0 );
}



/*
 * primate_contig
 * Set
 */
void primate_contig::Set( const basevector *clone_bases,
			  const vecbasevector *primate_bases,
			  const vecqualvector *primate_quals,
			  ostream *log )
{
  clone_bases_ = clone_bases;
  primate_bases_ = primate_bases;
  primate_quals_ = primate_quals;
  log_ = log;

  aligns_.resize( 0 );
  tiles_.Clear( );
}



/*
 * primate_contig
 * GenerateTiles
 */
void primate_contig::GenerateTiles( )
{
  vec<int> read_ids;
  vec<Bool> isRC;
  vec<align> plain_aligns;

  read_ids.reserve( aligns_.size( ) );
  isRC.reserve( aligns_.size( ) );
  plain_aligns.reserve( aligns_.size( ) );
  for (int ii=0; ii<(int)aligns_.size( ); ii++) {
    read_ids.push_back( aligns_[ii].query_id );
    isRC.push_back( aligns_[ii].rc1 );
    plain_aligns.push_back( aligns_[ii].a );
  }
  
  tiles_.SetPointers( clone_bases_, primate_bases_, primate_quals_ );
  tiles_.SetFromAligns( &read_ids, &plain_aligns, &isRC, log_ );

}



/*
 * primate_contig
 * Tiles
 */
const reads_tiling &primate_contig::Tiles( ) const
{
  return tiles_;
}



/*
 * primate_contig
 * AddAlignId
 */
void primate_contig::AddAlign( const look_align &new_align )
{
  aligns_.push_back( new_align );
}



/*
 * primate_contig
 * BulkBeginOnClone
 */
int primate_contig::BulkBeginOnClone( ) const
{
  if ( aligns_.size( ) < 1 )
    return 0;

  int leftmost = aligns_[0].a.pos2( );
  
  for (int ii=1; ii<(int)aligns_.size( ); ii++) {
    int al_begin = aligns_[ii].a.pos2( );
    if ( al_begin < leftmost )
      leftmost = al_begin;
  }

  return leftmost;
}



/*
 * primate_contig
 * BulkEndOnClone
 */
int primate_contig::BulkEndOnClone( ) const
{
  int rightmost = 0;

  for (int ii=0; ii<(int)aligns_.size( ); ii++) {
    int al_end = aligns_[ii].a.Pos2( );
    if ( al_end > rightmost )
      rightmost = al_end;
  }
  
  return rightmost;
}



/*
 * primate_contig
 * BulkLength
 */
int primate_contig::BulkLength( ) const
{
  int contig_end = this->BulkEndOnClone( );
  int contig_begin = this->BulkBeginOnClone( );

  return ( contig_end - contig_begin );
}



