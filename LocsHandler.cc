// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "Loader.h"
#include "LocsHandler.h"
#include "ReadLocation.h"
#include "String.h"

/*
 * lhandler
 * Constructor
 */
lhandler::lhandler( int n_reads, int n_contigs ) :
  n_reads_ ( n_reads ),
  n_contigs_ ( n_contigs )
{ }

/*
 * lhandler
 * Constructor
 */
lhandler::lhandler( int n_reads, int n_contigs, const String &locs_file ) :
  n_reads_ ( n_reads ),
  n_contigs_ ( n_contigs )
{
  this->LoadFromFile( locs_file );
}

/*
 * lhandler
 * LoadFromFile
 */
void lhandler::LoadFromFile( const String &locs_file )
{
  READX( locs_file, locs_ );
  sort( locs_.begin( ), locs_.end( ) );
  
  LoadFirstLocs( locs_, flocs_ );

  to_loc_.clear( );
  nreads_.clear( );
  to_loc_.resize( n_reads_ );
  nreads_.resize( n_contigs_, 0 );
  for (int ii=0; ii<(int)locs_.size( ); ii++) {
    to_loc_[locs_[ii].ReadId( )].push_back( ii );
    nreads_[locs_[ii].Contig( )] += 1;
  }
}

/*
 * lhandler
 * SetFromVec
 *
 * Warning: it will copy the original locs (i.e., not space efficient),
 * and it will deduce n_contigs_ from the given vec of read_locations.
 */
void lhandler::SetFromVec( const vec<read_location> &in_locs )
{
  locs_ = in_locs;
  sort( locs_.begin( ), locs_.end( ) );
  
  LoadFirstLocs( locs_, flocs_ );
  n_contigs_ = (int)flocs_.size( );

  to_loc_.clear( );
  nreads_.clear( );
  to_loc_.resize( n_reads_ );
  nreads_.resize( n_contigs_, 0 );
  for (int ii=0; ii<(int)locs_.size( ); ii++) {
    to_loc_[locs_[ii].ReadId( )].push_back( ii );
    nreads_[locs_[ii].Contig( )] += 1;
  }
}

/*
 * lhandler
 * PosInContig
 */
int lhandler::PosInContig( int read_id, int contig_id ) const
{
  int loc_id = -1;
  const vec<int> &to_loc = to_loc_[read_id];
  for (int ii=0; ii<(int)to_loc.size( ); ii++) {
    if ( locs_[to_loc[ii]].Contig( ) == contig_id ) {
      loc_id = to_loc[ii];
      break;
    }
  }
  if ( loc_id < 0 )
    return -1;

  int floc = flocs_[contig_id];
  if ( floc > loc_id )
    return -1;
  
  return ( loc_id - floc );
}

/*
 * lhandler
 * Coverage
 */
void lhandler::Coverage( longlong &total_rlen, longlong &total_clen ) const
{
  total_rlen = 0;
  for (int ii=0; ii<(int)locs_.size( ); ii++)
    total_rlen += locs_[ii].LengthOfRead( );

  total_clen = 0;
  for (int ii=0; ii<n_contigs_; ii++) {
    int floc = flocs_[ii];
    if ( floc > -1 )
      total_clen += locs_[floc].LengthOfContig( );
  }
}

/*
 * lhandler
 * GetPlacement
 */
const read_location *lhandler::GetPlacement( int read_id ) const
{
  if ( to_loc_[read_id].size( ) != 1 )
    return 0;
  return &( locs_[ to_loc_[read_id][0] ] );
}

