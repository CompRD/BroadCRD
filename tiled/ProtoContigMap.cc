// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include <map>

#include "system/System.h"
#include "Vec.h"
#include "tiled/ProtoContigId.h"
#include "tiled/ProtoContigMap.h"



/*
 * proto_contig_map
 * Constructor
 */
proto_contig_map::proto_contig_map( ) :
  checked_ ( false )
{ }



/*
 * proto_contig_map
 * Clear
 */
void proto_contig_map::Clear( ) 
{
  proto_id_.clear( );
  checked_ = false;
}



/*
 * proto_contig_map
 * Reserve
 */
void proto_contig_map::Reserve( int size ) 
{
  proto_id_.reserve( size );
}



/*
 * proto_contig_map
 * Add
 */
void proto_contig_map::Add( const proto_contig_id &proto_id )
{
  proto_id_.push_back( proto_id );
  checked_ = false;
}



/*
 * proto_contig_map
 * Add
 */
void proto_contig_map::Add( int cg_id, int chr_id, int pos )
{
  proto_id_.push_back( proto_contig_id( cg_id, chr_id, pos ) );
  checked_ = false;
}



/*
 * proto_contig_map
 * Find
 *
 * Find the proto_id with given chr_id and chr_pos.
 */
int proto_contig_map::Find( int chr_id, int pos )
{
  if ( !checked_ )
    this->CheckData( );
  
  order_protocontigid_ChrPos order;
  proto_contig_id temp_pci( -1, chr_id, pos );
  vec<proto_contig_id>::iterator lbound
    = lower_bound( proto_id_.begin( ), proto_id_.end( ), temp_pci, order );

  // Not found.
  if ( lbound == proto_id_.end( ) ||
       chr_id != lbound->ChrId( ) ||
       pos != lbound->ChrPos( ) )
    return -1;

  // Found.
  return distance( proto_id_.begin( ), lbound );

}



/*
 * proto_contig_map
 * Save
 */
void proto_contig_map::Save( const String &out_file )
{
  if ( !checked_ )
    this->CheckData( );
  
  WRITE( out_file, proto_id_ );
}



/*
 * proto_contig_map
 * Load
 */
void proto_contig_map::Load( const String &in_file )
{
  READX( in_file, proto_id_ );
  this->CheckData( );
}



/*
 * proto_contig_map
 * Size
 */
int proto_contig_map::Size( ) const
{
  return (int)proto_id_.size( );
}



/*
 * proto_contig_map
 * operator[]
 */
const proto_contig_id &proto_contig_map::operator[] ( int jj ) const
{
  return proto_id_[jj];
}



/*
 * proto_contig_map
 * CheckData
 *
 * Sort and unique.
 */
void proto_contig_map::CheckData( )
{
  if ( proto_id_.size( ) < 1 )
    return;
  
  order_protocontigid_ChrPos order;
  sort( proto_id_.begin( ), proto_id_.end( ), order );
  proto_id_.erase( unique( proto_id_.begin( ), proto_id_.end( ) ),
		   proto_id_.end( ) );
  
  checked_ = true;
}



