// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "String.h"
#include "tiled/PaddedSeq.h"
#include "tiled/Tile.h"



/*
 * tile
 * SetKnown
 */
void tile::SetKnown( int id, const padded_seq &pads )
{
  name_ = "k" + ToString( id );
  pads_ = pads;
}



/*
 * tile
 * SetContig
 */
void tile::SetContig( int id, bool isRC, const padded_seq &pads )
{
  String str_rc = isRC ? "_rc" : "";
  name_ = "c" + ToString( id ) + str_rc;
  pads_ = pads;
}



/*
 * tile
 * SetRead
 */
void tile::SetRead( int id, bool isRC, const padded_seq &pads )
{
  String str_rc = isRC ? "_rc" : "";
  name_ = "r" + ToString( id ) + str_rc;
  pads_ = pads;
}



/*
 * tile
 * Id
 *
 * Returns -1 if uninitialized.
 */
int tile::Id( ) const
{
  if ( name_.size( ) < 1 )
    return -1;
  
  String loc_name = name_.substr( 1, name_.size( ) - 1 );
  if ( loc_name.Contains( "_rc", -1 ) )
    loc_name = loc_name.Before( "_rc" );
  
  return loc_name.Int( );
}



/*
 * tile
 * PrettyName
 *
 * Generate printable name. If width > 0, then either add enough
 * spaces at the end of the name, or cut the name, in such a way the
 * returned String has size exactly width.
 */
String tile::PrettyName( int width ) const
{
  String read_name = name_;

  if ( read_name.Contains( "k", 0 ) )
    read_name.replace( 0, 1, "known_" );
  else if ( read_name.Contains( "c", 0 ) )
    read_name.replace( 0, 1, "contig_" );
  else if ( read_name.Contains( "r", 0 ) )
    read_name = read_name.After( "r" );
  else
    read_name = "";

  if ( width < 1 )
    return read_name;
  
  int old_size = read_name.size( );
  read_name.resize( width );
  for (int ii=old_size; ii<(int)read_name.size( ); ii++)
    read_name[ii] = ' ';

  return read_name;
}



/*
 * tile
 * operator<
 *
 * Sort by Begin.
 */
bool operator< ( const tile &left, const tile &right )
{
  return ( left.pads_.Begin( ) < right.pads_.Begin( ) );
}



/*
 * tile
 * operator>>
 */
istream &operator>> ( istream &in, tile &the_tile )
{
  in >> the_tile.name_
     >> the_tile.pads_;
  
  return in;
}



/*
 * tile
 * operator<<
 */
ostream &operator<< ( ostream &out, const tile &the_tile )
{
  out << the_tile.name_ << "\t"
      << the_tile.pads_ << "\n";
  
  return out;
}



