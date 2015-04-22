// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef TILE_H
#define TILE_H

#include "PackAlign.h"
#include "String.h"
#include "tiled/PaddedSeq.h"



/*
 * class tile
 *
 * A tile is a way to store informations about a padded_seq. There are
 * three types of tile: a known sequence tile, a consensus contig tile,
 * and a plain read tile.
 */
class tile {
  
public:
  
  tile( ) : name_ ( "" ) { }

  void SetKnown( int id, const padded_seq &pads );
  
  void SetContig( int id, bool isRC, const padded_seq &pads );
  
  void SetRead( int id, bool isRC, const padded_seq &pads );

  bool IsEmpty( ) const { return ( name_ == "" ); }

  bool IsKnown( ) const { return ( name_.Contains( "k", 0 ) ); }
  
  bool IsContig( ) const { return ( name_.Contains( "c", 0 ) ); }
  
  bool IsRead( ) const { return ( name_.Contains( "r", 0 ) ); }

  bool RC( ) const { return ( name_.Contains( "_rc", -1 ) ); }

  int Id( ) const;
  
  String PrettyName( int width = 18 ) const;

  padded_seq &PaddedSeq( ) { return pads_; }

  const padded_seq &PaddedSeq( ) const { return pads_; }

  friend bool operator< ( const tile &left, const tile &right );

  friend istream &operator>> ( istream &in, tile &the_tile );
  
  friend ostream &operator<< ( ostream &out, const tile &the_tile );
  
  
private:
  
  String name_;      // contains info about type, id, and orientation
  padded_seq pads_;  // the pads
  
};



#endif
