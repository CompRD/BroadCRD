// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "tiled/ChrMapper.h"


/*
 * chr_mapper
 * Constructor
 */
chr_mapper::chr_mapper( ) :
  chr_ ( -1 ),
  id_ ( -1 ),
  pos_ ( -1 ),
  begin_ ( -1 )
{ }



/*
 * chr_mapper
 * Constructor
 */
chr_mapper::chr_mapper( int chr, int id, int pos, int begin ) :
  chr_ ( chr ),
  id_ ( id ),
  pos_ ( pos ),
  begin_ ( begin )
{ }



/*
 * chr_mapper
 * Set
 */
void chr_mapper::Set( int chr, int id, int pos, int begin )
{
  chr_ = chr;
  id_ = id;
  pos_ = pos;
  begin_ = begin;
}



/*
 * chr_mapper
 * operator<
 *
 * Sort by chromosome id and by position on chromosome.
 */
bool operator< ( const chr_mapper &left, const chr_mapper &right )
{
  if ( left.chr_ == right.chr_ )
    return ( left.pos_ < right.pos_ );
  return ( left.chr_ < right.chr_ );
}



/*
 * chr_mapper
 * operator>>
 */
istream &operator>> ( istream &in, chr_mapper &cmap )
{
  in >> cmap.chr_
     >> cmap.id_
     >> cmap.pos_
     >> cmap.begin_;
  
  return in;
}



/*
 * chr_mapper
 * operator<<
 */
ostream &operator<< ( ostream &out, const chr_mapper &cmap )
{
  out << cmap.chr_ << "\t"
      << cmap.id_ << "\t"
      << cmap.pos_ << "\t"
      << cmap.begin_ << "\n";

  return out;
}



