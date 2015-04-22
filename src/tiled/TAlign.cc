// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "PackAlign.h"
#include "tiled/TAlign.h"
#include "tiled/Tile.h"



/*
 * t_align
 * Write
 */
void t_align::Write( ostream &out ) const
{
  out << id_ << "\t"
      << ( isRC_ ? int( 1 ) : int( 0 ) ) << "\t"
      << al_.pos1( ) << "\t"
      << al_.pos2( ) << "\t"
      << al_.Nblocks( ) << "\t";

  for (int ii=0; ii<al_.Nblocks( ); ii++)
    out << al_.Gaps( ii ) << "\t" << al_.Lengths( ii ) << "\t";

  out << "\n";
}



/*
 * TAlign
 * Read
 */
void t_align::Read( istream &in )
{
  int orient;
  int pos1;
  int pos2;
  int n_blocks;

  in >> id_
     >> orient
     >> pos1
     >> pos2
     >> n_blocks;

  isRC_ = ( orient == 1 ) ? True : False;

  avector<int> gaps;
  avector<int> lengths;
  gaps.Setsize( n_blocks );
  lengths.Setsize( n_blocks );
  for (int ii=0; ii<n_blocks; ii++)
    in >> gaps(ii) >> lengths(ii);

  al_.Set( pos1, pos2, gaps, lengths );
}



/*
 * t_align
 * WriteBinary
 */
void t_align::WriteBinary( ostream &out ) const
{
  int orient = isRC_ ? 1 : 0;
  int pos1 = al_.pos1( );
  int pos2 = al_.pos2( );
  int n_blocks = al_.Nblocks( );

  out.write( (char *) &( id_ ), sizeof( int ) );
  out.write( (char *) &( orient ), sizeof( int ) );
  out.write( (char *) &( pos1 ), sizeof( int ) );
  out.write( (char *) &( pos2 ), sizeof( int ) );
  out.write( (char *) &( n_blocks ), sizeof( int ) );

  for (int ii=0; ii<al_.Nblocks( ); ii++) {
    int gap = al_.Gaps(ii);
    int length = al_.Lengths(ii);
    out.write( (char *) &( gap ), sizeof( int ) );
    out.write( (char *) &( length ), sizeof( int ) );
  }
}



/*
 * t_align
 * ReadBinary 
 */
void t_align::ReadBinary( istream &in )
{
  int orient;
  int pos1;
  int pos2;
  int n_blocks;

  in.read( (char *) &( id_ ), sizeof( int ) );
  in.read( (char *) &( orient ), sizeof( int ) );
  in.read( (char *) &( pos1 ), sizeof( int ) );
  in.read( (char *) &( pos2 ), sizeof( int ) );
  in.read( (char *) &( n_blocks ), sizeof( int ) );

  isRC_ = ( orient == 1 ) ? True : False;

  avector<int> gaps;
  avector<int> lengths;
  gaps.Setsize( n_blocks );
  lengths.Setsize( n_blocks );
  for (int ii=0; ii<n_blocks; ii++) {
    in.read( (char *) &( gaps(ii) ), sizeof( int ) );
    in.read( (char *) &( lengths(ii) ), sizeof( int ) );
  }

  al_.Set( pos1, pos2, gaps, lengths );
}



