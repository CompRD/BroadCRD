/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "PackAlign.h"
#include "String.h"
#include "Vec.h"
#include "tiled/TAlign.h"
#include "tiled/VecTAligns.h"

/**
 * BinWriteTAligns
 */
void BinWriteTAligns( ostream &out, int id, const vec<t_align> &t_aligns )
{
  int n_tiles = t_aligns.size( );

  out.write( (char *) &( id ), sizeof( int ) );
  out.write( (char *) &( n_tiles ), sizeof( int ) );

  for (int ii=0; ii<n_tiles; ii++)
    t_aligns[ii].WriteBinary( out );
}

/**
 * BinReadTAligns
 */
void BinReadTAligns( istream &in, int &id, vec<t_align> &t_aligns )
{
  int n_tiles = -1;

  in.read( (char *) &( id ), sizeof( int ) );
  in.read( (char *) &( n_tiles ), sizeof( int ) );

  if ( !in ) {
    id = -1;
    t_aligns.clear( );
    return;
  }

  t_aligns.resize( n_tiles );
  for (int ii=0; ii<n_tiles; ii++)
    t_aligns[ii].ReadBinary( in );
}

/**
 * LoadAllTAlignPluses
 */
void LoadAllTAlignPluses( const String &in_file, vec<t_align_plus> &t_aplus )
{
  t_aplus.clear( );

  ifstream in( in_file.c_str( ) );
  while ( in ) {
    int contig_id;
    vec<t_align> aligns;
    BinReadTAligns( in, contig_id, aligns );
    if ( !in || contig_id < 0 )
      break;
    t_aplus.reserve( t_aplus.size( ) + aligns.size( ) );
    for (int ii=0; ii<(int)aligns.size( ); ii++) {
      t_align_plus tap( aligns[ii], contig_id );
      t_aplus.push_back( tap );
    }
  }
  in.close( );

}

