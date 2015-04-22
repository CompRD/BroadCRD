/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "STLExtensions.h"
#include "String.h"
#include "Vec.h"
#include "system/Assert.h"
#include "tiled/Tiling.h"
#include "tiled/Tiling.h"

/**
 * LoadVecTiling
 */
void LoadVecTiling( const String &tilings_file,
		    vec<tiling> &tiles,
		    const vec<int> *contig_ids )
{
  if ( contig_ids )
    ForceAssert( is_sorted( contig_ids->begin( ), contig_ids->end( ) ) );

  ifstream in_tiles( tilings_file.c_str( ) );
  while ( in_tiles ) {
    tiling input_tile;
    in_tiles >> input_tile;
    if ( !in_tiles )
      break;
    if ( contig_ids ) {
      int cid = input_tile.ContigId( );
      if ( ! binary_search( contig_ids->begin( ), contig_ids->end( ), cid ) )
	continue;
    }
    tiles.push_back( input_tile );
  }
  in_tiles.close( );
}

/**
 * SaveVecTiling
 */
void SaveVecTiling( const String &tilings_file, const vec<tiling> &tiles )
{
  ofstream out_tiles( tilings_file.c_str( ) );
  for (int ii=0; ii<(int)tiles.size( ); ii++)
    out_tiles << tiles[ii] << "\n";
  out_tiles.close( );
}



