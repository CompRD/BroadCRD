///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/TranslateAligns.h"

/**
 * TranslateAlignsToSupers
 */
void TranslateAlignsToSupers( const vec<superb> &supers,
			      vec<alignlet> &aligns )
{
  // Generate maps.
  vec<int> start_on_super;   // key: contig_id
  vec<int> to_super;         // key: contig_id
  vec<int> super_len;        // key: super_id

  int n_contigs = 0;
  for (int super_id=0; super_id<supers.isize( ); super_id++)
    n_contigs += supers[super_id].Ntigs( );

  start_on_super.resize( n_contigs, -1 );
  to_super.resize( n_contigs, -1 );
  super_len.reserve( supers.size( ) );
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    const superb &sup = supers[super_id];
    super_len.push_back( sup.FullLength( ) );

    int pos = 0;
    for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
      int contig_id = sup.Tig( cgpos );
      start_on_super[contig_id] = pos;
      to_super[contig_id] = super_id;
      pos += sup.Len( cgpos );
      if ( cgpos < sup.Ntigs( ) - 1 ) pos += sup.Gap( cgpos );
    }
  }

  // Translate aligns.
  for (size_t align_id=0; align_id<aligns.size( ); align_id++) {
    alignlet &al = aligns[align_id];
    int contig_id = al.TargetId( );
    bool fw = al.Fw1( );
    ForceAssert( to_super[contig_id] > -1 );
    int super_id = to_super[contig_id];
    al.SetTargetId( super_id );
    al.SetTargetLength( super_len[super_id], fw ); 
    al.Shift( start_on_super[contig_id] );
  }
  
}
