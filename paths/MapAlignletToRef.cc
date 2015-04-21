///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "paths/MapAlignletToRef.h"

/**
 * MapAlignletToRef
 */
void MapAlignletToRef( const alignlet &al_onC,
		       const vec<look_align_plus> &cg_hits,
		       const vec< vec<int> > &to_hits,
		       vec<alignlet> &als_onR,
		       vec<int> *cg_hit_ids )
{
  als_onR.clear( );
  if ( cg_hit_ids ) cg_hit_ids->clear( );

  // Contig not aligned.
  const int contig_id = al_onC.TargetId( );
  const vec<int> &hit_ids = to_hits[contig_id];
  if ( hit_ids.size( ) < 1 ) return;
  
  int clen = al_onC.TargetLength( );
  int len2 = cg_hits[ hit_ids[0] ].QueryLength( );
  ForceAssertEq( clen, len2 );

  // Loop over all alignments of contigs containing read.
  for (int hid=0; hid<hit_ids.isize( ); hid++) {
    const look_align_plus &hit = cg_hits[ hit_ids[hid] ];
    ForceAssertEq( hit.a.Gaps( 0 ), 0 );
      
    // Do we need to reverse the read-on-contig alignment?
    alignlet al_onC_rc;
    if ( hit.Rc1( ) ) {
      al_onC_rc = al_onC;
      al_onC_rc.Reverse( );
    }
    const alignlet theAl = hit.Rc1( ) ? al_onC_rc : al_onC;

    // Read out of aligned portion of contig.
    int cg_beg = hit.StartOnQuery( );
    int cg_end = hit.EndOnQuery( );
    if ( theAl.pos2( ) < cg_beg || theAl.Pos2( ) > cg_end ) continue;

    // Find start of read-on-reference alignment.
    int posRef = hit.StartOnTarget( );
    int posCg = cg_beg;
    for (int block_id=0; block_id<hit.a.Nblocks( ); block_id++) {
      int block_gap = hit.a.Gaps( block_id );
      if ( block_gap >= 0 )
	posRef += block_gap;
      else {
	if ( posCg <= theAl.pos2( ) && theAl.pos2( ) < posCg - block_gap )
	  break;
	posCg += - block_gap;
      }
      int block_len = hit.a.Lengths( block_id );
      if ( posCg <= theAl.pos2( ) && theAl.pos2( ) < posCg + block_len ) {
	posRef += ( theAl.pos2( ) - posCg );
	break;
      }
      posRef += block_len;
      posCg += block_len;
    }
    ForceAssertLt( posRef, hit.EndOnTarget( ) );

    // Build alignlet, add it to result.
    int new_pos2 = posRef;
    int new_Pos2 = posRef + ( theAl.Pos2( ) - theAl.pos2( ) );
    int new_tid = hit.target_id;
    int new_tlen = hit.target_length;
    int new_fw1 = hit.rc1 ? ! al_onC.Fw1( ) : al_onC.Fw1( );
    alignlet new_al( new_pos2, new_Pos2, new_tid, new_tlen, new_fw1 );

    als_onR.push_back( new_al );
    if ( cg_hit_ids ) cg_hit_ids->push_back( hit_ids[hid] );
  }

}

