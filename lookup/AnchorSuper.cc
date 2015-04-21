// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "LocsHandler.h"
#include "SupersHandler.h"
#include "lookup/AnchorSuper.h"
#include "lookup/LookAlign.h"

bool AnchorSuper( const look_align_plus &hit,
		  const lhandler &locs,
		  const shandler &supers,
		  read_location &out_loc )
{
  int read_id = hit.query_id;
  
  // Read is either unplaced or multiply placed.
  const read_location *loc = locs.GetPlacement( read_id );
  if ( !loc )
    return false;
  
  // Find which super the read belongs to.
  int contig_id = loc->Contig( );
  int super_id = supers.ToSuper( contig_id );
  int super_len =  supers[super_id].FullLength( );
  
  // Orientation of super on target.
  int on_contig = loc->OrientationOnContig( ) == ReverseOr ? 1 : 0;
  int on_target = hit.rc1 ? 1 : 0;
  orientation orient = on_target == on_contig ? ForwardOr : ReverseOr;
  
  // Implied start of super on target.
  int rd_start = loc->StartOnContig( ) + supers.StartOnSuper( contig_id );
  int rd_end = super_len - rd_start - hit.query_length;
  int actual_start = orient == ForwardOr ? rd_start : rd_end;
  int implied_start = hit.a.pos2( ) - hit.a.pos1( ) - actual_start;
  
  // Fill out_loc and return.
  out_loc.SetReadId( super_id );
  out_loc.SetLengthOfRead( super_len );
  out_loc.SetContig( hit.target_id );
  out_loc.SetStartOnContig( implied_start );
  out_loc.SetOrientationOnContig( orient );
  out_loc.SetLengthOfContig( hit.target_length );
  out_loc.SetOnSupers( False );

  return true;
}
