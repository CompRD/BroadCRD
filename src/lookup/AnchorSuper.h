// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Anchor a supercontig to a target genome, by using a given look_align.

#ifndef ANCHOR_SUPER_H
#define ANCHOR_SUPER_H

#include "LocsHandler.h"
#include "SupersHandler.h"
#include "lookup/LookAlign.h"

/*
 * AnchorSuper
 *
 * The given look_align (hit) represents the alignment of a read onto a
 * target genome. If the read has been assembled, and it is uniquely
 * placed, then it is possible to anchor the supercontig to the target
 * genome by using the implied positioning given by hit.
 *
 * Returns true iff the read at hit is uniquely placed in the assembly
 * described by locs and supers (in which case it will also fill out_loc
 * with the implied position of the supercontig onto the target genome).
 *
 * Notice that "read" in out_loc stores the supercontig data, and "contig"
 * the target_id data (e.g. out_loc.ReadId( ) is the id of the supercontig
 * the read at hit belongs to).
 *
 * hit: alignment of a read onto a target genome
 * locs: locs handler (see LocsHandler.h)
 * supers: supers handler (see SupersHandler.h)
 * out_loc: output, as a read_location
 */
bool AnchorSuper( const look_align_plus &hit,
		  const lhandler &locs,
		  const shandler &supers,
		  read_location &out_loc );

#endif
