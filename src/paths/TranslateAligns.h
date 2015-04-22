///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef TRANSLATE_ALIGNS_H
#define TRANSLATE_ALIGNS_H

#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"

/**
 * TranslateAlignsToSupers
 *
 * Load alignments and supers, and translate alignments from contig
 * coordinates to super coordinates.
 */
void TranslateAlignsToSupers( const vec<superb> &supers,
			      vec<alignlet> &aligns );

#endif
