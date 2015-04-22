///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__MONOMER_UTILS_TALENS__H
#define BTL__MONOMER_UTILS_TALENS__H

#include "String.h"
#include "Vec.h"
#include "lookup/LookAlign.h"

/**
 * ShortName
 *
 * It turns monomerHD into HD, and terminHD into HD.t
 */
String ShortName( const String &name );

/**
 * MonomersOverlap
 *
 * Decide if the two aligns overlap on target.
 */
bool MonomersOverlap( const look_align_plus &al1,
		      const look_align_plus &al2 );

/**
 * AlignsToChains
 *
 * Build a vector of aligns (chains).
 */
void AlignsToChains( const vec<look_align_plus> &hits,
		     const vec<bool> &is_termin,
		     vec<look_align_plus> &chains );

/**
 * PrintChains
 */
void PrintChains( const vec<look_align_plus> &chains,
		  const vec<bool> &is_termin,
		  const vecString &parts_ids,
		  ostream &log );

/**
 * PrintError
 */
void PrintError( const vec< pair<int,int> > &to_putative,
		 const vecString &parts_ids,
		 const String &descrip,
		 const int qq,
		 const int pos,
		 ostream &log );

#endif
