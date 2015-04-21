///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__LOADER_TALENS__H
#define BTL__LOADER_TALENS__H

#include "Basevector.h"
#include "String.h"
#include "Vec.h"

/**
 * LoadParts
 *
 * Merge parts into a single fastb file, and generate a parallel
 * (finished grade) qualb file (not really a loader module). It
 * returns the name of the parts vecbvec (or "" on error).
 */
String LoadParts( const String &PARTS_DIR,
		  ostream *log = 0 );

/**
 * LoadProducts
 *
 * Load names of products from a given dir, by looking at all the
 * matching sets of .ab1 files. This now means that for each <name>
 * there must be files <name>-{F1,R1}.ab1.
 *
 * NB: if R2 is true, look for <name>-R2.ab1 (not for <name>-R1.ab1).
 */
void LoadProducts( const String &AB1_DIR,
		   const bool R2,
		   vec<String> &prods,
		   ostream &log );

/**
 * LoadPrimers
 *
 * Parse input fastsavector file, and load as fastb.
 */
void LoadPrimers( const String &in_fastavec_file,
		  vecbvec &primers );

/**
 * LoadReferences
 *
 * Parse input references file (if != "skip"), and generate references
 * in sync with the given product names (it returns false on error).
 */
bool LoadReferences( const String &in_references_file,
		     const vec<String> &prod_names,
		     vec<String> &references,
		     ostream *log = 0 );

#endif
