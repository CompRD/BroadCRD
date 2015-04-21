///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef TEXPS_H
#define TEXPS_H

// Go through the DISCOVAR lines, generating a data structure texps, that
// we will use for mapping the lines to the BNG maps.
//
// This code also prints out a human-readable representation of the 
// translation of each line to "map space".
//
// Example:
//
// 40630:
// {{11426,2888,7665,382,540,2717,1527,375,2627,723}}
// 
// This is the simplest case: the line is unambiguously expanded as a list
// of restrictions-site-free intervals.  Note that the first and last invervals
// are incomplete and thus special.
//
// (Identical documentation in .cc.)

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

void Texps( const HyperBasevectorX& hb, const vec<vec<vec<vec<int>>>>& lines,
     const String& cut, const int min_line, const int max_ignored_indel,
     vec<vec<vec<int>>>& texps, const Bool verbose, ostream& out );

#endif
