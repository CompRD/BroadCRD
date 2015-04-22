///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SVG_COLORS_H
#define SVG_COLORS_H

// Compute distance between two RGBPC colors according to the CIEDE2000 formula.
// Return value between 0 and somewhat over 107.  Figuring out the exact max
// would require testing 10^12 possibilities (which is doable).

double ColorDist( const vec<int>& rgbpc1, const vec<int>& rgbpc2 );

#endif
