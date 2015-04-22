/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// DisplayMatrix generates a png file consisting of a matrix of colored letters,
/// together with lines, which are defined by integer coordinates relative to
/// the matrix.  The upper left-hand corner is (0,0).  Uses CourierBold font.

#ifndef DISPLAY_MATRIX_H
#define DISPLAY_MATRIX_H

#include "CoreTools.h"
#include "graphics/Color.h"

class line {

     public:

     int x1, y1;
     int x2, y2;

     line( const int x1, const int y1, const int x2, const int y2 )
          : x1(x1), y1(y1), x2(x2), y2(y2) { }

};

void DisplayMatrix( const vec< vec< pair<char,color> > > letters,
     const vec<line> lines, const double fontsize, const String& outfile );
 
#endif
