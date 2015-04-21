/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_CONTIG_PLOTTER_H
#define C_CONTIG_PLOTTER_H

#include "Vec.h"
#include "graphics/CDataPlotter.h"
#include "graphics/CTic.h"

/**
 * clss CContigPlotter
 *
 * Plot features on a contig. Features are data points expressed in a
 * coordinate system in which the contig is represented by the interval
 * [0, contig_length). Output is generated as an eps file.
 */
class CContigPlotter {

public:
  
  CContigPlotter( int length );
  
  // Add a new set of features.
  void AddData( const CDataPlotter &newdata );
  
  // Add horizontal tics.
  void AddHorizontalTics( const vec<CTic> &tics );

  // Add vertical tics.
  void AddVerticalTics( const vec<CTic> &tics );

  // Generate the figure.
  void GoFigure( const String &outfile, bool aspoints = false );
  
  
private:
  
  // Clean up from old runs and create data representation for contig.
  void CleanUp( );
  
  // If ytics_ is not empty, generate y axis.
  void GenerateYAxis( double min, double max );
  
  
private:
  
  int length_;              // length of the "contig"
  vec<CDataPlotter> data_;  // sets of contig (data_[0]) and data points
  vec<CTic> xtics_;         // decoration of x axis (optional)
  vec<CTic> ytics_;         // decoration of y axis (optional)

};

#endif // C_CONTIG_PLOTTER_H
