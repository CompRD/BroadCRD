/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_DATA_PLOTTER_H
#define C_DATA_PLOTTER_H

#include "String.h"
#include "graphics/Whiteboard.h"

using namespace ns_whiteboard;

/**
 * class CDataPlotter
 *
 * It handles the data points for one set of data of a CContigPlotter. 
 * Scaling and shift can be performed on the two axes independently.
 */
class CDataPlotter {

public:

  CDataPlotter( color col = red, double width = 1.0 );
  
  ~CDataPlotter( );
  
  void Reserve( int size );
  
  void AddPoint( double xpt, double ypt );

  void SetScaling( double hscal, double vscal, double hshift, double vshift );

  void SetHorizontalScaling( double scaling ) { hscaling_ = scaling; }
  
  void SetVerticalScaling( double scaling ) { vscaling_ = scaling; }
  
  void SetHorizontalShift( double shift ) { hshift_ = shift; }
  
  void SetVerticalShift( double shift ) { vshift_ = shift; }
  
  void Sort( ) { sort( data_.begin( ), data_.end( ) ); }

  int Size( ) const { return data_.size( ); }
  
  double HorizontalScaling( ) const { return hscaling_; }
  
  double VeriticalScaling( ) const { return vscaling_; }

  double HorizontalShift( ) const { return hshift_; }

  double VeriticalShift( ) const { return vshift_; }
  
  // Min of data (it asserts if data_ is empty).
  double MinData( ) const;
  
  // Max of data (it asserts if data_ is empty).
  double MaxData( ) const;
  
  // Plot data as a line (it asserts if data_ is empty).
  void PlotAsLines( whiteboard &board, bool terminate = false ) const;
  
  // Plot data as points (it asserts if data_ is empty).
  void PlotAsPoints( whiteboard &board ) const;
  
  
private:

  // Create a pointer to a point.
  point *Point( double xpt, double ypt ) const;
 
  // Create a pointer to a terminal line (little vertical segment at point).
  line *TerminalLine( double xpt, double ypt ) const;
  
  // Create a pointer to a line.
  line *Line( double xpt1, double ypt1, double xpt2, double ypt2 ) const;
  
  
private:

  color color_;                       // color of line
  double width_;                      // width of line
  vec< pair<double, double> > data_;  // (x,y) coordinates of data points
  
  double hscaling_;  // horizontal scaling
  double vscaling_;  // vertical scaling
  double hshift_;    // horizontal shift
  double vshift_;    // vertical shift

  mutable vec< line* > lines_;    // lines to be plotted
  mutable vec< point* > points_;  // points to be plotted

};

#endif // C_DATA_PLOTTER_H
