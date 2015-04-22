/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "math/Functions.h"
#include "String.h"
#include "graphics/CDataPlotter.h"
#include "graphics/Whiteboard.h"

/**
 * CDataPlotter
 * Constructor
 */
CDataPlotter::CDataPlotter( color col, double width ) :
  color_ ( col ),
  width_ ( width ),
  hscaling_ ( 1.0 ),
  vscaling_ ( 1.0 ),
  hshift_ ( 0 ),
  vshift_ ( 0 )
{ }

/**
 * CDataPlotter
 * Destructor
 */
CDataPlotter::~CDataPlotter( )
{
  for (int ii=0; ii<(int)lines_.size( ); ii++)
    delete ( lines_[ii] );
  for (int ii=0; ii<(int)points_.size( ); ii++)
    delete ( points_[ii] );
}

/**
 * CDataPlotter
 * Reserve
 */
void CDataPlotter::Reserve( int size )
{
  data_.reserve( size );
}

/**
 * CDataPlotter
 * AddPoint
 */
void CDataPlotter::AddPoint( double xpt, double ypt )
{
  data_.push_back( make_pair( xpt, ypt ) );
}

/**
 * CDataPlotter
 * SetScaling
 */
void CDataPlotter::SetScaling( double hscal,
			       double vscal,
			       double hshift,
			       double vshift )
{
  hscaling_ = hscal;
  vscaling_ = vscal;
  hshift_ = hshift;
  vshift_ = vshift;
}

/**
 * CDataPlotter
 * MinData
 */
double CDataPlotter::MinData( ) const
{
  if ( data_.size( ) < 1 )
    return 0;
    
  double min = data_[0].second;
  for (int ii=1; ii<(int)data_.size( ); ii++)
    min = Min( min, data_[ii].second );
  return min;
}

/**
 * CDataPlotter
 * MaxData
 */
double CDataPlotter::MaxData( ) const
{
  if ( data_.size( ) < 1 )
    return 0;

  double max = data_[0].second;
  for (int ii=1; ii<(int)data_.size( ); ii++)
    max = Max( max, data_[ii].second );
  return max;
}

/**
 * CDataPlotter
 * PlotAsLines
 *
 * If terminate=True it will tag with small segments first and last points.
 */
void CDataPlotter::PlotAsLines( whiteboard &board, bool terminate ) const
{
  if ( data_.size( ) < 1 )
    return;
  
  if ( terminate ) {
    board.Add(this->TerminalLine( data_[0].first, data_[0].second ));
    board.Add(this->TerminalLine( data_.back( ).first, data_.back( ).second ));
  }

  board.Add( this->Point( data_[0].first, data_[0].second ) );
  for (int ii=1; ii<(int)data_.size( ); ii++) {
    double x1 = data_[ii-1].first;
    double y1 = data_[ii-1].second;
    double x2 = data_[ii].first;
    double y2 = data_[ii].second;
    board.Add( this->Line( x1, y1, x2, y2 ) );
  }
}

/**
 * CDataPlotter
 * PlotAsPoints
 */
void CDataPlotter::PlotAsPoints( whiteboard &board ) const
{
  if ( data_.size( ) < 1 )
    return;

  for (int ii=0; ii<(int)data_.size( ); ii++)
    board.Add( this->Point( data_[ii].first, data_[ii].second ) );
}

/**
 * CDataPlotter
 * Point
 */
point *CDataPlotter::Point( double xpt, double ypt ) const
{
  double x1 = hshift_ + xpt * hscaling_;
  double y1 = vshift_ + ypt * vscaling_;
  xy_coords p1 = make_pair( x1, y1 );
  
  point *newpoint = new point( p1, width_, color_ );
  points_.push_back( newpoint );
  return newpoint;
}

/**
 * CDataPlotter
 * TerminalLine
 */
line *CDataPlotter::TerminalLine( double xpt, double ypt ) const
{
  double thickness = 4.0;

  double xcenter = hshift_ + xpt * hscaling_;
  double ycenter = vshift_ + ypt * vscaling_;
  xy_coords p1 = make_pair( xcenter, ycenter - thickness );
  xy_coords p2 = make_pair( xcenter, ycenter + thickness );
  
  line *newline = new line( p1, p2 , width_, color_ );
  lines_.push_back( newline );
  return newline;
}

/**
 * CDataPlotter
 * Line
 */
line *CDataPlotter::Line( double xpt1, double ypt1,
			  double xpt2, double ypt2 ) const
{
  double x1 = hshift_ + xpt1 * hscaling_;
  double y1 = vshift_ + ypt1 * vscaling_;
  xy_coords p1 = make_pair( x1, y1 );
  
  double x2 = hshift_ + xpt2 * hscaling_;
  double y2 = vshift_ + ypt2 * vscaling_;
  xy_coords p2 = make_pair( x2, y2 );
  
  line *newline = new line( p1, p2 , width_, color_ );
  lines_.push_back( newline );
  return newline;
}

