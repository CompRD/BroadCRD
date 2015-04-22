/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "graphics/CContigPlotter.h"
#include "graphics/CDataPlotter.h"
#include "graphics/CTic.h"

/**
 * CContigPlotter
 * Constructor
 */
CContigPlotter::CContigPlotter( int length ) :
  length_ ( length )
{
  this->CleanUp( );
}

/**
 * CContigPlotter
 * AddData
 */
void CContigPlotter::AddData( const CDataPlotter &newdata )
{
  data_.push_back( newdata );
}

/**
 * CContigPlotter
 * AddHorizontalTics
 */
void CContigPlotter::AddHorizontalTics( const vec<CTic> &tics )
{
  xtics_ = tics;
}

/**
 * CContigPlotter
 * AddVerticalTics
 */
void CContigPlotter::AddVerticalTics( const vec<CTic> &tics )
{
  ytics_ = tics;
}

/**
 * CContigPlotter
 * GoFigure
 */
void CContigPlotter::GoFigure( const String &outfile, bool aspoints )
{
  // These define the actual dimensions of the figure (in pixels).
  double width = 800;
  double height = 400;
  double border = 0;
  
  // Leave some empty space between the plot and the edge of the figure.
  double slack_width = 100;
  double slack_height = 50;
  
  // Raw width, height (with some override if data are on a flat line).
  double raw_min = data_[0].MinData( );
  double raw_max = data_[0].MaxData( );
  for (int ii=1; ii<(int)data_.size( ); ii++) {
    raw_min = Min( raw_min, data_[ii].MinData( ) );
    raw_max = Max( raw_max, data_[ii].MaxData( ) );
  }
  if ( raw_min == raw_max ) {
    raw_min = raw_min - 10.0;
    raw_max = raw_max + 10.0;
  }
  
  // If there are ytics, generate a y axis.
  this->GenerateYAxis( raw_min, raw_max );

  // Set scalings and shifts.
  double hscaling = ( width - 2 * slack_width ) / length_;
  double vscaling = ( height - 2 * slack_height ) / ( raw_max - raw_min );
  double hshift = slack_width;
  double vshift = Max( - raw_min + slack_height, slack_height );
  
  // Generate a whiteboard (plot x-axis last).
  whiteboard board;
  
  for (int ii=0; ii<(int)xtics_.size( ); ii++) {
    xtics_[ii].SetScaling( hscaling,  vscaling, hshift, vshift );
    xtics_[ii].PlotXTic( board );
  }
  
  for (int ii=0; ii<(int)ytics_.size( ); ii++) {
    ytics_[ii].SetScaling( hscaling,  vscaling, hshift, vshift );
    ytics_[ii].PlotYTic( board );
  }
  
  for (int ii=0; ii<(int)data_.size( ); ii++) {
    data_[ii].SetScaling( hscaling,  vscaling, hshift, vshift );
    data_[ii].PlotAsLines( board );
  }
  data_[0].SetScaling( hscaling,  vscaling, hshift, vshift );
  data_[0].PlotAsLines( board, true );
  
  // Open stream and print.
  ofstream out( outfile.c_str( ) );
  ps_display display( out, width, height, 0 );
  board.DisplayOn( &display );
  out.close( );
}

/**
 * CContigPlotter
 * CleanUp
 */
void CContigPlotter::CleanUp( )
{
  // Clean from old runs.
  data_.clear( );
  xtics_.clear( );
  ytics_.clear( );
  
  // The first entry in data_ always represents the contig.
  double thickness = 2.0;
  CDataPlotter contig( black, thickness );
  contig.AddPoint( 0, 0 );
  contig.AddPoint( length_, 0 );
  this->AddData( contig );
}

/**
 * CContigPlotter
 * GenerateYAxis
 */
void CContigPlotter::GenerateYAxis( double min, double max )
{
  if ( ytics_.size( ) < 1 )
    return;
  
  CDataPlotter yaxis( gray );
  yaxis.AddPoint( 0, min );
  yaxis.AddPoint( 0, max );
  this->AddData( yaxis );
}
