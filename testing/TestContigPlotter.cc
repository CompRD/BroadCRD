/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "graphics/CContigPlotter.h"

/**
 * TestContigPlotter
 *
 * Save eps outfile on OUTFILE.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( OUTFILE );
  EndCommandArguments;
  
  // Contig length.
  int len = 354624;

  // Create a plotter.
  CContigPlotter plotter( len );
  
  // Generate two data sets (but wait before adding the data to the plotter).
  CDataPlotter data1( red );
  data1.Reserve( 4 );
  data1.AddPoint( 3000, 2.5 );
  data1.AddPoint( 25000, 5.000 );
  data1.AddPoint( 125000, 4.500 );
  data1.AddPoint( 300000, 6.500 );

  CDataPlotter data2( blue );
  data2.Reserve( 6 );
  data2.AddPoint( 30000, 8.500 );
  data2.AddPoint( 150000, 9.500 );
  data2.AddPoint( 185000, 9.000 );
  data2.AddPoint( 280000, 2.500 );
  data2.AddPoint( 300000, 1.500 );
  data2.AddPoint( 310000, 2.200 );
  
  // Plot a grid (gray, thickness = 0.3 ).
  float thickness = 0.3;

  for (int ii=1; ii<=10; ii++) {
    CDataPlotter data( gray, thickness );
    data.AddPoint( 0, ii );
    data.AddPoint( len, ii );
    plotter.AddData( data );
  }

  for (int ii=0; ii<8; ii++) {
    CDataPlotter data( gray, thickness );
    data.AddPoint( ii * 50000, 0 );
    data.AddPoint( ii * 50000, 10.0 );
    plotter.AddData( data );
  }

  // Adding data after the grid will make the lines appear above the grid.
  plotter.AddData( data1 );
  plotter.AddData( data2 );
  
  // Generate figure.
  plotter.GoFigure( OUTFILE.c_str( ) );
    
}
