/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef TILING_PRINTER_H
#define TILING_PRINTER_H

#include "PackAlign.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "Vec.h"
#include "tiled/PaddedSeq.h"
#include "tiled/Tiling.h"
#include "tiled/TilingAnalyzer.h"

/**
 * class tiling_printer
 *
 * It manages a tiling for printing. Printing options (a default value is
 * provided for all options in ResetOptions) are:
 *
 * width_: defines the width of the printed chunks. There will be width
 *  bases and/or quals per line (plus the space at the beginning reserved
 *  to the name of the read).
 *
 * interval [to_print_begin_, to_print_end_): defines the window on the master
 *  sequence to be printed. The default value is (-1, -1), which means print
 *  all. Notice that begin and end refer to the padded master, i.e. lengths
 *  are counted including pads.
 *
 * print_wings_: if true, and if a read extends the master on the left or on
 *  the right, then print also extending bases. If false, then print only up
 *  to (and not over) to_print_begin_, to_print_end_.
 *
 * print_quals_: print also quality scores. They will be printed in a
 *  short hand format, i.e. you will see only the first digit of the
 *  phred quality score.
 *
 * print_discrep_: if true, print also an extra line on top of each chunk,
 *  tagging all columns for which there is at least one discrepancy.
 */
class tiling_printer {
  
public:
  
  tiling_printer( tiling_analyzer *analyzer );

  void ResetOptions( );

  void SetWidth ( int width );

  void SetWindowToPrint( int begin, int end );

  void SetPrintWings( bool print_wings );
 
  void SetPrintQuals( bool print_quals );

  void SetPrintDiscrep( bool print_discrep );
  
  void ToStream( ostream &out, vec<int> *read_pos = 0 ) const;

  void ToStream454( ostream &out ) const;
  
  
private:

  void ActualWindow( pair<int, int> &win ) const;
  
  // begin, end in padded coordinates.
  void PrintCoordinates( ostream &out, int begin, int end ) const;
  
  // Generate the line with tags for discrepancies.
  vec<char> DiscrepTag( const vec< vec<char> > &rbases, bool tagKC ) const;
  
  
private:
  
  mutable tiling_analyzer *analyzer_;

  int width_;
  int to_print_begin_;
  int to_print_end_;
  bool print_wings_;
  bool print_quals_;
  bool print_discrep_;

};

#endif
