/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_TIC_H
#define C_TIC_H

#include "math/Functions.h"
#include "String.h"
#include "graphics/Whiteboard.h"

/**
 * class CTic
 *
 * Used to decorate the axes of a CContigPlotter.
 */
class CTic {

public:

  CTic( double pos, const String &name );
    
  CTic( double pos, double size, const String &font, const String &name );
  
  ~CTic( );

  void SetScaling( double hscal, double vscal, double hshift, double vshift );
  
  void SetHorizontalScaling( double scaling ) { hscaling_ = scaling; }
  
  void SetVerticalScaling( double scaling ) { vscaling_ = scaling; }
  
  void SetHorizontalShift( double shift ) { hshift_ = shift; }
  
  void SetVerticalShift( double shift ) { vshift_ = shift; }

  void SetPos( double pos ) { pos_ = pos; }

  void SetSize( double size ) { size_ = size; }

  void SetFont( const String &font ) { font_ = font; }

  void SetName( const String &name ) { name_ = name; }
  
  double HorizontalScaling( ) const { return hscaling_; }
  
  double VeriticalScaling( ) const { return vscaling_; }

  double HorizontalShift( ) const { return hshift_; }

  double VeriticalShift( ) const { return vshift_; }

  double Pos( ) const { return pos_; }

  double Size( ) const { return size_; }

  const String &Font( ) const { return font_; }

  const String &Name( ) const { return name_; }
  
  // Plot tic on the x axis.
  void PlotXTic( ns_whiteboard::whiteboard &board ) const;
  
  // Plot tic on the y axis.
  void PlotYTic( ns_whiteboard::whiteboard &board ) const;
  
  
private:

  double pos_;    // position on the axis
  double size_;   // font size
  String font_;   // font name
  String name_;   // decoration (what is printed)

  double hscaling_;    // horizontal scaling
  double vscaling_;    // vertical scaling
  double hshift_;      // horizontal shift
  double vshift_;      // vertical shift

  mutable vec< ns_whiteboard::line* > lines_;   // lines to be plotted
  mutable vec< ns_whiteboard::text* > texts_;   // texts to be plotted
  
};

#endif // C_TIC_H
