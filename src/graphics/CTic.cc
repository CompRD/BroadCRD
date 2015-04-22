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
#include "graphics/CTic.h"

using namespace ns_whiteboard;

/**
 * CTic
 * Constructor
 */
CTic::CTic( double pos, const String &name ) :
  pos_ ( pos ),
  name_ ( name )
{
  size_ = 8;
  font_ = "Arial";
  this->SetScaling( 1.0, 1.0, 0.0, 0.0 );
}

/**
 * CTic
 * Constructor
 */
CTic::CTic( double pos, double size, const String &font, const String &name ) :
  pos_ ( pos ),
  size_ ( size ),
  font_ ( font ),
  name_ ( name )
{
  this->SetScaling( 1.0, 1.0, 0.0, 0.0 );
}

/**
 * CTic
 * Destructor
 */
CTic::~CTic( )
{
  for (int ii=0; ii<(int)lines_.size( ); ii++)
    delete ( lines_[ii] );
  for (int ii=0; ii<(int)texts_.size( ); ii++)
    delete( texts_[ii] );
}

/**
 * CTic
 * SetScaling
 */
void CTic::SetScaling( double hscal,
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
 * CTic
 * PlotXTic
 */
void CTic::PlotXTic( whiteboard &board ) const
{
  double thickness = 4.0;
  double gap = 8;

  double xcenter = hshift_ + pos_ * hscaling_;
  double ycenter = vshift_;
  xy_coords p1 = make_pair( xcenter, ycenter - thickness );
  xy_coords p2 = make_pair( xcenter, ycenter + thickness );
  
  line *newline = new line( p1, p2 , 1.0, gray );
  lines_.push_back( newline );
  board.Add( newline );
  
  xy_coords p3 = make_pair( xcenter, ycenter - gap );
  
  text *newtext = new text( p3, name_, black, size_, font_, 270.0 );
  newtext->SetHorizAlign( align_left );
  newtext->SetVertAlign( align_middle );
  texts_.push_back( newtext );
  board.Add( newtext );
}

/**
 * CTic
 * PlotYTic
 */
void CTic::PlotYTic( whiteboard &board ) const
{
  double thickness = 4.0;
  double gap = 8;

  double xcenter = hshift_;
  double ycenter = vshift_ + pos_ * vscaling_;
  xy_coords p1 = make_pair( xcenter - thickness, ycenter );
  xy_coords p2 = make_pair( xcenter + thickness, ycenter );
  
  line *newline = new line( p1, p2 , 1.0, gray );
  lines_.push_back( newline );
  board.Add( newline );
  
  xy_coords p3 = make_pair( xcenter - gap, ycenter );
  text *newtext = new text( p3, name_, black, size_, font_ );
  newtext->SetHorizAlign( align_right );
  newtext->SetVertAlign( align_middle );
  texts_.push_back( newtext );
  board.Add( newtext );
}

