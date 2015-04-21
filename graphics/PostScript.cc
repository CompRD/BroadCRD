// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "graphics/PostScript.h"

std::ostream& ns_psplot::operator<<( std::ostream& o, const ns_psplot::text& t )
{    
  o << t.Col( ) << ns_psplot::ps_string( t.Str( ) ) << " newpath "
    << "/Times-Roman findfont " << t.FontSize( ) << " scalefont setfont ";
  if ( t.origin_ == ns_psplot::left ) 
    o << "0 0";
  else if ( t.origin_ == ns_psplot::center )
    o << ns_psplot::ps_string( t.Str( ) ) << " stringwidth pop 2 div neg 0";
  else if ( t.origin_ == ns_psplot::right )
    o << ns_psplot::ps_string( t.Str( ) ) << " stringwidth pop neg 0";
  o << " moveto show\n";
  return o;
}

std::ostream& ns_psplot::operator<<( std::ostream& o, const ns_psplot::segment& seg )
{
  float delta_x = seg.x2_ - seg.x1_;
  float delta_y = seg.y2_ - seg.y1_;

  return o << seg.col_
	   << ns_psplot::startpath( seg.x1_, seg.y1_ )
	   << ns_psplot::rlineto( delta_x, delta_y)
	   << " stroke\n";
}

// Need this (especially so gcc 4.2 like it) so in proper namespace
// Got rid of leading ns_psplot:: in front of "operator" below because of this.
namespace ns_psplot {
bool operator==(const ns_psplot::segment &seg1, const ns_psplot::segment &seg2)
{
  return ( seg1.X1() == seg2.X1() && seg1.Y1() == seg2.Y1() &&	   
	   seg1.X2() == seg2.X2() && seg1.Y2() == seg2.Y2() );
}
}
  

