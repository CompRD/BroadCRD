// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "assembly/display/PsRead.h"

#include "graphics/Color.h"

void
PsRead::PrintHeader( ostream &out )
{
  out << "/printread {\n";
  out << "  /repetitive exch def\n";
  out << "  /filled exch def\n";
  out << "  /tailXdiff exch def\n";
  out << "  /originY exch def\n";
  out << "  /originX exch def\n";
  out << "  /readHeight " << s_readHeight << " def\n";
  out << "  \n";
  out << "  newpath\n";
  out << "    originX originY moveto\n";
  out << "    0 readHeight rlineto\n";
  out << "    tailXdiff readHeight neg rlineto\n";
  out << "  filled 1 eq\n";
  out << "  { gsave\n";
  out << "    fill\n";
  out << "    grestore }\n";
  out << "  if\n";
  out << "  stroke\n";
  out << "  newpath\n";
  out << "    originX originY moveto\n";
  out << "    tailXdiff 0 rlineto\n";
  out << "  repetitive 1 eq\n";
  out << "  { gsave\n";
  out << "    0.2 setgray\n";
  out << "    stroke\n";
  out << "    grestore }\n";
  out << "  { stroke }\n";
  out << "  ifelse\n";
  out << "} def\n";
  out << "\n";
  out << "/printreadname {\n";
  out << "  /tailXdiff exch def\n";
  out << "  /originY exch def\n";
  out << "  /originX exch def\n";
  out << "  /name exch def\n";
  out << "  \n";
  out << "  originX originY moveto\n";
  out << "  tailXdiff name stringwidth pop sub 2 div 0 rmoveto\n";
  out << "  name show\n";
  out << "} def\n";
  out << "\n";
}

void
PsRead::Print( ostream &out, const bool filled ) const
{
  out << m_originX << " " 
      << m_originY << " " 
      << m_tailXDiff << " "
      << ( filled ? "1" : "0" ) << " "
      << ( m_repetitive ? "1" : "0" ) << " "
      << "printread\n";
}

void
PsRead::PrintId( ostream &out ) const
{
  out << "(" << m_read.GetId() << ") "
      << m_originX << " " 
      << m_originY << " " 
      << m_tailXDiff << " printreadname\n";
}

void
PsRead::PrintName( ostream &out ) const
{
  out << "(" << m_read.GetName() << ") "
      << m_originX << " " 
      << m_originY + 1.0 << " " 
      << m_tailXDiff << " printreadname\n";
}

