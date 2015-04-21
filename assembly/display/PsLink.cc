// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology


#include "assembly/display/PsLink.h"
#include "graphics/Color.h"

void
PsLink::PrintHeader( ostream &out )
{
  out << "/drawlink {\n";
  out << "  /linkHeight exch def\n";
  out << "  /linkHalfDist exch def\n";
  out << "  /linkStartY exch def\n";
  out << "  /originY exch def\n";
  out << "  /originX exch def\n";
  out << "\n";
  out << "  newpath\n";
  out << "    originX originY moveto\n";
  out << "    linkHalfDist abs 10 gt\n";
  out << "    { matrix currentmatrix\n";
  out << "      originX linkHalfDist add linkStartY translate\n";
  out << "      1 linkHeight linkHalfDist abs div scale\n";
  out << "      linkHalfDist 0 ge\n";
  out << "      { 0 0 linkHalfDist 180 90 arcn }\n";
  out << "      { 0 0 linkHalfDist neg 0 90 arc }\n";
  out << "      ifelse\n";
  out << "      setmatrix }\n";
  out << "    { originX linkStartY lineto\n";
  out << "      0 linkHeight 0 linkHeight linkHalfDist linkHeight rcurveto }\n";
  out << "    ifelse\n";
  out << "  stroke\n";
  out << "} def\n";
  out << "\n";
  out << "/drawoffsupermark {\n";
  out << "  /linkHeight exch def\n";
  out << "  /linkHalfDist exch def\n";
  out << "  /linkStartY exch def\n";
  out << "  /originY exch def\n";
  out << "  /originX exch def\n";
  out << "  /superId exch def\n";
  out << "\n";
  out << "/Helvetica findfont " << 2 * PsRead::s_readHeight << " scalefont setfont\n";
  out << "  originX linkHalfDist add originY linkHeight add moveto\n";
  out << "  linkHalfDist 0 lt\n";
  out << "  { superId stringwidth pop neg 0 rmoveto }\n";
  out << "  if\n";
  out << "  superId show\n";
  out << "} def\n";
  out << "\n";
}


void
PsLink::Print( ostream& out ) const
{
  if ( ! mp_read1 && ! mp_read2 )
    return;

  const float baseColorVal = 0.9;
  
  // Figure out link color and link height.
  color linkColor( red );
  float linkHeight;
  float linkHalfDist;

  bool isSuperStretched = false;
  bool isOrphaned = ( mp_read2 == 0 );

  if ( isOrphaned )
  {
    linkColor = color( 1.0, 0.7, 0.0 );
    linkHeight = m_expectedSize / 2.0;
    linkHalfDist = m_expectedSize / 2.0;
    if ( mp_read1->IsReverse() )
      linkHalfDist = -linkHalfDist;
  }
  else
  {
    bool isLogical = ( mp_read1->IsForward() != mp_read2->IsForward() );
    
    float observedDist = mp_read2->GetOriginX() - mp_read1->GetOriginX();

    float observedSize = observedDist;
    if ( isLogical && mp_read1->IsReverse() )
      observedSize = -observedSize;

    isSuperStretched = ( fabs( m_expectedSize - observedSize ) > 100000.0 );

    if ( ! isLogical )
    {
      linkColor = color( 1.0, 0.6, 0.6 );
      linkHeight = m_expectedSize / 2.0;
      linkHalfDist = observedDist / 2.0;
    }
    else
    {
      const float maxStretch = 10.0;
      const float minStretch = -10.0;
          
      float cappedStretch = ( observedSize - (float)m_expectedSize ) / (float)m_expectedStdDev;
      cappedStretch = min( cappedStretch, maxStretch );
      cappedStretch = max( cappedStretch, minStretch );
      
      float keyColorVal = baseColorVal - baseColorVal * fabs( cappedStretch ) / maxStretch;
      
      if ( cappedStretch > 0.0 )
        linkColor = color( keyColorVal, baseColorVal, keyColorVal );
      else
        linkColor = color( keyColorVal, keyColorVal, baseColorVal );

      linkHeight = m_expectedSize / 2.0;
      linkHalfDist = observedDist / 2.0;
    }

    if ( isSuperStretched  ) {
      // raise them up over the fray
      linkHeight = 50000;
      // fade them
      float r = linkColor.R(), g = linkColor.G(), b = linkColor.B();
      r = r + (1.-r)/2.;
      g = g + (1.-g)/2.;
      b = b + (1.-b)/2.;
      linkColor = color( r, g, b );
    }

  }
   
  // Draw reads.

  out << "% link with importance " << m_importance << "\n";
  out << linkColor << "\n";

  float maxReadY = 0;
  if ( mp_read1 )
  {
    maxReadY = max( maxReadY, mp_read1->GetOriginY() );
    mp_read1->PrintFilled( out );
  }
      
  if ( mp_read2 )
  {
    maxReadY = max( maxReadY, mp_read2->GetOriginY() );
    mp_read2->PrintFilled( out );
  }

  if ( fabs(linkHalfDist) < 1.0 )
    if ( linkHalfDist > 0.0 )
      linkHalfDist = 1.0;
    else
      linkHalfDist = -1.0;
  // Draw link.
  
  if ( mp_read1 )
  {
    if ( isOrphaned )
      out << "[250 750] 0 setdash" << "\n";
    else if ( isSuperStretched )
      out << "[1000 3000] 0 setdash" << "\n";

    out << mp_read1->GetOriginX() << " " 
        << mp_read1->GetOriginY() << " "
        << maxReadY << " "
        << linkHalfDist << " "
        << linkHeight << " "
        << "drawlink\n";

    if ( isOrphaned || isSuperStretched )
      out << "[] 0 setdash" << "\n";
  }

  if ( mp_read2 )
  {
    linkHalfDist = -linkHalfDist;

    if ( isSuperStretched )
      out << "[1000 3000] 0 setdash" << "\n";

    out << mp_read2->GetOriginX() << " " 
        << mp_read2->GetOriginY() << " "
        << maxReadY << " "
        << linkHalfDist << " "
        << linkHeight << " "
        << "drawlink\n";

    if ( isSuperStretched )
      out << "[] 0 setdash" << "\n";
  }
  else
  {
    out << black << "\n";
    out << "(" << m_super2.GetId() << ") "
        << mp_read1->GetOriginX() << " " 
        << mp_read1->GetOriginY() << " "
        << maxReadY << " "
        << linkHalfDist << " "
        << linkHeight << " "
        << "drawoffsupermark\n";
  }
}
