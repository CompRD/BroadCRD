/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "graphics/Color.h"
#include "graphics/DisplayMatrix.h"

void DisplayMatrix( const vec< vec< pair<char,color> > > letters,
     const vec<line> lines, const double fontsize, const String& outfile )
{
     ForceAssert( outfile.Contains( ".png", -1 ) );
     String outhead = outfile.Before( ".png" );
     Ofstream( out, outhead + ".ps" );

     int nrows = letters.size( );

     const double fwmult = 0.7;
     const double fwmultsub = 0.05;
     double fontwidth = fontsize * fwmult;

     out << "200 50 translate\n";
     out << "1 1 scale\n";
     out << "/CourierBold findfont " << fontsize << " scalefont setfont\n";

     for ( int row = 0; row < nrows; row++ )
     {    for ( int i = 0; i < letters[row].isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < letters[row].isize( ); j++ )
                    if ( letters[row][j].second != letters[row][i].second ) break;
               out << letters[row][i].second << " ";

               for ( int k = i; k < j; k++ )
               {    out << "newpath " << k * fontwidth << " "
                         << (nrows-row) * fontsize << " moveto ";
                    out << "(" << letters[row][k].first << ")";
                    out << " show\n";    }

               i = j - 1;    }    }

     out << "0 0 0 setrgbcolor\n" << "0.25 setlinewidth\n";
     for ( int i = 0; i < lines.isize( ); i++ )
     {    const line& L = lines[i];
          out << "newpath " << fontwidth * ( L.x1 - fwmultsub ) << " "
               << ( nrows - L.y1 + 1 - 0.25 ) * fontsize << " moveto "
               << fontwidth * ( L.x2 - fwmultsub ) << " "
               << ( nrows - L.y2 + 1 - 0.25 ) * fontsize << " lineto stroke\n";    }

     out << "showpage" << endl;

     Bool to_png = True;
     const double POSTSCALE = 0.3;

     if (to_png)
     {    String command1 = "ps2epsi " + outhead + ".ps " + outhead + ".eps";
          int status1 = System(command1);
          if ( status1 == 0 ) Remove( outhead + ".ps" );
          else
          {    cout << "failed to run:\n" << command1 << "\n";
               cout << "Abort.\n";
               TracebackThisProcess( );    }
          String command2 = "pstopnm -portrait -xmax 8000 -ymax 8000 -stdout "
               + outhead + ".eps | pnmcrop | "
               + "pnmpad -white -left=25 -right=25 -top=25 -bottom=25 | ";
          if ( POSTSCALE != 1.0 )
               command2 += "pamscale " + ToString(POSTSCALE,3) + " | ";
          command2 += "pnmtopng > " + outhead + ".png";
          int status2 = System(command2);
          if ( status2 == 0 ) Remove( outhead + ".ps" );
          else
          {    cout << "failed to run:\n" << command2 << "\n";
               cout << "Abort.\n";
               TracebackThisProcess( );    }    }    }
