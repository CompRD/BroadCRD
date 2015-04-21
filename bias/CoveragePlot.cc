/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// CoveragePlot.  Plot coverage from alignments of reads to reference.
          
#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "graphics/BasicGraphics.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"

void CoveragePlot( const vecbasevector& reads, const vec<look_align>& aligns,
     const vecbasevector& ref, const int TIG, int START, int STOP,
     const int WINDOW, const int GRANULARITY, const double COVERAGE_DIVIDER,
     const String& OUT )
{
     // Hardwired graphics parameters.

     const double LINEWIDTH = 1.0;
     const double SCALE = 1.0;
     const double XTRANS = 200;
     const double YTRANS = 50;
     const Bool SHOW = True;

     // Compute read start points.

     int N = ref[TIG].size( );
     vec< vec<int> > starts(2);
     for ( int i = 0; i < 2; i++ )
          starts[i].resize( N, 0 );
     for ( int j = 0; j < aligns.isize( ); j++ )
     {    const look_align& la = aligns[j];
          int id = la.query_id, id2 = la.target_id;
          if ( id2 != TIG ) continue;
          align a = la.a;
          const basevector& rd1 = reads[la.query_id];
          if (la.rc1) a.ReverseThis( rd1.size( ), ref[id2].size( ) );
          ForceAssertGe( a.pos2( ), 0 );
          ForceAssertLe( static_cast<unsigned>(a.Pos2()), ref[id2].size( ) );
          ++starts[ la.rc1 ? 1 : 0 ][ a.pos2( ) ];    }

     // Find number of start points in windows.

     vec<int> count( N, 0 );
     for ( int i = 0; i < WINDOW; i++ )
          count[0] += starts[0][i] + starts[1][N-i-1];
     for ( int i = 1; i < WINDOW/2; i++ )
          count[i] = count[0];
     for ( int j = 0; j < N; j++ )
     {    if ( j >= WINDOW/2 && j < N - WINDOW/2 )
          {    count[j] = count[j-1];
               count[j] += starts[0][j+WINDOW/2] + starts[1][N-(j+WINDOW/2)-1];
               count[j] 
                    -= starts[0][j-WINDOW/2] + starts[1][N-(j-WINDOW/2)-1];    }    }

     // Set up to plot.

     if ( STOP < 0 ) STOP = N;
     Bool to_png = OUT.Contains( ".png", -1 );
     ForceAssert(to_png);
     vec<graphics_primitive> points, lines;
     lines.push_back( SetLineWidth(LINEWIDTH) );

     // Plot.

     double xm = 0.0, ym = 0.0;
     for ( int i = START; i < STOP; i += WINDOW/GRANULARITY )
     {    double x = i, y = double(count[i])/COVERAGE_DIVIDER;
          if ( i > START ) lines.push_back( Segment( xm, ym, x, y ) );
          points.push_back( Point( x, y, 1.0 ) );
          xm = x, ym = y;    }

     double x = MinX(points), X = MaxX(points);
     double y = MinY(points), Y = MaxY(points);
     vec<graphics_primitive> x_axis = AxisX( x, X );
     vec<graphics_primitive> y_axis = AxisY( x, X, y, Y );
     lines.append(y_axis);

     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(x_axis),   heights.push_back(30);
     stack.push_back(lines),    heights.push_back(200);

     String outhead = OUT.Before( ".png" );
     double width = 400;
     {    Ofstream( out, outhead + ".ps" );
          VerticalDisplay( stack, heights, width, out, SCALE, XTRANS,
               YTRANS, SHOW );    }
     if (to_png)
     {    String command
               = "pstopnm -portrait -xmax 8000 -ymax 8000 -stdout "
               + outhead + ".ps | pnmcrop | "
               + "pnmpad -white -left=25 -right=25 -top=25 -bottom=25 | "
               + "pnmtopng > " + outhead + ".png";
          int status = System(command);
          if ( status == 0 ) Remove( outhead + ".ps" );
          else cout << "failed to run:\n" << command << "\n";     }    }
