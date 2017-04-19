///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Generate a dotplot as a png or ps or pdf file OUT, depending upon its extension.  
// Input: a file with entries:
// x1 y1
// x2 y2
// ...
// Arguments PR, PG, PB control point color.
//
// MAX_Y: cap Y values by MAX_Y and use it for the highest value on the y-axis
//
// COLOR_BY_POINT=True: instead the input file has entries 
// x1 y1 r1 g1 b1
// ...
// specifying the color of each point.
//
// If CONNECT=True, connect points with segments of width LINEWIDTH.
//
// If VERT_FILE is non-null, it reads the file of x values at which points
// a blue vertical line is drawn.  This can be used, for example, to denote
// exon boundaries in a coverage plot.
//
// If SEGMENT_FILE is non-null, it reads the file of x1 and x2 values at
// which points a green then magenta vertical lines are drawn.  This has
// similar applications to VERT_FILE.  Will always draw lines if there is 
// just one segment

#include "graphics/BasicGraphics.h"
#include "MainTools.h"
#include "math/Functions.h"
#include <map>

const double undefined = 1000000000;

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     CommandArgument_Double_OrDefault(POINTSIZE, 1.0);
     CommandArgument_Double_OrDefault(SCALE, 1.0);
     CommandArgument_Double_OrDefault(POSTSCALE, 1.0);
     CommandArgument_Double_OrDefault(PR, 0);
     CommandArgument_Double_OrDefault(PG, 0);
     CommandArgument_Double_OrDefault(PB, 0);
     CommandArgument_String_OrDefault(TITLE, "");
     CommandArgument_Double_OrDefault(XTRANS, 200);
     CommandArgument_Double_OrDefault(YTRANS, 50);
     CommandArgument_Bool_OrDefault(SHOW, True);
     CommandArgument_Double_OrDefault(MIN_Y, undefined);
     CommandArgument_Double_OrDefault(MAX_Y, undefined);
     CommandArgument_UnsignedInt_OrDefault(TITLE_FONTSIZE, 15);
     CommandArgument_Bool_OrDefault(COLOR_BY_POINT, False);
     CommandArgument_Bool_OrDefault(CONNECT, False);
     CommandArgument_Double_OrDefault(LINEWIDTH, 1.0);
     CommandArgument_Bool_OrDefault(POINTS_LINES_TOGETHER, True);
     CommandArgument_Int_OrDefault(X_AXIS_OFFSET, 30);
     CommandArgument_Double_OrDefault(X_AXIS_EXTEND, 0.05);
     CommandArgument_Double_OrDefault(Y_AXIS_EXTEND, 0.05);
     CommandArgument_Bool_OrDefault(X_AXIS_NAKED, False);
     CommandArgument_Bool_OrDefault(Y_AXIS_NAKED, False);
     CommandArgument_Double_OrDefault(XMIN, undefined);
     CommandArgument_Double_OrDefault(XMAX, undefined);
     CommandArgument_String_OrDefault(VERT_FILE, "");
     CommandArgument_String_OrDefault(SEGMENT_FILE, "");
     CommandArgument_Int_OrDefault(YHEIGHT, 200);
     CommandArgument_Bool_OrDefault(Y_MINOR_TICS, True);
     CommandArgument_Bool_OrDefault(PIPEFAIL, True);
     EndCommandArguments;

     CheckValidGraphicsSuffix( OUT, SHOW );

     vec<graphics_primitive> points, lines;
     lines.push_back( SetLineWidth(LINEWIDTH) );

     if ( !COLOR_BY_POINT )
     {    color pointcolor( PR, PG, PB );
          points.push_back( SetColor(pointcolor) );    }

     map<double, double> antisegments;
     if (SEGMENT_FILE != "") {
        Ifstream( sfin0, SEGMENT_FILE );
        double last = 0.0;
        while (!sfin0.eof()) {
           double x1, x2; 
           sfin0 >> x1 >> x2;
           antisegments[last] = x1;
           last = x2;
        }
        sfin0.close();
     }

     Ifstream( in, IN );
     Bool first = True;
     double xm = 0.0, ym = 0.0;
     float rm = -1, gm = -1, bm = -1;
     while(1)
     {    double x, y;
          in >> x >> y;
          float r, g, b;
          if (COLOR_BY_POINT) in >> r >> g >> b;
          if ( !in ) break;
          if ( MAX_Y != undefined && y > MAX_Y ) y = MAX_Y;
          if (COLOR_BY_POINT) points.push_back( SetColor( color( r, g, b ) ) );
          points.push_back( Point( x, y, POINTSIZE ) );
          if ( CONNECT && !first ) 
          {    if ( !COLOR_BY_POINT || ( r == rm && g == gm && b == bm ) )
               {    if (POINTS_LINES_TOGETHER) {
                       if (SEGMENT_FILE == "" || antisegments.size() == 2)
                          points.push_back( Segment( xm, ym, x, y ) );
                       else {
                          // Here we make sure no lines are drawn across segments
                          Bool found = false;
                          for (map<double,double>::iterator miter = antisegments.begin(); miter != antisegments.end(); miter++) {
                             if (xm <= miter->first && x >= miter->second)
                                found = true;
                           }
                           if (!found)
                              points.push_back( Segment( xm, ym, x, y ) ); 
                       }
                    } else {
                            if (SEGMENT_FILE == "" || antisegments.size() == 2) 
                                lines.push_back( Segment( xm, ym, x, y ) );    
                            else {
                                // Here we make sure no lines are drawn across segments
                                Bool found = false;
                                for (map<double,double>::iterator miter = antisegments.begin(); miter != antisegments.end(); miter++) {
                                   if (xm <= miter->first && x >= miter->second)
                                      found = true;
                                }
                                if (!found)
                                   lines.push_back( Segment( xm, ym, x, y ) ); 
                            }
                    }    
               }
          }
          xm = x;
          ym = y;
          rm = r;
          gm = g;
          bm = b;
          first = False; 
     }

     points.push_back( SetColor(black) );
     points.append(lines);

     double x = MinX(points), X = MaxX(points);
     if ( XMIN != undefined ) x = Min( x, XMIN );
     if ( XMAX != undefined ) X = Max( X, XMAX );
     double y = MinY(points), Y = MaxY(points);
     if ( MIN_Y != undefined ) y = MIN_Y;
     if ( MAX_Y != undefined ) Y = MAX_Y;


     vec<graphics_primitive> vert_lines;
     if (VERT_FILE != "") {
        vert_lines.push_back( SetLineWidth(LINEWIDTH) );
        Ifstream( vin, VERT_FILE );
        while (!vin.eof()) {
           double vx;
           vin >> vx;
           vert_lines.push_back( Segment(vx, y, vx, Y, blue) );
        }
        vin.close();
        points.append(vert_lines);
     }

     // I know, the same file is being opened twice. 
     // Artifact of code development.  Fix some day...
     vec<graphics_primitive> seg_lines;
     if (SEGMENT_FILE != "") {
        seg_lines.push_back( SetLineWidth(LINEWIDTH) );
        Ifstream( sfin, SEGMENT_FILE );
        while (!sfin.eof()) {
           double vx1, vx2;
           sfin >> vx1 >> vx2;
           seg_lines.push_back( Segment(vx1, y, vx1, Y, green) );
           seg_lines.push_back( Segment(vx2, y, vx2, Y, magenta) );
        }
        sfin.close();
        points.append(seg_lines);
     }


     vec<graphics_primitive> x_axis 
          = AxisX( x, X, 1.0, True, "", X_AXIS_EXTEND, X_AXIS_NAKED );
     x_axis.push_front( SetLineWidth(LINEWIDTH) );
     vec<graphics_primitive> y_axis 
          = AxisY( x, X, y, Y, Y_MINOR_TICS, "", Y_AXIS_EXTEND, Y_AXIS_NAKED );
     vec<graphics_primitive> points2(y_axis);
     points2.append(points);
     points = points2;

     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(x_axis),    heights.push_back(X_AXIS_OFFSET);
     stack.push_back(points),    heights.push_back(YHEIGHT);
     if ( TITLE != "" )
     {    vec<graphics_primitive> title;
          title.push_back( TextCenter( TITLE, (X+x)/2.0, 0, 0, 
               ArialBold(TITLE_FONTSIZE) ) );
          stack.push_back(title);
          heights.push_back(20);    }
     RenderGraphics( OUT, stack, heights, SCALE, XTRANS, YTRANS, SHOW,
          POSTSCALE, False, PIPEFAIL );    }
