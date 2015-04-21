/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// PlotGC.  Plot GC content of a genomic region.
//
// GENOME = fastb file for genome
// TIG = index of genome contig in fastb file
// START = start on contig
// STOP = stop on contig
// WINDOW = window size to average over
// OUT = graphics output file name

#include "Basevector.h"
#include "graphics/BasicGraphics.h"
#include "MainTools.h"
#include "math/Functions.h"

const double undefined = 1000000000;

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(GENOME);
     CommandArgument_UnsignedInt(TIG);
     CommandArgument_UnsignedInt(START);
     CommandArgument_UnsignedInt(STOP);
     CommandArgument_UnsignedInt(WINDOW);
     CommandArgument_Int_OrDefault(GRANULARITY, 1);
     CommandArgument_String(OUT);
     CommandArgument_Double_OrDefault(SCALE, 1.0);
     CommandArgument_Double_OrDefault(PR, 0);
     CommandArgument_Double_OrDefault(PG, 0);
     CommandArgument_Double_OrDefault(PB, 0);
     CommandArgument_String_OrDefault(TITLE, "");
     CommandArgument_Double_OrDefault(XTRANS, 200);
     CommandArgument_Double_OrDefault(YTRANS, 50);
     CommandArgument_Bool_OrDefault(SHOW, True);
     CommandArgument_Double_OrDefault(MAX_Y, undefined);
     CommandArgument_UnsignedInt_OrDefault(TITLE_FONTSIZE, 15);
     CommandArgument_Double_OrDefault(LINEWIDTH, 1.0);
     CommandArgument_String_OrDefault(TRANSPARENT, "");
     EndCommandArguments;

     Bool to_png = OUT.Contains( ".png", -1 );
     Bool to_ps = OUT.Contains( ".ps", -1 );
     Bool to_pdf = OUT.Contains( ".pdf", -1 );
     ForceAssert( to_png || to_ps || to_pdf );
     if ( to_png && !SHOW )
     {    cout << "For a png file, you must specify SHOW=True.\n";
          exit(-1);    }

     vec<graphics_primitive> points, lines;
     lines.push_back( SetLineWidth(LINEWIDTH) );

     color pointcolor( PR, PG, PB );
     points.push_back( SetColor(pointcolor) );

     vecbasevector genome(GENOME);
     ForceAssertLt( TIG, genome.size( ) );
     ForceAssert( START >= WINDOW/2 );
     const basevector& g = genome[TIG];
     double xm = 0.0, ym = 0.0;
     ForceAssert( STOP <= g.size( ) - WINDOW/2 );
     for ( unsigned int p = START; p < STOP; p += WINDOW/GRANULARITY )
     {    int gc = 0;
          for ( unsigned int j = p - WINDOW/2; j < p + WINDOW/2; j++ )
               if ( g[j] == 1 || g[j] == 2 ) ++gc;
          double x = p;
          double y = 100.0 * double(gc) / double(WINDOW);
          if ( MAX_Y != undefined && y > MAX_Y ) y = MAX_Y;
          if ( p > START ) lines.push_back( Segment( xm, ym, x, y ) );
          xm = x;
          ym = y;    }

     points.append(lines);
     points.push_back( SetColor(black) );

     double x = MinX(points), X = MaxX(points);
     double y = MinY(points), Y = MaxY(points);
     if ( MAX_Y != undefined ) Y = MAX_Y;
     vec<graphics_primitive> x_axis = AxisX( x, X );
     vec<graphics_primitive> y_axis = AxisY( x, X, y, Y );
     points.append(y_axis);

     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(x_axis);    
     heights.push_back(30);
     stack.push_back(points),    heights.push_back(200);
     if ( TITLE != "" )
     {    vec<graphics_primitive> title;
          title.push_back( TextCenter( TITLE, (X+x)/2.0, 0, 0, 
               TimesBold(TITLE_FONTSIZE) ) );
          stack.push_back(title);
          heights.push_back(20);    }
     String outhead;
     if (to_ps) outhead = OUT.Before( ".ps" );
     else if (to_pdf) outhead = OUT.Before( ".pdf" );
     else outhead = OUT.Before( ".png" );
     double width = 400;
     {    Ofstream( out, outhead + ".ps" );
          VerticalDisplay( stack, heights, width, out, SCALE, XTRANS, 
               YTRANS, SHOW );    }
     if (to_pdf)
     {    System( "ps2pdf " + outhead + ".ps" );
          Remove( outhead + ".ps" );    }
     if (to_png)
     {    String command                            
               = "pstopnm -portrait -xmax 8000 -ymax 8000 -stdout "
               + outhead + ".ps | pnmcrop | "
               + "pnmpad -white -left=25 -right=25 -top=25 -bottom=25 | "
               + "pnmtopng";
          if ( TRANSPARENT == "WHITE" ) command += " -transparent =white";
          command += " > " + outhead + ".png";
          int status = System(command);
          if ( status == 0 ) Remove( outhead + ".ps" );
          else cout << "failed to run:\n" << command << "\n";     }    }
