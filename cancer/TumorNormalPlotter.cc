/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// TumorNormalPlotter.  Plot data from tumor and matched normal shotgun data.
// Under development.
//
// IN1, IN2 = outputs of CancerCopy with different arguments for WINDOW
//            (filename suffixes must be .Ni, where Ni = window size, i = 1, 2)
//
// CHR.START-STOP = interval on genome to be plotted

#include "graphics/BasicGraphics.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "math/Functions.h"

const double undefined = 1000000000;

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN1);
     CommandArgument_String(IN2);
     CommandArgument_String(OUT);
     CommandArgument_Double_OrDefault(POINTSIZE, 1.0);
     CommandArgument_UnsignedInt_OrDefault(TITLE_FONTSIZE, 15);
     CommandArgument_String(CHR);
     CommandArgument_String(START);
     CommandArgument_String(STOP);
     CommandArgument_Bool_OrDefault(SHOW_SOURCE, True);
     CommandArgument_String_OrDefault(REARRANGEMENTS, "");
     CommandArgument_Double_OrDefault(POSTSCALE, 1.0);
     EndCommandArguments;

     CheckValidGraphicsSuffix(OUT);

     vec<double> rearrangements;
     if ( REARRANGEMENTS != "" ) ParseDoubleSet( REARRANGEMENTS, rearrangements );

     int start = int( floor( START.Double( ) ) ); 
     int stop = int ( floor( STOP.Double( ) ) );
     int WIDTH1 = IN1.substr( IN1.PosRev( "." ) + 1 ).Int( );
     int WIDTH2 = IN2.substr( IN2.PosRev( "." ) + 1 ).Int( );
     ForceAssert( (start % WIDTH2) == 0 );
     ForceAssert( (stop % WIDTH2) == 0 );
     ForceAssertLt( start, stop );
     ForceAssert( WIDTH2 % WIDTH1 == 0 );
     double width1 = double(WIDTH1)/1000000.0; 
     double width2 = double(WIDTH2)/1000000.0; 

     String START0 = START;
     if ( START0.Contains( "M", -1 ) ) START0 = START0.Before( "M" );
     String TITLE = START0 + " - " + STOP + " on chromosome " + CHR;

     vec<graphics_primitive> points, points2;
     points.push_back( SetLineWidth(1.0) );

     color pointcolor( 0, 0, 0 );
     points.push_back( SetColor(pointcolor) );

     String line;
     static vec<String> tokens;
     for ( int pass = 1; pass <= 2; pass++ )
     {    Ifstream( in, ( pass == 1 ? IN1 : IN2 ) );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "chr", 0 ) ) continue;
               if ( !line.Contains( "chr" + CHR + ":", 0 ) ) continue;
               line.GlobalReplaceBy( "-", "" );
               line.ReplaceBy( ":", ": " );
               Tokenize( line, tokens );
               double x = tokens[1].Before( "," ).Double( );
               if ( x < start || x >= stop ) continue;
               if (tokens[9] == "nan" || tokens[10] == "nan" || tokens[11] == "nan,")
                    continue;
               if (tokens[9] == "inf" || tokens[10] == "inf" || tokens[11] == "inf,")
                    continue;
               x /= 1000000.0;
               double y1 = tokens[9].Double( );
               double y2 = tokens[10].Double( );
               double y3 = tokens[11].Before( "," ).Double( );
               if ( pass == 1 )
               {    points.push_back( RectangleBaseCentered( 
                         x + width1/2.0, y1, y3-y1, width1, lightergray, False ) );
                    points2.push_back( 
                         Segment( x, y2, x + width1, y2, darkergray ) );    }
               else
               {    points.push_back( RectangleBaseCentered( 
                         x + width2/2.0, y1, y3-y1, width2, lighterpink, False ) );
                    points.push_back( 
                         Segment( x, y2, x + width2, y2, red ) );    }    }    }

     points.push_back( SetColor(black) );
     points.append(points2);

     double x = MinX(points), X = MaxX(points);
     double y = MinY(points), Y = MaxY(points);

     vec<graphics_primitive> vlines;
     for ( int i = 0; i < rearrangements.isize( ); i++ )
          vlines.push_back( Segment( rearrangements[i]/1000000.0, 0,
               rearrangements[i]/1000000.0, Y, blue ) );
     vlines.append(points);
     points = vlines;

     vec<graphics_primitive> x_axis0 = AxisX( x, X, 1.0, True, "", 0.0 );
     vec<graphics_primitive> x_axis;
     x_axis.push_back( SetLineWidth(0.25) );
     x_axis.append(x_axis0);
     vec<graphics_primitive> y_axis = AxisY( x, X, 0, Y, True, "", 0.0 );
     points.push_back( SetLineWidth(0.25) );
     points.append(y_axis);

     vec< vec<graphics_primitive> > stack;
     vec<double> heights;

     if (SHOW_SOURCE)
     {    vec<graphics_primitive> source;
          source.push_back( TextCenter( "generated by " + command.TheCommand( ), 
               (X+x)/2.0, 0, 0, TimesRoman(5) ) );
          stack.push_back(source),    heights.push_back(30);    }

     stack.push_back(x_axis),    heights.push_back(0);
     stack.push_back(points),    heights.push_back(200);
     if ( TITLE != "" )
     {    vec<graphics_primitive> title;
          title.push_back( TextCenter( TITLE, (X+x)/2.0, 0, 0, 
               TimesBold(TITLE_FONTSIZE) ) );
          stack.push_back(title);
          heights.push_back(20);    }
     RenderGraphics( OUT, stack, heights, 1.0, 200, 50, True, POSTSCALE );    }
