// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// This code generates a figure from insert information for a library.  It is
// being used to generate a specific figure for the finished human genome paper.
// It also functions as a test code for BasicGraphics (which is why it is in this 
// directory).

// More specifically, this program takes as input a file of ordered pairs, each
// pair representing the starting and stopping positions of an insert on a 
// chromosome.  It generates a figure from this.  (More explanation to be added 
// later.)

// This file is still under construction.

/*
InsertPlotter PAIRS=graphics/data_BW3.1 SITE=46969320 CHR=10 \
     TITLE="BW3.1, 5.8 kb deletion, common polymorphism"
InsertPlotter PAIRS=graphics/data_BW17 CHR=11 TITLE="BW17, no evidence of deletion"
*/

#include "graphics/BasicGraphics.h"
#include "MainTools.h"
#include "math/Functions.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PAIRS);
     CommandArgument_String(TITLE);
     CommandArgument_String(CHR);
     CommandArgument_UnsignedInt_OrDefault(SITE, 0);
     EndCommandArguments;

     vec<int> start, stop;
     Ifstream( in, PAIRS );
     while(1)
     {    int a, b;
          in >> a >> b;
          if ( !in ) break;
          start.push_back(a), stop.push_back(b);    }
     int n = start.size( ), min_start = Min(start);
     for ( int i = 0; i < n; i++ )
     {    start[i] -= min_start, stop[i] -= min_start;    }

     const int insert_mean_standard = 39924;
     const int insert_dev_standard = 2735;

     float dl = SITE - min_start;
          
     vec<graphics_primitive> inserts;
     for ( int i = 0; i < n; i++ )
          inserts.push_back( Segment( start[i], 100 * i, stop[i], 100 * i ) );
     if ( SITE != 0 )
          inserts.push_back( Segment( dl, -100, dl, 100 * (n+2), cyan ) );
     inserts.push_back( SetTimesRoman(10) );
     for ( int i = 0; i < n; i++ )
     {    int len = stop[i] - start[i];
          String delta_s = ToString( len - insert_mean_standard );
          if ( isdigit( delta_s[0] ) ) delta_s = "+" + delta_s;
          inserts.push_back( TextToRight( delta_s, stop[i], 100 * i, 5 ) );    }
     if ( SITE != 0 )
     {    inserts.push_back( 
               TextToRight( "deletion location", dl, 100 * (n+2), 5, 
                    TimesRoman(12) ) );    }

     double x = MinX(inserts), X = MaxX(inserts);
     vec<graphics_primitive> x_axis = AxisX( x, X, 1000 );
     x_axis.push_back( TextToRight( "kb", MaxX(x_axis), 0, 5, TimesRoman(12) ) );

     vec<graphics_primitive> devs;
     double last_p = -1.0, last_delta = -1.0;
     // for ( double p = x; p <= X; p += 500 )
     for ( double p = stop.front( ); p <= start.back( ); p += 500 )
     {    vec<int> vc;
          for ( int i = 0; i < n; i++ )
          {    if ( start[i] > p || stop[i] < p + 500 ) continue;
               vc.push_back( stop[i] - start[i] );    }
          if ( vc.size( ) >= 2 )
          {    sort( vc.rbegin( ), vc.rend( ) );
               vc.resize( vc.size( ) - 1 );
               double delta = -(Mean(vc) - insert_mean_standard)
                    / ( insert_dev_standard / sqrt( double(vc.size( ) ) ) );
               if ( last_p >= 0.0 )
                    devs.push_back( Segment( last_p, last_delta, p, delta ) );
               last_p = p, last_delta = delta;    }    }
     devs.push_back( DottedSegment( x, 3.5, X, 3.5, red ) );
     devs.push_back( TextToRight( "3.5 dev", X, 3.5, 5, TimesRoman(12) ) );

     double x1 = Min( MinX(inserts), MinX(devs) ); 
     double X2 = Max( MaxX(inserts), MaxX(devs) );
     double y1 = MinY(devs), Y2 = MaxY(devs);
     devs.append( AxisY( x1, X2, y1, Y2, False ) );

     vec<graphics_primitive> title, title2, title3;
     title.push_back( TextToRight( TITLE, MinX(inserts), 0, 0, TimesRoman(15) ) );
     title2.push_back( TextToRight( 
          "chromosome " + CHR + ", window starts at " + ToString(min_start),
          MinX(inserts), 0, 0, TimesRoman(15) ) );
     title3.push_back( TextToRight( " ", MinX(inserts), 0, 0, TimesRoman(15) ) );

     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(devs),    heights.push_back(100);
     stack.push_back(x_axis),  heights.push_back(50);
     stack.push_back(inserts), heights.push_back( 11 * n );
     stack.push_back(title3),  heights.push_back(30);
     stack.push_back(title2),  heights.push_back(20);
     stack.push_back(title),   heights.push_back(20);
     double w = 300;
     {    Ofstream( out, "InsertPlotter.ps" );
          VerticalDisplay( stack, heights, w, out );    }
     // System( "ghostview InsertPlotter.ps" );    
     System( "ps2pdf InsertPlotter.ps" );
          }
