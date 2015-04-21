/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <numeric>
#include <strstream>

#include "graphics/BarGraph.h"
#include "math/Functions.h"
#include "Vec.h"
#include "system/Assert.h"

using namespace ns_psplot;

ostream& ns_psplot::operator<<( ostream& o, const tic& t )
{    o << translate(0, t.T( ));

     // Spacing.  I don't understand it.  Hardcoded to work for angle 0 and 90.

     float spacex, spacey;
     float F = t.Label( ).FontSize( );
     if ( t.Angle( ) == 0 )
     {    spacex = 0;
          spacey = F/3;    }
     else 
     {    spacex = F/3;
          spacey = F/2;    }

     // make shorter length for tics with no label printed:
     int tLength = ( t.Label().Str().size()==0 ) ? TicLength/2 : TicLength;
     if ( t.dir_ == left )
          o << startpath( -tLength, 0 ) << rlineto( tLength, 0 ) << " stroke"
               << translate( -tLength - TicLabelSep, 0 ) 
               << " " << t.Angle( ) << " rotate "
               << translate( -spacex, -spacey )
               << text( t.Label( ), right )
               << translate( spacex, spacey )
               << " " << -t.Angle( ) << " rotate "
               << translate( tLength + TicLabelSep, 0 ) 
                    << "\n";
     else if ( t.dir_ == right )
          o << startpath(0, 0) << rlineto( tLength, 0 ) << " stroke"
               << translate( tLength + TicLabelSep, 0 )
               << " " << t.Angle( ) << " rotate "
               << translate( -spacex, -spacey )
               << text( t.Label( ), left )
               << translate( spacex, spacey )
               << " " << -t.Angle( ) << " rotate "
               << translate( -tLength - TicLabelSep, 0 ) << "\n";
     o << translate(0, -t.T( ));
     return o;    }

ostream& ns_psplot::operator<<( ostream& o, const bar& b )
{
  if ( 0 == b.y_ )
    return o;
  float w = b.width_ - EdgeFudge, y = b.y_ - EdgeFudge;
  return o << b.c_
	   << startpath( b.x_ - w/2, 0 )
	   << rlineto(w, 0)
	   << rlineto(0, y)
	   << rlineto(-w, 0)
	   << rlineto(0, -y)
	   << " fill";
}


ostream& ns_psplot::operator<<( ostream& o, const bargraph& g)
{
  //  extra_space: a few pixels extra space below the x axis (left of y axis.)
  int extra_space = 10;

  // Sides of the eps BoundingBox. Remark: this is only an approximated
  // value, since I could not determine exactly the size of a given printed
  // string. I use the following heuristics: the width of the string '00'
  // (two zeroes,) is approximately equal to FontSize.
  int n_char_left = g.left_tics_.back( ).Label( ).Str( ).size( );
  int left_factor = ( n_char_left + 2 ) / 2;
  int n_char_bottom = g.bottom_tics_.back( ).Label( ).Str( ).size( );
  int bottom_factor = ( n_char_bottom + 2 ) / 2;
  int tic_label_length_bottom
    = bottom_factor * g.bottom_tics_[0].Label( ).FontSize( );
  int tic_label_length_left
    = left_factor * g.left_tics_[0].Label( ).FontSize( );
  int left_boundbox = -extra_space -TicLength -TicLabelSep -tic_label_length_left;
  
  int right_boundbox = int(g.Width()) + (n_char_left/2);
  
  int bottom_boundbox = -extra_space -TicLength -TicLabelSep -tic_label_length_bottom;
  
  int top_boundbox = int(g.Height()) + 1;
  
  // Sante --- Fri Jun 22 15:29:47 EDT 2001
  // If g.use_eps_ is true, then output is an eps file.
  if ( g.use_eps_ ) {
    
    o << "%!PS-Adobe-3.0 EPSF-3.0\n"
      << "%%BoundingBox: "
      << "0 0 "
      << right_boundbox - left_boundbox << " " << top_boundbox - bottom_boundbox << "\n";
    
    // Set line width.
    o << "0.2 setlinewidth\n";
    
    // Translate boundbox to origin
    o << translate( -left_boundbox, -bottom_boundbox ) << "\n";

  } // Sante --- End "if ( g.use_eps_ )".x
  
  o << translate( -g.LeftBoundary( ), 0 ) << " -90 rotate";
  
  for ( unsigned int i = 0; i < g.bottom_tics_.size( ); i++ ) {
    o << black << tic( g.bottom_tics_[i], right );
  }

  o << " 90 rotate";
  o << translate( g.LeftBoundary( ), 0 );
  
  // Sante --- Mon Jun 25 10:35:44 EDT 2001
  // Plot tics and horizontal "rules."
  for ( unsigned int i = 0; i < g.left_tics_.size( ); i++ )
  {
    o << black << tic( g.left_tics_[i], left );
    
    o << lightgray
      << startpath(0, g.left_tics_[i].T())
      << rlineto(g.Width( ), 0)
      << " stroke\n";
  }
  // End_Sante
  
  o << translate( -g.LeftBoundary( ), 0 );
  for ( unsigned int i = 0; i < g.bars_.size( ); i++ )
    o << g.bars_[i];
  
  o << translate( g.LeftBoundary( ) + g.Width( )/2, g.BarTopBoundary( ) );
  o << translate( 0, LabelSep );     
  for ( int i = g.labels_.size( ) - 1; i >= 0; i-- )
    o << translate( 0, LabelSep ) << text( g.labels_[i], center )
      << translate( 0, g.labels_[i].FontSize( ) );
  
  for ( int i = g.labels_.size( ) - 1; i >= 0; i-- )
    o << translate( 0, -LabelSep - g.labels_[i].FontSize( ) );
  
  o << translate( - g.Width( )/2, -LabelSep - g.BarTopBoundary( ) );

  o << black << startpath(0, 0) << rlineto( g.Width( ), 0 )
    << rlineto( 0, g.Height( ) ) << rlineto( -g.Width( ), 0 )
    << rlineto( 0, -g.Height( ) ) << " stroke\n";
  
  return o;    }


// Sante --- Tue Jun 12 13:48:08 EDT 2001
ostream& ns_psplot::operator<<( ostream& o, const plotgraph& plotg )
{
  // Legenda:
  //  extra_space: a few pixels extra space below the x axis (left of y axis.)
  int extra_space = 0;

  // Sides of the eps BoundingBox. Remark: this is only an approximated
  // value, since I could not determine exactly the size of a given printed
  // string. I use the following heuristics: the width of the string '00'
  // (two zeroes,) is approximately equal to FontSize.
  int n_char_left = plotg.left_tics_[plotg.left_tics_.size()-1].Label( ).Str().size();
  int left_factor = (n_char_left+2)/2;
  
  int n_char_bottom
    = plotg.bottom_tics_[plotg.bottom_tics_.size()-1].Label( ).Str().size();
  int bottom_factor = (n_char_bottom+2)/2;

  int tic_label_length_bottom
    = bottom_factor * plotg.bottom_tics_[0].Label( ).FontSize( );
  int tic_label_length_left
    = left_factor * plotg.left_tics_[0].Label( ).FontSize( );

  int left_boundbox = -extra_space -TicLength -TicLabelSep -tic_label_length_left;

  int right_boundbox = int(plotg.Width()) + (n_char_left/2);

  int bottom_boundbox = -extra_space -TicLength -TicLabelSep -tic_label_length_bottom;

  int top_boundbox = int(plotg.Height()) + 1;

  o << "%!PS-Adobe-3.0 EPSF-3.0\n"
    << "%%BoundingBox: "
    << "0 0 "
    << right_boundbox - left_boundbox << " " << top_boundbox - bottom_boundbox << "\n";

  // Set line width.
  o << "0.2 setlinewidth\n";

  // Frame.
  o << black 
    << startpath( 0, 0 ) 
    << rlineto( right_boundbox-left_boundbox, 0 )
    << rlineto( 0, top_boundbox-bottom_boundbox ) 
    << rlineto( -right_boundbox+left_boundbox, 0 )
    << rlineto( 0, -top_boundbox+bottom_boundbox ) << " stroke\n";

  // Translate boundbox to origin
  o << translate( -left_boundbox, -bottom_boundbox );

  // Plot tics and rules.
  o << translate( -plotg.LeftBoundary( ), 0 ) << " -90 rotate";
  
  o << translate (+extra_space, plotg.pixels_per_unit_x_);
  for ( unsigned int i = 0; i < plotg.bottom_tics_.size( ); i++ )
    {
      o << black << tic(plotg.bottom_tics_[i], right);
      if ( i == 0 )
	o << black;
      else
	o << lightgray;
      o << startpath(-extra_space -plotg.PlotTopBoundary() , plotg.bottom_tics_[i].T())
	<< rlineto(+extra_space +plotg.PlotTopBoundary(), 0)
	<< " stroke\n";
    }
  o << translate (-extra_space, -plotg.pixels_per_unit_x_);

  o << " 90 rotate";
  o << translate( plotg.LeftBoundary( ), 0 );

  o << translate (-extra_space, 0);  
  for ( unsigned int i = 0; i < plotg.left_tics_.size( ); i++ )
    {
      o << black << tic( plotg.left_tics_[i], left ); 
      if ( i == 0 )
	o << black;
      else
	o << lightgray;
      o << startpath(0, plotg.left_tics_[i].T())
	<< rlineto(extra_space +plotg.Width(), 0)
	<< " stroke\n";
    }
  o << translate (+extra_space, 0);
  
  // Plot segments.
  o << translate( -plotg.LeftBoundary( ) +plotg.pixels_per_unit_x_, 0 );

  o << plotg.segments_[0].Color() << "\n";
  o << startpath( plotg.segments_[0].X1(), plotg.segments_[0].Y1() ) << "\n";

  for ( unsigned int i = 0; i < plotg.segments_.size(); i++ )
    {
      o << lineto( plotg.segments_[i].X2(), plotg.segments_[i].Y2() ) << "\n";
    }
  o << "stroke\n";

  o << translate( -plotg.pixels_per_unit_x_, 0 );

  // Labels.
  o << translate( plotg.LeftBoundary( ) + plotg.Width( )/2, plotg.PlotTopBoundary( ) );
  o << translate( 0, LabelSep );     
  for ( int i = plotg.labels_.size( ) - 1; i >= 0; i-- )
    o << translate( 0, LabelSep ) << text( plotg.labels_[i], center )
      << translate( 0, plotg.labels_[i].FontSize( ) );
  
  for ( int i = plotg.labels_.size( ) - 1; i >= 0; i-- )
    o << translate( 0, -LabelSep - plotg.labels_[i].FontSize( ) );
  
  // Return.
  return o;
}



// Sante --- Thu Jun 14 13:38:36 EDT 2001
//  zero_tic: set it to true if you want to print also 0 among the tics.
vec<freetic> 
MakeTics( float min_value, float max_value, float max_bar_height, color c, 
	  short font_size, bool zero_tic = false, const int digits_to_right = 4 )
{
  const float default_delta = 0.1;  // used only if graph is empty.

  float delta = min_value < max_value ? max_value - min_value : default_delta;
  float top = pow( 10, floor(log10(delta)) );
  float atop = top;
  int first_digit = int(round(delta/top));
  
  if ( first_digit == 1 )
    {
      first_digit = 10;
      atop /= 10;
    }
  
  int tics = first_digit;
  while( tics * atop <= delta )
    ++tics;
  
  if ( !zero_tic )
    --tics;
  
  vec<freetic> answer( tics );

  if ( !zero_tic )
    for ( int i = 1; i <= tics; i++ )
      {
	float d_tic = ( min_value + i * atop ) / delta * max_bar_height;
	String sz_tic = ToString( min_value + i*atop, digits_to_right );
	answer[i-1] = freetic( d_tic, sz_tic, c, font_size );      
      }
  else
    {
      float d_tic = min_value;
      String sz_tic = ToString( min_value, 4 );
      answer[0] = freetic( d_tic, sz_tic, c, font_size );      
      
      for ( int i = 1; i < tics; i++ )
	{
	  float d_tic = ( min_value + i * atop ) / delta * max_bar_height;
	  String sz_tic = ToString( min_value + i*atop, 4);
	  answer[i] = freetic( d_tic, sz_tic, c, font_size );      
	}
    }
  
  return answer;
}

// Sante --- Fri Jun 22 14:53:27 EDT 2001
//  use_eps: if true creates eps output (default is false.)
bargraph ns_psplot::BarPlot(vec<float> y, float low, float high, 
			    float bar_width, float bar_sep, color c,
			    short font_size, float max_bar_height,
			    vec<freetext> labels, vec<float> f_bottom_tics,
			    bool use_eps, int digits_to_right, bool bottom_tics_horizontal )
{    vec<bar> bars;
     Assert( y.size( ) > 0 );
     for ( unsigned int i = 0; i < y.size( ); i++ )
     {    Assert( y[i] >= 0 );
          bars.push_back( bar( i * (bar_width + bar_sep), 
               max_bar_height * y[i]/Max(y), bar_width, c ) );    }
     vec<freetic> bottom_tics( y.size( ) ), left_tics;

     // Plot user-given bottom_tics, if given.
     if ( f_bottom_tics.size() == y.size() ) {
       for (int ii=0; ii<(int)y.size(); ++ii) {
	  string x = ToString( f_bottom_tics[ii], digits_to_right );	  
	  bottom_tics[ii] = freetic( ii * (bar_width + bar_sep),
				     x, black, font_size );
       }
     } else {
       float skip = 0;
       for ( unsigned int i = 0; i < y.size( ); i++ ) {
	   string x = ToString( low + (high - low) * float(i) 
				/ float(y.size( ) - 1) , digits_to_right );
	   skip -= (bar_width+bar_sep);
	   if ( skip > 0 ) {
	     bottom_tics[i] = freetic( i * (bar_width + bar_sep), "", 
				       black, font_size );
	   } else {
	     if ( bottom_tics_horizontal ) {
	       bottom_tics[i] = freetic( i * (bar_width + bar_sep), x, 
					 black, font_size,90 );
	       skip = font_size*x.size()+2;
	     } else {
	       bottom_tics[i] = freetic( i * (bar_width + bar_sep), x, 
					 black, font_size );
	       skip = 0;
	     }
	   }
       }
     } 
     left_tics = MakeTics( 0, Max(y), max_bar_height, black, font_size,
          false, digits_to_right  );
     return bargraph( bars, labels, bottom_tics, left_tics, use_eps );    }


/*
 * Sante --- Tue Jun 12 11:40:12 EDT 2001
 *
 * Plots a function graph y=f(x), for two given vectors x[ii], and y[ii],
 * i.e. it plots segments joining (x[i], y[i]) with (x[i+1], y[i+1]).
 * Drawback: it can only plot graphs in the first quadrant!
 *
 * Legenda: 
 *  x: vector of x points;
 *  y: vector of y points;
 *  col: color;
 *  font_size: tic labels font size;
 *  labels: graph labels;
 *  graph_width: width of the graph;
 *  graph_height: hieght of the graph;
 *  plot_dots: if true, plots a bullet on each point (default = false.
 *             Not implemented yet.)
 */
plotgraph ns_psplot::GraphPlot(vec<float>& x, vec<float>&y, color col,
			       short font_size, vec<freetext> labels,
			       int graph_width, int graph_height,
			       bool plot_dots)
{
  // Vector sizes.
  Assert( y.size( ) == x.size() );
  int n_points = (int)x.size();
  
  // Sante --- Wed Mar 13 17:42:55 EST 2002
  //  I change strictly increasing in just increasing.
  // x must be a strictly increasing sequence, and x[0]>=0.
  Assert ( x[0] >= 0 );  
  for (int ii=0; ii<n_points-1; ++ii)
    Assert ( x[ii] <= x[ii+1]);  

  float max_x = x[n_points-1];
  float width_multiplier = x[n_points-1] - x[0];
  float pixels_per_unit_x = ( 1.0 / width_multiplier ) * float(graph_width);

  // y[ii] must be >= 0, and y cannot be a constant vector.
  float max_y = y[0];
  float min_y = y[0];
  for (int ii=0; ii<n_points; ++ii)
    {
      Assert ( y[ii] >= 0 );

      if ( y[ii] < max_y )
	min_y = y[ii];
      
      if ( y[ii] > max_y )
	max_y = y[ii];
    }

  float y_order = pow(10, floor(log10(max_y)));
  if ( max_y / y_order < 1.1 && y_order-1 > 0 )
    y_order--;
  max_y = floor( ( max_y / y_order ) + 1 ) * y_order;

  Assert ( max_y > min_y );
  float pixels_per_unit_y = ( 1.0 / max_y ) * float(graph_height);

  // Create segments.
  vec<segment> segments;

  // For large collections of data, we don't want enormous numbers of
  // segments, so we need to reduce the granularity of the data.

  // First, we define the minimum perceptible difference for both axes in the units of that axis.

  // Define the minimum difference as 1/3 of a pixel, i.e. ~1/200 of an inch.
  float min_x_difference = 1.0 / 3.0;
  float min_y_difference = 1.0 / 3.0;

  // Next, we add consecutive segments together until we pass that
  // threshold on one axis or the other.

  // We start at the first data point.
  segment current_segment( pixels_per_unit_x * x[0],
			   pixels_per_unit_y * y[0], 
			   pixels_per_unit_x * x[0], 
			   pixels_per_unit_y * y[0], col );

  for ( int ii=1; ii<n_points; ++ii )
  {
    // We then try to extend the current segment by merging it with the next.
    segment next_segment( pixels_per_unit_x * x[ii-1],
			  pixels_per_unit_y * y[ii-1],
			  pixels_per_unit_x * x[ii],
			  pixels_per_unit_y * y[ii], col );

    // We know we'll be able to merge them, so we don't check the return value.
    current_segment.merge(next_segment);

    // We check if the segment is long enough to use.
    if ( current_segment.deltaX() > min_x_difference ||
	 current_segment.deltaY() > min_y_difference )
    {
      // If it is, we store it, then start a new segment at the end of the last one.
      segments.push_back( current_segment );
      current_segment = segment( pixels_per_unit_x * x[ii],
				 pixels_per_unit_y * y[ii],
				 pixels_per_unit_x * x[ii],
				 pixels_per_unit_y * y[ii], col );
    }
    // If it isn't, we move on to the next segment.
  }

  // Tics.
  vec<freetic> left_tics;
  vec<freetic> bottom_tics;

  left_tics = MakeTics( 0, max_y, graph_height, black, font_size, true );
  bottom_tics = MakeTics( x[0], max_x, graph_width, black, font_size, true );

  // Return.
  return plotgraph(segments, labels, bottom_tics, left_tics, 
		   pixels_per_unit_x, pixels_per_unit_y, max_y, plot_dots);
}



bargraph ns_psplot::BiBarPlot( vec<float> y1, vec<float> y2, float low, 
			       float high, float bar_width, float bar_sep,
			       color c1, color c2, short font_size,
			       float max_bar_height, vec<freetext> labels, bool use_eps)
{    vec<bar> bars;
     Assert( y1.size( ) > 0 && y1.size( ) == y2.size( ) );
     float M = Max( Max(y1), Max(y2) );
     Assert( M > 0 );
     for ( unsigned int i = 0; i < y1.size( ); i++ )
     {    Assert( y1[i] >= 0 && y2[i] >= 0 );
          bars.push_back( bar( i * (2*bar_width + bar_sep) - bar_width/2, 
               max_bar_height * y1[i]/M, bar_width, c1 ) );
          bars.push_back( bar( i * (2*bar_width + bar_sep) + bar_width/2, 
               max_bar_height * y2[i]/M, bar_width, c2 ) );    }
     vec<freetic> bottom_tics( y1.size( ) ), left_tics;
     for ( unsigned int i = 0; i < y1.size( ); i++ )
     {    string x = ToString( low + (high - low) * float(i) / float(y1.size( )-1) , 4);
          bottom_tics[i] = freetic( i * (2*bar_width + bar_sep), x, black, 
               font_size );    }
     left_tics = MakeTics( 0, M, max_bar_height, black, font_size );
     return bargraph( bars, labels, bottom_tics, left_tics, use_eps );    }
  
  

