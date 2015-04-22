// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//

#include "graphics/BarGraph.h"
#include "math/Functions.h"
#include "graphics/GraphPlotsForDummies.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"

using namespace ns_psplot;



/*
 * GenerateGraphPlot
 */
int GenerateGraphPlot( const vec<freetext> &labels,
		       const vec<float> &x_points,
		       const vec<float> &y_points,
		       const String &file_name )
{
  // Default setting.
  //  min_points: min number of points in the vectors x and y;
  //  max_points: max number of points in x and y.  
  //  graph_width: graph width;
  //  graph_height: height;
  //  col: color for the bars (see BarGraph.h;)
  //  font_size: font size for the tags (tic labels.)
  int min_points = 5;
  // int max_points = 1000000;
  int graph_width = 400;
  int graph_height = 250;
  color col = black;
  short font_size = 7;

  // Check vector sizes.
  if ( x_points.size() != y_points.size() )
    return 0;

  if ( (int)x_points.size() < min_points )
       // || (int)x_points.size() > max_points )
    return 0;
  
  // x must be increasing, and x[0]>=0.
  if ( x_points[0] < 0 )
    return 0;
  for (int ii=0; ii<(int)x_points.size()-1; ++ii)
    if ( x_points[ii] > x_points[ii+1])
      return 0;

  // y must be >= 0, and y cannot be a constant vector.
  float max_y = y_points[0];
  float min_y = y_points[0];
  for (int ii=0; ii<(int)y_points.size(); ++ii) {
    if ( y_points[ii] < 0 )
      return 0;

    if ( y_points[ii] < max_y )
      min_y = y_points[ii];
    
    if ( y_points[ii] > max_y )
      max_y = y_points[ii];
  }
  if ( max_y <= min_y )
    return 0;
  
  // Create the file.
  using namespace ns_psplot;

  String save_file = file_name + ".eps";
  ofstream f_fig( save_file.c_str() );  

  // Local copy of (const) x and y points.
  vec<float> x_pts = x_points;
  vec<float> y_pts = y_points;

  plotgraph graph = GraphPlot( x_pts, y_pts, col, font_size, labels,
			       graph_width, graph_height );
			       
  
  f_fig << graph << "\n";
  
  //  All right, can return.
  return 1;
}



