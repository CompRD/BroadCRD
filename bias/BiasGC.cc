/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "bias/BiasGC.h"



/// Generates GC plot and saves it into a file.
///
/// @param file_name name of the file to save the plot into
/// @param curves vector of precomputed curves to plot (x=GC percentage,
/// y=normalized frequency of reads); each curve's title (if it is set)
/// will be shown in the plot's legend
/// @param x_label_text Text to print for the X axis label (if specified)
/// @param colors if specified, colors[i] color will be used for curves[i]
/// curve; if not specified up to 8 default colors will be used (can not
/// pass more than 8 curves in this case!)
void GenerateGCPlot(String & file_name,
	     const vec<Curve> & curves,
	     const char * x_label_text,
	     const vec<color> * colors,
             const double TITLE_FONTSIZE ) {

    vec<graphics_primitive> points;
    vec<color> cols;
    if ( colors == 0 ) {
        // set 8 default colors
        cols.push_back( red, blue, green );
        cols.push( 1, 0.25, 1 );
        cols.push( 1, 0.85, 0 );
        cols.push( 0, 0.8, 0.8 );
	cols.push( 0.65, 0.65, 0.65 );
	cols.push( 0.7, 0.3, 0.4 );
    } else {
        cols = (*colors); // if colors are provided, copy them
    }
     
    if ( curves.size() > cols.size() ) {
        PRINT2(curves.size(),cols.size());
	FatalErr("ERROR::too many curves, too few colors");
    }

    const double POINTSIZE = 1.0;

    for ( unsigned int i = 0 ; i < curves.size() ; i++ ) {
        points.push_back( SetColor( cols[i] ) );
	for ( unsigned int j = 0; j < curves[i].NPoints(); j++ ) {
	    points.push_back( Point( curves[i].x[j], 
				     curves[i].y[j], 
				     POINTSIZE ) );    
	}
    }
    
    vec<graphics_primitive> title;
    vec<String> heads0;
    for ( unsigned int i = 0; i < curves.size( ); i++ ) {
        heads0.push_back( curves[i].title + " " );
    }
    cols.resize( curves.size( ) );
    title.push_back( RainbowTextCenter( heads0, cols, 100.0/2.0, 0, 0,
          TimesBold(TITLE_FONTSIZE) ) );

     // Generate graph.

     double LINEWIDTH = 0.5;
     double POSTSCALE = 0.2;
     points.push_back( SetColor(black), SetLineWidth(LINEWIDTH) );
     vec<graphics_primitive> x_axis = AxisX( 0, 100, 1.0, True, "%", 0.0 );
     x_axis.push_front( SetLineWidth(LINEWIDTH) );
     points.append( AxisY( 0, 100, 0, 2.3, True, "", 0.0 ) );
     vec<graphics_primitive> xlabel;
     xlabel.push_back( TextCenter( x_label_text,
          100.0/2.0, 0, 0, TimesBold(11) ) );
     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(xlabel);
     heights.push_back(65);
     stack.push_back(x_axis);    
     heights.push_back(0);
     stack.push_back(points);    
     heights.push_back(200);
     stack.push_back(title);     
     heights.push_back(20);
     RenderGraphics( file_name, stack, heights, 1.0, 200, 50, True, POSTSCALE );    
     
}     

