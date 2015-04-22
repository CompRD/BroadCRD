/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This file contains tools for generating bargraphs.  It produces postscript
// output.  One way to view them is to run ps2pdf on the output and then view the
// pdf file with acroread.
//
// An instance of class bargraph is constructed using
//
// bargraph( vec<bar> bars, vec<freetext> labels, vec<freetic> bottom_tics,
//           vec<freetic> left_tics, bool use_eps )
//
// where:
//
// bar( x, y, width, color ) describes a bar
//
// freetext( string, color, fontsize ) describes alabel
//
// freetic( t, string, color, fontsize, angle ) 
//          describes a labeled tic mark at position t
// (angles not well implemented)
//
// Labels go across the top of the graph, centered, one line per entry in "labels".
//
// If you use the ostream operator on a bargraph you get postscript output.

#ifndef BARGRAPH
#define BARGRAPH

#include "CoreTools.h"
#include "graphics/PostScript.h"
#include "math/Functions.h"

// Remark: the tools now produce eps bar and function graphs (with some
// limitations.)
//
// Namespace:
//  ns_psplot.
//
// Principal classes:
//  bargraph (used by BarPlot and BiBarPlot;)
//  plotgraph (used by GraphPlot.)
//
// Principal functions:
//  BarPlot;
//  BiBarPlot;
//  GraphPlot.

namespace ns_psplot {

const int TicLength = 5, TicLabelSep = 3;

class freetic {

     public:

     freetic( float t, freetext label, float angle = 0 ) 
          : t_(t), label_(label), angle_(angle) { }

     freetic( float t, string str, color col, short font_size, float angle = 0 )
          : t_(t), label_(str, col, font_size), angle_(angle) { }

     freetic( ) { }

     float T( ) const { return t_; }
     freetext Label( ) const { return label_; }
     float Angle( ) const { return angle_; }

     private:

     float t_;
     freetext label_;
     float angle_;

};

class tic : public freetic {

     public:

     tic( float t, place dir, freetext label ) : freetic(t, label), dir_(dir)
     {    Assert( dir == left || dir == right );    }

     tic( freetic f, place dir ) : freetic(f), dir_(dir)
     {    Assert( dir == left || dir == right );    }

     tic( ) { }

     friend ostream& operator<<( ostream& o, const tic& t );

     private:

     float t_;
     place dir_;
     freetext label_;

};

const float EdgeFudge = 0.5;

class bar {

     public:

     bar( float x, float y, float width, color c )
          : x_(x), y_(y), width_(width), c_(c) { }

     bar( ) { }

     float X( ) const { return x_; }
     float Y( ) const { return y_; }
     float Width( ) const { return width_; }

     friend ostream& operator<<( ostream& o, const bar& b );

     private:

     float x_, y_, width_;
     color c_;

};

const int LabelSep = 5, Margin = 8;

class bargraph {

     public:
  
     bargraph( vec<bar> bars, vec<freetext> labels, vec<freetic> bottom_tics,
	       vec<freetic> left_tics, bool use_eps ) :
       bars_(bars),
       labels_(labels),
       bottom_tics_(bottom_tics),
       left_tics_(left_tics),
       use_eps_(use_eps)
     {    Assert( bars.size( ) > 0 );    }

     friend ostream& operator<<( ostream& o, const bargraph& g );

     private:

     float LeftBoundary( ) const
     {    float m = 1000000;
          for ( unsigned int i = 0; i < bars_.size( ); i++ )
               m = Min( m, bars_[i].X( ) - bars_[i].Width( )/2 );
          return m - Margin;    }

     float RightBoundary( ) const
       {
	 float ret_value = 0;
	 
	 float M = -1000000;
	 for ( unsigned int i = 0; i < bars_.size( ); i++ )
	   M = Max( M, bars_[i].X( ) + bars_[i].Width( )/2 );
	 
	 // Approximate length of labels (roughly speaking, two characters
	 // are ~ FontSize(), in length.)
	 int label_length = 0;
	 for (int ii=0; ii<(int)labels_.size(); ++ii)
	   {
	     int n_char = labels_[ii].Str().size();
	     int loc_length = ( n_char/2 ) * labels_[ii].FontSize( );
	     label_length = (loc_length > label_length ) ? loc_length : label_length;
	   }

	 ret_value = ( label_length > M + Margin ) ? label_length : M + Margin;
	 ret_value = float( int(ret_value) + 1 );

	 return ret_value;
       }

     float Width( ) const { return RightBoundary( ) - LeftBoundary( ); }

     float BarTopBoundary( ) const
     {    float M = -1000000;
          for ( unsigned int i = 0; i < bars_.size( ); i++ )
               M = Max( M, bars_[i].Y( ) );
          return M;    }

     float Height( ) const
     {    float M = BarTopBoundary( );
          for ( unsigned int i = 0; i < labels_.size( ); i++ )
               M += labels_[i].FontSize( ) + LabelSep;
          return M + Margin;    }

     vec<bar> bars_;
     vec<freetext> labels_;
     vec<freetic> bottom_tics_, left_tics_;
     bool use_eps_;
};

class plotgraph {
 public:
  
     plotgraph(vec<segment> segments,
	       vec<freetext> labels,
	       vec<freetic> bottom_tics,
	       vec<freetic> left_tics,
	       float pixels_per_unit_x,
	       float pixels_per_unit_y,
	       float max_y,
	       bool plot_dots) :
       segments_(segments),
       labels_(labels),
       bottom_tics_(bottom_tics),
       left_tics_(left_tics),
       pixels_per_unit_x_(pixels_per_unit_x),
       pixels_per_unit_y_(pixels_per_unit_y),
       max_y_(max_y),
       plot_dots_(plot_dots) { Assert( segments.size( ) > 0 ); }
     
     friend ostream& operator<<( ostream& o, const plotgraph& plotg );

 private:
     float LeftBoundary( ) const { return segments_[0].X1(); }
     
     float RightBoundary( ) const
       {
	 float ret_value = 0;

	 ret_value = segments_[segments_.size()-1].X2() + pixels_per_unit_x_;

	 // Approximate length of labels (roughly speaking, two characters
	 // are ~ FontSize(), in length.)
	 int label_length = 0;
	 for (int ii=0; ii<(int)labels_.size(); ++ii)
	   {
	     int n_char = labels_[ii].Str().size();
	     int loc_length = ( n_char/2 ) * labels_[ii].FontSize( );
	     label_length = (loc_length > label_length ) ? loc_length : label_length;
	   }

	 ret_value = ( label_length > ret_value ) ? label_length : ret_value;

	 return ret_value;
       }

     float Width( ) const { return RightBoundary( ) - LeftBoundary( ); }

     float PlotTopBoundary( ) const
       {
	 return max_y_*pixels_per_unit_y_;
       }

     float Height( ) const
       {
	 float M = PlotTopBoundary( );
	 for ( unsigned int i = 0; i < labels_.size( ); i++ )
	   M += labels_[i].FontSize( ) + LabelSep;
	 return M + Margin;
       }

     vec<segment> segments_;
     vec<freetext> labels_;
     vec<freetic> bottom_tics_;
     vec<freetic> left_tics_;
     float pixels_per_unit_x_;
     float pixels_per_unit_y_;
     float max_y_;
     bool plot_dots_;
};

bargraph BarPlot( vec<float> y, float low, float high, float bar_width, 
		  float bar_sep, color c, short font_size,
		  float max_bar_height, vec<freetext> labels,
		  vec<float> f_bottom_tics = vec<float>(0), bool use_eps = false,
                  int digits_to_right = 4, bool bottom_tics_horizontal = false );

bargraph BiBarPlot( vec<float> y1, vec<float> y2, float low, float high, 
     float bar_width, float bar_sep, color c1, color c2, short font_size, 
     float max_bar_height, vec<freetext> labels, bool use_eps = false );

plotgraph GraphPlot( vec<float>& x, vec<float>&y, color col,
		     short font_size, vec<freetext> labels,
		     int graph_width, int graph_height,
		     bool plot_dots = false );

}

#endif
