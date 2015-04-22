// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "graphics/PartitionedLineGraph.h"
#include <strstream>
#include <numeric>

using namespace ns_psplot;

class rlineto {
     public: rlineto( float x, float y ) : x_(x), y_(y) { }
     friend ostream& operator<<( ostream& o, const rlineto& r )
     {    return o << " " << r.x_ << " " << r.y_ << " rlineto";    }
     private: float x_, y_;
};

class lineto {
     public: lineto( float x, float y ) : x_(x), y_(y) { }
     friend ostream& operator<<( ostream& o, const lineto& r )
     {    return o << " " << r.x_ << " " << r.y_ << " lineto";    }
     private: float x_, y_;
};

class translate {
     public: translate( float x, float y ) : x_(x), y_(y) { }
     friend ostream& operator<<( ostream& o, const translate& t )
     {    return o << " " << t.x_ << " " << t.y_ << " translate";    }
     private: float x_, y_;
};

class startpath {
     public: startpath( float x, float y ) : x_(x), y_(y) { }
     friend ostream& operator<<( ostream& o, const startpath& s )
     {    return o << " newpath " << s.x_ << " " << s.y_ << " moveto";    }
     private: float x_, y_;
};


// Partitioned line graph

partitioned_line_graph::partitioned_line_graph( const vec<float>& chunks, const float total, 
						const int height, const int width,
						const vec<freetext> labels,
						const vec<unsigned int>& top_Ns,
						const color unused_portion_color,
						const color frame_color )
  : height_(height), width_(width),
    labels_(labels),
    top_Ns_(top_Ns),
    NXs_(9),
    unused_portion_col_(unused_portion_color), 
    frame_col_(frame_color) 
{
  partial_sum( chunks.begin(), chunks.end(),
	       back_inserter( slices_ ) );
  
  int slice_idx = 0;
  for ( int N = 0; N < 9; ++N )
  {
    float target = (float)(N+1)/10.0 * slices_.back();
    
    while ( slices_[slice_idx] < target )
      ++slice_idx;

    NXs_[ N ] = chunks[slice_idx];
  }

  float pixels_per_unit = (float) width / total;
  transform( slices_.begin(), slices_.end(),
	     slices_.begin(), 
	     bind1st(multiplies<float>(), pixels_per_unit) );

  float too_small_to_see = total / (float) width * 0.25;
  vec<float>::const_iterator first_invisible_chunk 
    = lower_bound( chunks.begin(), chunks.end(),
		   too_small_to_see,
		   greater<float>() );
  
  first_invisible_slice_index_ = distance( chunks.begin(), first_invisible_chunk );
}


ostream& ns_psplot::operator<<( ostream& o, const partitioned_line_graph& plg )
{

  const int legend_fontsize = 6;

  vec<short int> label_heights;
  transform( plg.labels_.begin(), plg.labels_.end(),
	     back_inserter(label_heights),
	     mem_fun_ref( &freetext::FontSize ) );

  short int total_label_height = 
    accumulate( label_heights.begin(), label_heights.end(), 0 ) + LabelSep * plg.labels_.size();

  int total_top_n_bracket_height = plg.top_Ns_.size() * (legend_fontsize+2);

  const int width = plg.width_;
  const int height = plg.height_;
  const int left_horiz_buffer = 75;
  const int right_horiz_buffer = 75;
  const int top_vert_buffer = 20 + total_label_height;
  const int bottom_vert_buffer = 20 + total_top_n_bracket_height;

  o << "%!PS-Adobe-3.0 EPSF-3.0\n"
    << "%%BoundingBox: "
    << "0 0 "
    << left_horiz_buffer + width + right_horiz_buffer << " " 
    << bottom_vert_buffer + height + top_vert_buffer << endl;

  o << left_horiz_buffer << " " << bottom_vert_buffer << " translate" << endl;

  // Make sure we always have at least two decimal places.
  int pixel_precision = 2 + (int) ceilf(log10(max<float>(height,width)));
  int old_precision = o.precision( pixel_precision );

  // Set line width.
  o << "0.6 setlinewidth" << endl;

  // Draw frame.
  o << plg.frame_col_
    << startpath( 0, height/2 + 0.3) 
    << lineto( width, height/2 + 0.3)
    << " stroke" << endl;

  // Draw unused portion.
  float last_slice = plg.slices_.back();

  o << startpath( last_slice , height/2 + 0.3) 
    << lineto( width, height/2 + 0.3)
    << plg.unused_portion_col_ << " stroke " << plg.frame_col_ << endl;

  // Set line width.
  o << "0.2 setlinewidth" << endl;

  // Draw caps.
  o << startpath( 0, 0 ) 
    << lineto( 0, height )
    << " stroke" << endl;

  o << startpath( width, 0 ) 
    << lineto( width, height )
    << " stroke" << endl;

  // Draw N* hashes
  for ( unsigned int N = 0; N < 9; ++N )
  {
    float offset = float(N+1)/10 * last_slice;
    o << "newpath" << endl;
    o << offset << " " << 0 << " moveto" << endl;
    o << 0 << " " << legend_fontsize*2+2 << " rlineto" << endl;
    o << "0.5 setgray" << endl;
    o << "stroke" << endl;
    o << "0.0 setgray" << endl;
    o << "gsave" << endl;
    o << offset-1 << " " << legend_fontsize+2 << " translate" << endl;
    o << text( freetext( "N" + ToString( (N+1)*10 ), black, legend_fontsize ), right ) << endl;
    o << 2 << " " << 0 << " translate" << endl;

    float value = plg.NXs_[N];
    String unit;
    if ( value > 1000.0 )
    {
      value /= 1000.0;
      unit = " Kbp";
    }
    if ( value > 1000.0 )
    {
      value /= 1000.0;
      unit = " Mbp";
    }
    if ( value > 1000.0 )
    {
      value /= 1000.0;
      unit = " Gbp";
    }

    int digits;
    if ( unit.empty() )
    {
      unit = " bp";
      digits = 0;
    }
    else
    {
      if ( value > 100.0 )
        digits = 0;
      else if ( value > 10.0 )
        digits = 1;
      else 
        digits = 2;
    }

    o << text( freetext( ToString( value, digits ) + unit, black, legend_fontsize ), left ) << endl;
    o << "grestore" << endl;
  }

  // Draw visible slices.
  for ( unsigned int i = 0; i < plg.first_invisible_slice_index_; i++ )
  {
    o << startpath( plg.slices_[i], 0 )
      << rlineto( 0, height ) << " stroke" << endl;
  }

  // Draw block representing invisible slices.
  float last_visible_slice;
  if ( plg.first_invisible_slice_index_ == 0 )
    last_visible_slice = 0;
  else
    last_visible_slice = plg.slices_[plg.first_invisible_slice_index_ - 1];
  
  o << startpath( last_visible_slice, 0 )
    << lineto( last_slice, 0 )
    << lineto( last_slice, height )
    << lineto( last_visible_slice, height )
    << lineto( last_visible_slice, 0 )
    << " closepath" << endl;
    
  o << " fill" << endl;

  // Calculate and indicate n largest, where n is 10, 25, 50, 75, 100, 250...

  o << "/bracket  % params: points_up vertical_offset vertical_size horizontal_size" << endl
    << "{ /width exch def   % define width" << endl
    << "  /height exch def  % define height" << endl
    << "  /offset exch def  % define vertical offset" << endl
    << "  { /abs_or_neg { neg } def } % define abs_or_neg as abs if points_up is true" << endl
    << "  { /abs_or_neg { abs } def } % define abs_or_neg as neg if points_up is false" << endl
    << "  ifelse" << endl
    << "  newpath                            % start bracket" << endl
    << "  0 offset 1 add abs_or_neg moveto   % move to 0, +-(offset+1)" << endl
    << "  0 height abs_or_neg rlineto        % draw a line up or down height pixels" << endl
    << "  width 0 rlineto                    % draw a line over length pixels" << endl
    << "  0 height abs_or_neg neg rlineto    % draw a line back down or up" << endl
    << "  0 setgray stroke                   % stroke the bracket in black" << endl
    << "  newpath                            % add tracking lines to bracket" << endl
    << "  0 1 abs_or_neg moveto              % move to 0, +-1" << endl
    << "  0 offset abs_or_neg rlineto  % draw a line to left end of bracket" << endl
    << "  width 1 abs_or_neg moveto      % move to width, +-1" << endl
    << "  0 offset abs_or_neg rlineto  % draw a line to right end of bracket" << endl
    << "  0.65 setgray stroke                % draw light gray tracking lines" << endl
    << "  0 setgray                % reset color to black" << endl
    << "} def" << endl;

  o << "/namedbracket  % params: on_bottom vertical_offset vertical_size" << endl
    << "               %         horizontal_length name percentage" << endl
    << "{ /percent exch def        % define percentage" << endl
    << "  /name exch def           % define name" << endl
    << "  /length exch def         % define horizontal_length" << endl
    << "  /size exch def           % define vertical_size" << endl
    << "  /offset exch def         % define vertical_offset" << endl
    << "  /on_bottom exch def      % define on bottom" << endl
    << "  on_bottom                % if on_bottom" << endl
    << "  { /abs_or_neg { neg } def }     % define abs_or_neg as neg if on_bottom is true" << endl
    << "  { 0 " << height << " translate  % translate above graph and" << endl
    << "    /abs_or_neg { abs } def }     % define abs_or_neg as abs if on_bottom is false" << endl
    << "  ifelse" << endl
    << "  on_bottom offset size length bracket  % draw bracket" << endl
    << "  newpath                  % start a new path" << endl
    << "  name stringwidth pop 2 add neg  % get -(string width + 2)" << endl
    << "  offset abs_or_neg 2 add         % get 2+-(offset)" << endl
    << "  on_bottom { size 2 add sub } if   % subtract size+2 if on bottom" << endl
    << "  moveto                   % move there" << endl
    << "  name show                % show the name" << endl
    << "  length 4 add 0 rmoveto   % move to end of bracket" << endl
    << "  percent show             % show the percentage" << endl
    << "  stroke         % stroke the name in black" << endl
    << "  on_bottom not            % if not on_bottom" << endl
    << "  { 0 " << height << " neg translate }  % translate back below graph" << endl
    << "  if" << endl
    << "} def" << endl;

  o << "/Times-Roman findfont " << legend_fontsize << " scalefont setfont" << endl;

  string prefix;
  unsigned int num_Ns = plg.top_Ns_.size();
  for ( unsigned int offset = 0; offset < num_Ns; offset++ )
  {
    if ( offset == num_Ns-1 )
      prefix = "longest ";

    unsigned int n = plg.top_Ns_[offset];
    o << "true " // on_bottom
      << (num_Ns - offset - 1) * (legend_fontsize+2) << " "
      << legend_fontsize+1 << " "
      << plg.slices_[n-1] 
      << " (" << prefix << n << ") " 
      << setprecision( 3 )
      << " (" << plg.slices_[n-1] / last_slice * 100 << "%)"
      << setprecision( pixel_precision )
      << " namedbracket" << endl;
  }    

  // Labels.
  o << translate( width/2, legend_fontsize+2 + 2*LabelSep );
  for ( int i = plg.labels_.size( ) - 1; i >= 0; i-- )
    o << translate( 0, LabelSep ) << text( plg.labels_[i], center )
      << translate( 0, plg.labels_[i].FontSize( ) );
  
  // Reset precision.
  o.precision( old_precision );

  // Return.
  return o;
}
  

partitioned_line_graph ns_psplot::PartitionedLinePlot( const vec<float>& chunks, 
						       const float total, 
						       const int height, const int width,
						       const vec<freetext> labels,
						       const Bool display_top_n_chunks,
						       const color unused_portion_color,
						       const color frame_color )
{
  float largest( chunks.front() );
  float smallest( chunks.back() );
  if ( smallest < 2000 )
    smallest = 2000;

  vec<unsigned int> top_Ns;
  if ( display_top_n_chunks ) 
  {
    unsigned int n = 10;
    unsigned int n_magnitude = 10;
 
    unsigned int num_chunks = chunks.size();
    while ( n <= num_chunks )
    {
      top_Ns.push_back(n);
      
      if ( n == n_magnitude)
      {
        n_magnitude *= 10;
        n = n_magnitude / 4;
      }
      else
        n *= 2;
    }
  }

  return partitioned_line_graph( chunks, total,
				 height, width,
				 labels,
				 top_Ns, 
				 unused_portion_color, frame_color );
}
    
  
  

