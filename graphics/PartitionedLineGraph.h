// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef PARTITIONEDLINEGRAPH
#define PARTITIONEDLINEGRAPH

#include "graphics/BarGraph.h"


namespace ns_psplot {

class partitioned_line_graph 
{
public:
  partitioned_line_graph( const vec<float>& chunks, const float total, 
			  const int height, const int width,
			  const vec<freetext> labels,
			  const vec<unsigned int>& top_Ns,
			  const color unused_portion_color = red,
			  const color frame_color = black );

  friend ostream& operator<<( ostream& o, const partitioned_line_graph& cbg );

private:

  float LeftBoundary( ) const
    { 
      return 0;    
    }

  float RightBoundary( ) const
    {
      return width_;
    }
  
  float Width( ) const { return RightBoundary( ) - LeftBoundary( ); }
  
  float Height( ) const
    {
      return height_;    
    }

  int height_, width_;
  vec<float> slices_;
  unsigned int first_invisible_slice_index_;
  vec<freetext> labels_;
  vec<unsigned int> top_Ns_;
  vec<float> NXs_;
  float used_portion_size_;
  color unused_portion_col_;
  color frame_col_;
};

// PartitionedLinePlot requires that the vec<int> chunks be in sorted order
// from smallest to largest, e.g. sort( chunks.rbegin(), chunks.rend() )

partitioned_line_graph PartitionedLinePlot( const vec<float>& chunks, const float total, 
					    const int height, const int width,
					    const vec<freetext> labels,
					    const Bool display_top_n_chunks,
					    const color unused_portion_color = red,
					    const color frame_color = black );

}

#endif
