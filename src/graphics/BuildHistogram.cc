// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include <algorithm>
#include <math.h>

#include "Vec.h"
#include "graphics/BuildHistogram.h"



/*
 * class: float_histogram
 * constructor
 */
float_histogram::float_histogram( vec<float>& lengths ) :
  num_bars_( 35 ),
  count_outliers_(true),
  manual_limits_(false),
  lengths_( lengths )
{
  sort(lengths_.begin(), lengths_.end(), greater<float>() ); // sort in descending order
    
    // min and max values
  if ( lengths_.size() > 0 ) {
    xmax_ = lengths_.front();
    xmin_ = lengths_.back();
  } else {
    xmax_ = xmin_ = 0;
  }
  AdjustBarWidth();
}


void float_histogram::AdjustBarWidth() {
  
    bar_width_ = (xmax_ - xmin_)/static_cast<float>(num_bars_);
    PRINT(bar_width_);
    if ( manual_limits_ ) { 
      // don't adjust limits if they were set manually, just stay with raw bar_width whatever it turns out to be
      return;
    }

    // if automatic limits are used, we will try to compute some reasonable
    // and nicely looking bar_width:
    
    // Check data. Since lengths_ is sorted, the last entry is the
    // shortest length.
    if ( 0 == lengths_.size() ) return;
    
    // "raw" bar width can be any weird number, let's adjust barwidth up towards
    //  the closest of (1,2,3,4,5,10)*1eN, where N is an arbitrary power.
    float sep_log = log10f(bar_width_);
    float sep_power = floorf(sep_log);
    float sep_base = bar_width_ / powf(10.0f, sep_power);

    if ( sep_base <= 5.0f ) bar_width_ = ceilf(sep_base)*powf(10.0f,sep_power);
    else bar_width_ = powf(10.0f, sep_power+1);
    
    // adjust xmin_ down to a multiple of bar_width
    // (we do not want, e.g.,  xmin=0.1123 if "thick" bars of width 0.1 are used,
    // it would be much more natural to use 0.1-0.2-0.3-... instead!)
    xmin_ = floorf(xmin_/bar_width_)*bar_width_; 

    if ( xmin_ + bar_width_*static_cast<float>(num_bars_) > xmax_ ) {
      // ooops, our adjustment of xmin_ resulted in some values falling out on the right!
      // recompute/readjust bar_width (another iteration)
      bar_width_ = (xmax_-xmin_)/static_cast<float>(num_bars_);

      sep_log = log10f(bar_width_);
      sep_power = floorf(sep_log);
      sep_base = bar_width_ / powf(10.0f, sep_power);
      if ( sep_base <= 5.0f ) bar_width_ = ceilf(sep_base)*powf(10.0f,sep_power);
      else bar_width_ = powf(10.0f, sep_power+1);
    }

    PRINT3(xmin_, xmax_, bar_width_);
}


/*
 * class: float_histogram
 * ProduceHistogram
 *
 * Produces the histogram data. y_bars[ii] is the number of entries in
 * lengths_ which are in the interval ( x_bars[ii-1], x_bars[ii] ],
 * for ii>1, and in the interval [0, x[ii] ], for ii=0. Remark, lengths_
 * must be sorted, and non negative.
 *
 * Legenda:
 *  x_bars: where x coordinates will be stored (output;)
 *  y_bars: where y coordinates will be stored (output;)
 *
 * Return:
 *  1 for success, 0 for failure.
 */
int float_histogram::ProduceHistogram( vec<float> &x_bars, vec<float> &y_bars, 
                                       bool generate_x_bars )
{
  if ( !generate_x_bars && x_bars.size() < 1 ) {
    cout << "Use of external breakpoints is requested, but empty vector passed" << endl;
    ForceAssert( 1 == 0 );
  }

  // Check if x_bars is given in input or must be generated.
  if ( !generate_x_bars ) {
    num_bars_ = x_bars.size();
    // this is ugly, but conformant to the historic interface: xbars contain only left
    // side of each requested bar, so that we have no idea where the rightmost bar ends
    // (i.e. where its right side is located). Hence we can not discard the "outliers"
    if ( ! count_outliers_ ) {
      cout << "Can not ignore outliers when use of external breakpoints is requested" << endl;
      return 0;
    }
  }
  else {
    PRINT(num_bars_);
    // set x_bars:
    x_bars.resize(num_bars_);
    x_bars[0] = xmin_;
    for (int ii=1; ii<num_bars_; ++ii)
      x_bars[ii] = x_bars[ii-1] + bar_width_;
  }    

  // y_bars.
  y_bars.clear();
  y_bars.resize(num_bars_, 0);
  
  // now let's count values into appropriate bars
  // (remember that values are sorted in descending order!):
  
  int jj_level = num_bars_ - 1;
  int ii = 0;
  if ( ! count_outliers_ ) {
    // skip all the values that are greater than the rightmost bar: 
    for ( ; lengths_[ii] > x_bars[jj_level]+bar_width_; ii++ );
  }

  for ( ; ii<(int)lengths_.size(); ++ii) {
    if ( jj_level > 0 ) {
      while ( lengths_[ii] < x_bars[jj_level] ) {
	--jj_level;	
	if ( 0 == jj_level )
	  break;
      }
    }
    // if we discard outliers and current value is below lowest bar - it's time to stop counting
    if ( ! count_outliers_ && lengths_[ii] < x_bars[0] ) break; 
    ++y_bars[jj_level];
  }
  
  // All right.
  return 1;
}
  
