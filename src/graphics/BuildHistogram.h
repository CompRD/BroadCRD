// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef BUILD_HISTOGRAM_H
#define BUILD_HISTOGRAM_H

#include "Vec.h"



/*
 * Class histogram
 *
 * Utility class for preparing and handling data for histograms of positive floating point numbers.
 * The responsibility of this class is to take in the list values from which a histogram is to be 
 * created and to generate breaking points and heights of the histogram's bars. First configure the
 * histogram using Set... methods, then call ProduceHistogram() to generate bars.
 */
class float_histogram {

public:

  float_histogram( vec<float>& lengths );

  // set the desired number of bars in the histogram (optional,
  // a default value is set).
  void SetBarNumber( int num_bars ) { 
    int min_bars = 3;
    int max_bars = 100;
    if ( num_bars < min_bars || num_bars > max_bars ) {
      cout << "Allowed number of bars is between " << min_bars << " and " << max_bars << endl;
      cout << "Request to set number of bars to " << num_bars << " was ignored" << endl;
      return;
    }
    num_bars_ = num_bars;
    AdjustBarWidth();
  }


  // set the limits for the region of data values to create the histogram
  // bars for; by default, the values falling out of the specified region are
  // counted into the leftmost and
  // rightmost bars, respectively (see SetCountOutliers()). 
  // If limits are not set explicitly with this method,
  // they will be set automatically to the whole range spanned by the data values.
  void SetXLimits(float xmin, float xmax) { 
    xmin_ = xmin; 
    xmax_ = xmax; 
    manual_limits_ = true;
    AdjustBarWidth();
  };

  // If set to \c true (default), then data values below and above of the plotting range are
  // counted into the leftmost and rightmost bars, respectively; if set to \c false, then these
  // out-of-range values are simply discarded.
  void SetCountOutliers(bool count_outliers) { count_outliers_ = count_outliers; }

  // If you want to use your own x_bars, then pass a nonempty vector,
  //  and set generate_x_bars = false (otherwise x_bars are generated).
  int ProduceHistogram( vec<float> &x_bars, vec<float> &y_bars,
			bool generate_x_bars = true );
  
  
private:
  
  // computes bar width from the data
  // (hence used when bars are generated automatically);
  // may also adjust xmin down a bit in order to have
  // more naturally looking bar boundaries
  void AdjustBarWidth() ;

  int num_bars_;          // number of histogram bars;
  float xmin_, xmax_;
  float bar_width_;
  bool count_outliers_ ;
  bool manual_limits_; ///< state variable; used to remember if limits were set manually; 
                       /// affects behavior of AdjustBarWidth()
  vec<float>& lengths_;   // data; will be kept sorted in descending order

};



#endif
