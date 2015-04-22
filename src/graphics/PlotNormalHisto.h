// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology



/*
 * PlotNormalHisto
 *
 * Functions to plot very basic (normal) histogram graphs.
 */
#ifndef PLOT_NORMAL_HISTO_H
#define PLOT_NORMAL_HISTO_H



/*
 * RawPlot
 * 
 * This is the core plot function.
 *
 * max_bars: max number of bars to be plotted
 * plot_implied: if true will bi-plot also the normal N(mu, sigma)
 */
void RawPlot( const String &eps_file,
	      const String &title,
	      const vec<int> data,
	      float mean,
	      float stdev,
	      int max_bars = 50,
	      bool plot_implied = true );



/*
 * BasicNormalHistoPlot
 * 
 * It bar-plots a set of n points (i.e. an histogram where the x-bars are
 * the points [0, 1, ..., n-1]), with the option of plotting also the
 * density function of the normal distribution implied by the given data.
 * It generates an eps file.
 */
void BasicNormalHistoPlot( const String &eps_file,
			   const String &title,
			   const vec<int> data,
			   int max_bars = 50,
			   bool plot_implied = true );



/*
 * FatTailHistoPlot
 * 
 * Pretty much as in BasicNormalHistoPlot, but with a twist. The idea
 * is that data is believed to have a fat tail, so that mean and standard
 * deviation are computed by only looking at "half" the data. Notice that
 * this will work only if there are enough data points, since mean and
 * standard deviation are computed on the subset data[1], ...., data[n],
 * where: (a) data[0] has obviously been excluded; (b) data[1] < data[2]
 * data[3] < ... < data[n]; and (c) data[n] > data[n+1]. (These points
 * are then mirror-imaged to get the remaining data-ponts). max_bars
 * sets a limit to the number of bars to be plotted (data may be very large).
 *
 * It returns false on failure.
 */
bool FatTailHistoPlot( const String &eps_file,
		       const String &title,
		       const vec<int> data,
		       int max_bars = 50,
		       bool plot_implied = true );



#endif
