///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

// This template defines a class that encapsulates the generation of histograms.
// The type T must have the following operations defined:
// operator< ( const T& )
// operator+= ( const T& )
// operator*= ( const T& )
// operator* ( const T&, const T& )
//
// Note that in this implementation, neither vecData nor vecBinLimits
// needs to be sorted.

#include "CoreTools.h"
#include "PrettyPrintTable.h"
#include "graphics/BarGraph.h"
#include "math/Functions.h"

#include <map>

template <typename T>
class histogram {

 public:
  // Convenience typedef.
  typedef map<T,longlong> bin_map;
  // the only kind of iterator we support is const iterator of the bin_map
  typedef typename bin_map::const_iterator iterator;
  typedef typename bin_map::const_iterator const_iterator;

  histogram( )
    : mUnderflowCount( 0 ),
      mShowUnderflow( false ),
      mIsEmpty( true )
  { }

  /// Reset the bins' counts to zero.
  void ResetCounts(longlong value = 0)
  {
    for ( typename bin_map::iterator binIter = mBins.begin();
          binIter != mBins.end(); ++binIter )
      binIter->second = value;
    mUnderflowCount = 0;
  }

  /// Add a bin, reset counts if not empty.
  void AddBin( const T &binStart )
  {
    mBins.insert( make_pair( binStart, 0 ) );
    if ( ! mIsEmpty )
      ResetCounts();
  }

  /// Adds a set of bins starting at the given values.
  void AddBins( const vec<T> &vecBinStarts )
  {
    for ( typename vec<T>::const_iterator binStartIter = vecBinStarts.begin();
          binStartIter != vecBinStarts.end(); ++binStartIter )
      AddBin( *binStartIter );
  }

  /// This adds a set of bins whose starts follow a linear progression.

  void AddLinearProgressionOfBins( const T &lowestBinStart,
                                   const T &binSize,
                                   const int numBins )
  {
    T binStart = lowestBinStart;

    for ( int binNum = 0; binNum < numBins; ++binNum, binStart += binSize )
      AddBin( binStart );
  }

  /// This adds a set of bins whose starts follow a quasi-logarithmic
  /// progression.  No bin is added that is lower than the given lowest
  /// bin, and no bin is added that is higher than the given highest
  /// bin.  We set the current bin magnitude to the given minimum, then
  /// walk through the set of multiple progressions, adding a bin for
  /// each multiple equal to the product of the multiple and the
  /// current magnitude.  After each pass through the progression set,
  /// we multiple the current magnitude by the given bin magnitude
  /// jump.
  ///
  /// For example, to generate the logarithmic progression of bins (1,
  /// 10, 100, 1000), we would call this method with the parameters
  /// lowestBin=1, highestBin=1000, minBinMagnitude=1,
  /// binMagnitudeJump=10, and vecBinMultipleProgression=(1).
  ///
  /// To generator the quasi-logarithmic progression (10, 20, 50, 100, 200,
  /// 500, 1000), we would call this method with the parameters
  /// lowestBin=1, highestBin=100, minBinMagnitude=10,
  /// binMagnitudeJump=10, and vecBinMultipleProgression=(1,2,5).

  void AddLogProgressionOfBins( const T &lowestBin,
                                const T &highestBin,
                                const T &minBinMagnitude,  // must be > 0
                                const T &binMagnitudeJump, // must be > 1
                                const vec<T> &vecBinMultipleProgression )
  {
    ForceAssertGt( minBinMagnitude, 0 );
    ForceAssertGt( binMagnitudeJump, 1 );

    T binMagnitudeLimit = max( -lowestBin, highestBin );

    for ( T binMagnitude = minBinMagnitude;
          binMagnitude <= binMagnitudeLimit;
          binMagnitude *= binMagnitudeJump )
      for ( typename vec<T>::const_iterator multipleIter = vecBinMultipleProgression.begin();
            multipleIter != vecBinMultipleProgression.end(); ++multipleIter )
      {
        int binStart = *multipleIter * binMagnitude;
        if ( -binStart >= lowestBin )
          AddBin( -binStart );
        if ( binStart <= highestBin )
          AddBin( binStart );
      }
  }

  /// Renew this histogram from the data with default values that cover
  /// the complete range of the data.
  template <class T_iter>
  void SetAllFromData(T_iter begin, T_iter end, int nbins = 20) {
    T min = *min_element(begin, end);
    T max = *max_element(begin, end);
    T binSize = (max - min) / nbins;
    SetAllFromData(begin, end, min, binSize, nbins);
  }

  /// Renew this histogram from the data with bins set as in
  /// AddLinearProgressionOfBins
  template <class T_iter>
  void SetAllFromData(T_iter begin, T_iter end,
		      const T & min, const T & binSize, int nbins) {
    mBins.clear();
    mUnderflowCount = 0;
    mShowUnderflow = mIsEmpty = false;
    AddLinearProgressionOfBins(min, binSize, nbins);
    AddData(begin,end);
  }


  /// Get the number of bins.
  int size() const { return mBins.size(); }

  /// Get the total of all bins
  longlong CountData() const {
    longlong tot = mUnderflowCount;
    for ( typename bin_map::const_iterator binIter = mBins.begin();
          binIter != mBins.end(); ++binIter )
      tot += binIter->second;
    return tot;
  }

  /// Add a datum to the histogram, incrementing the appropriate bin
  void AddDatum( const T &datum ) { AddIdenticalData(datum,1); }

  /// Add N copies of the datum
  void AddIdenticalData( const T &datum, const longlong N )
  {
    mIsEmpty = false;

    // Find the bin after the one that should contain this datum.
    typename bin_map::iterator binIter = mBins.upper_bound( datum );

    // If we're pointing to the first bin, it's an underflow.
    if ( binIter == mBins.begin() )
      mUnderflowCount += N;

    // Otherwise, we move back one bin and increment its counter.
    else
    {
      --binIter;
      binIter->second += N;
    }
  }

  /// Get value in specified bin
  longlong operator[]( const T &datum )
  {
    typename bin_map::iterator binIter = mBins.upper_bound( datum );

    if ( binIter == mBins.begin() )
      return 0;
    else
    {
      --binIter;
      return binIter->second;
    }
  }

  /// Increment the bin counts with the given data.
  template <class T_iter>
  void AddData( T_iter dataBegin, T_iter dataEnd )
  {
    for ( ; dataBegin != dataEnd; ++dataBegin )
      this->AddDatum( *dataBegin );
  }

  /// Get the number of elements which were less than the lowest bin limit.
  longlong GetUnderflowCount() const { return mUnderflowCount; }

  /// Should we print the underflow count?
  void ShowUnderflow() { mShowUnderflow = true; }
  /// Should we print the underflow count?
  void HideUnderflow() { mShowUnderflow = false; }

  /// Get constant iterators into the bin map.
  typename bin_map::const_iterator begin() const { return mBins.begin(); }
  typename bin_map::const_iterator end()   const { return mBins.end(); }

  /// Print the histogram.
  friend
  ostream& operator<<( ostream &out, const histogram &h )
  {
    if ( h.mBins.empty() ) {
      out << "No bins in histogram.\n";
      return out;
    }
    if ( h.mShowUnderflow )
      out << "< " << h.mBins.begin()->first << "\t" << h.mUnderflowCount << "\n";
    for ( typename bin_map::const_iterator binIter = h.mBins.begin();
	  binIter != h.mBins.end(); ++binIter )
      out << binIter->first << "\t" << binIter->second << "\n";
    return out;
  }

  vec<longlong> PrintDataToVec( ) {
    vec<longlong> bins( 0 );
    if ( mBins.empty() ) return bins;
    if ( mShowUnderflow ) bins.push_back( mUnderflowCount );
    for ( typename bin_map::const_iterator binIter = mBins.begin();
	  binIter != mBins.end(); ++binIter )
      bins.push_back( binIter->second );

    return bins;
  }

  ///Print the bin starts as a line
  void PrintBinsAsLine(ostream & out) {
    if ( mBins.empty() ) {
      out << "No bins in histogram.\n";
      return;
    }

    char c[9];
    if ( mShowUnderflow )
      out << "under\t";
    for ( typename bin_map::const_iterator binIter = mBins.begin();
	  binIter != mBins.end(); ++binIter ) {
      sprintf ( c, "%-8d", static_cast<int>(binIter->first) );
      out << c;
    }
    out << "\n";
  }

  /// Print the histogram as a line.
  void PrintDataAsLine(ostream & out) {
    vec<longlong> data = PrintDataToVec( );

    if ( data.size( ) == 0 ) out << "No bins in histogram.";
    for ( int i = 0; i < data.isize( ); i++ )
    {    String s = ToString(data[i]);
         int blanks = Max( 1, 8 - s.isize( ) );
         for ( int j = 0; j < blanks; j++ )
              out << " ";
         out << s;    }
    out << "\n";
  }

  /// Print the histogram as a line, in percentages of total counts.
  void PrintDataPercentAsLine(ostream & out) {
    vec<longlong> data = PrintDataToVec( );
    char c[8];

    if ( data.size( ) == 0 ) out << "No bins in histogram.";
    double total = CountData();
    for ( unsigned int i = 0; i < data.size( ); i++ ) {
      sprintf ( c, "%.3f%%", 100.0 * data[i] / total );
      while ( strlen( c ) < 8 )	c[strlen( c )] = ' '; // pad with spaces
      out << c;
    }
    out << "\n";
  }
  
  /**
   * PrintAsColumuns
   *
   * On each line: "bin  data  cumulative<optional>  pc_total<optional>"
   *
   * pc_total: toggle to print percent of the total
   * rjust: right align columns
   * cumulative: print cumulative (if true, pc_total are cumulative too)
   */
  void PrintAsColumns(ostream &out,
		      bool pc_total = true,
		      bool rjust = false,
		      bool cumulative = true) {
    double total = CountData( );
    if ( mBins.empty() ) {
      out << "No bins in histogram." << "\n";
      return;
    }

    // Table with data.
    vec< vec<String> > table;

    // First line is the legend.
    vec<String> row;
    row.push_back( "bin" );
    row.push_back( "data" );
    if ( cumulative ) row.push_back( "cumul" );
    if ( pc_total )
      if ( cumulative ) row.push_back( "cumul %" );
      else row.push_back( "total %" );
    table.push_back( row );

    // The underflow line.
    if ( mShowUnderflow ) {
      row.clear( );
      row.push_back( "underflow" );
      row.push_back( ToString( mUnderflowCount ) );
      if ( cumulative ) row.push_back( "na" );
      if ( pc_total ) row.push_back( "na" );
      table.push_back( row );
    }

    // Generic line.
    T cumul = 0;
    for ( typename bin_map::const_iterator binIter = mBins.begin();
	  binIter != mBins.end(); ++binIter ) {
      row.clear( );
      T bin = binIter->first;
      T val = binIter->second;
      row.push_back( ToString( bin ) );
      row.push_back( ToString( val ) );
      if ( cumulative ) {
	cumul += T( val );
	row.push_back( ToString( cumul ) );
      }
      if ( pc_total ) {
	double numerator = cumulative ? cumul : val;
	double pc = 100.0 * (numerator / total);
	row.push_back( ToString( pc, 1 ) );
      }
      table.push_back( row );
    }

    // Format and print.
    int ncols = 2 + ( pc_total ? 1 : 0 ) + ( cumulative ? 1 : 0 );
    String str_rjust;
    if ( rjust )
      for (int ii=0; ii<ncols; ii++) str_rjust += "r";
    PrintTabular( out, table, 3, str_rjust );
  }
  
  /// Generate an eps plot (wrapper around BarPlot). If sep and dev
  ///  are not null, the plot will show in orange the implied normal
  ///  distribution.
  void PrintAsEps( ostream &out, 
		   const vec<ns_psplot::freetext> labels,
		   const int digits_to_right = 1,
		   const int *sep = 0,
		   const int *dev = 0 ) const {
    using namespace ns_psplot;
    
    // Define x and y bars.
    int num_bars = distance( mBins.begin( ), mBins.end ( ) );
    vec<float> x_bars;
    vec<float> y_bars;
    
    x_bars.reserve( num_bars );
    y_bars.reserve( num_bars );
    typename bin_map::const_iterator binIter;    
    for ( binIter = mBins.begin(); binIter != mBins.end(); ++binIter ) {
      x_bars.push_back( float( binIter->first ) );
      y_bars.push_back( float( binIter->second ) );
    }
    
    // Define y2 (the normal distribution implied by sep and dev).
    vec<float> y2_bars;
    color col2 = orange;
    if ( sep && dev ) {
      float fsep = float ( *sep );
      float fdev = float ( *dev );
      float deltax = x_bars[1] - x_bars[0];
      NormalDistribution nd( fsep, fdev );

      float totalsum = 0;
      for (int ii=0; ii<y_bars.isize( ); ii++)
	totalsum += y_bars[ii];
      for (int ii=0; ii<y_bars.isize( ); ii++)
	y_bars[ii] = 100.0 * ( y_bars[ii] / totalsum );
    
      y2_bars.reserve( num_bars );
      for ( const_iterator it=mBins.begin( ); it!=mBins.end( ); it++ ) {
	float x2 = float( it->first );
	float y2 = 100.0 * deltax * nd.ProbabilityDensity( x2 );
	y2_bars.push_back( y2 );
      }
    }
    
    // Generate figure.
    using namespace ns_psplot;
    
    float max_bar_height = 150;   // hight (in pixels) of the tallest bar
    float bar_width = 6;          // width of each bar (in pixels)
    float bar_sep = 1.5;          // distance between bars (in pixels)
    color col = gray;             // color for the bars (see BarGraph.h)
    short font_size = 7;          // font size for the tags (tic labels)
    
    float low = x_bars[0];
    float high = x_bars[x_bars.size()-1];
    vec<float> bot_tics;
    for ( typename bin_map::const_iterator binIter = mBins.begin();
	  binIter != mBins.end();
	  ++binIter )
      bot_tics.push_back( float( binIter->first ) );
    
    bargraph graph
      = ( sep && dev )
      ? BiBarPlot( y_bars, y2_bars, low, high, bar_width, bar_sep, col,
		   col2, font_size, max_bar_height, labels, true )
      : BarPlot( y_bars, low, high, bar_width, bar_sep, col, font_size,
		 max_bar_height, labels, bot_tics, true, digits_to_right );
    out << graph << endl;
  }
  
  
 private:
  bin_map mBins;
  longlong mUnderflowCount;
  bool mShowUnderflow;
  bool mIsEmpty;
};

template <class T>
inline void PrintSimpleHistogram( const vec< T >& v ) {
  histogram< T > h;
  h.SetAllFromData( v.begin(), v.end() );
  cout << h << endl;
}

#endif
