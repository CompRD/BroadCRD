/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_PLOT_COVERAGE_H
#define C_PLOT_COVERAGE_H

#include "Qualvector.h"
#include "SeqInterval.h"
#include "String.h"
#include "Vec.h"

/**
 * Class CPlotCoverage
 *
 * Tool to plot various coverage related figures. Given a reference
 * genome (or assembly), and given any kind of coverage for the same,
 * it will plot one figure per chromosome (or contig), and one overall
 * histogram for the distribution of coverages on the sliding windows.
 */
class CPlotCoverage {

public:

  CPlotCoverage(int win, float cap, const String &name, const String &logfile);

  ~CPlotCoverage( );
  
  void SetWinSize( int window_size ) { win_size_ = window_size; }

  void SetMult( float mult ) { mult_ = mult; }

  // Use output from CoverageAnalyzer::GetAllCoverages( ) to fill cov_.
  void SetFromSeqIntervals( const vec<int> &tlens, 
			    const vec<int> &fcov,
			    const vec<seq_interval> &coverages );
  
  // Use output from GenerateCoverageVecvec to fill cov_.
  void SetFromVecvec( const vecqualvector &vcov );
  
  // Plot only overall histo if histo_only = true (out_dir must exist)
  void PlotTargets( const String &out_dir, bool histo_only );

  int WinSize( ) const { return win_size_; }

  float Mult( ) const { return mult_; }

  String Name( ) const { return name_; }
  
  
private:

  // The histogram of averages is saved on histo_file.
  float GlobalMode( const String &histo_file );

  void ResizeAverages( const vec<int> &tlens );

  void GenerateAverages( int target_id, const vec<float> &raw );
  

private:

  int win_size_;     // size of sliding (hopping) window
  int hop_size_;     // hop amount
  float cap_;        // maximum coverage in the histogram plots
  float mult_;       // multiplicative factor, to tag high coverage windows
  ofstream log_;     // output stream
  String name_;      // descriptor used to decorate plots

  float mode_;                  // global mode (max of distribution)
  vec<int> tlens_;              // target lengths
  vec< vec<float> > averages_;  // averages on windows (not sorted)
  
};

#endif
