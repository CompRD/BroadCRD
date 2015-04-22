// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef COVERAGE_ANALYZER_H
#define COVERAGE_ANALYZER_H

#include <fstream>
#include <list>
#include <vector>

#include "Vec.h"
#include "SeqInterval.h"
#include "system/Types.h"

using namespace std;

/**
 * \class CoverageAnalyzer
 *
 * Analyzes the coverage of intervals in sequences (e.g. the read coverage
 * of a contig, the read coverage of a supercontig, the insert coverage
 * of a supercontig).
 *
 * Intervals are (simply connected) subsets of genomic sequences. Think of
 * this example: intervals=reads, sequences=contigs (when looking for the 
 * read coverage of the contigs).
 *
 */

class CoverageAnalyzer
{
public:
  /// create an empty analyzer.  If a blank one is created, then you'll
  /// have to call CreateCoverages() after creation.  If you specify an
  /// ostream, then any error messages will be printed out via this stream.
  CoverageAnalyzer(ostream* pOStream = NULL);

  /// create a analyzer and have it create coverages for the intervals
  /// specified in vecIntervals.  The constructor simply calls
  /// CreateCoverages() below, so see that for usuage.  If you specify an
  /// ostream, then any error messages will be printed out via this stream.
  CoverageAnalyzer(const vec<seq_interval>& vecIntervals,
		   const vec<int>* pVecSeqLengths = NULL,
		   ostream* pOStream = NULL);

  ~CoverageAnalyzer() {
    if ( m_pSeqLengths )
      delete m_pSeqLengths;
  }

  /// The first thing CreateCoverages() does is to delete coverages that
  /// it created previously.  If you want to keep those coverages, then
  /// you should probably create another coverage analyzer object.
  ///
  /// CreateCoverages() creates coverages for all the intervals specified
  /// in vecIntervals.  If you do not specify a set of seq lengths, then
  /// the coverage regions will have no bounds.  If you specify the seq
  /// lengths, any interval that is outside of the area defined by
  /// [ 0, (*pVecSeqLengths)[seq_id] ), then it is trimmed if parts of the
  /// interval overlaps the bounds, or dropped if it is outside.  If you 
  /// do not specify pVecSeqLengths, see below for its effect on the various
  /// Get*() methods.
  void CreateCoverages(const vec<seq_interval>& vecIntervals,
		       const vec<int>* pVecSeqLengths = NULL);

  // The following methods will not return windows with coverage 0
  // unless the pVecSeqLengths parameter was nonzero in the
  // constructor or the last call to CreateCoverages().  If the
  // pVecSeqLengths parameter is zero, the code can't figure out where
  // the sequence(s) whose coverage is being checked ends, only the
  // last base that is covered, and so it cannot deduce where the
  // final window (possibly of zero coverage) ends.

  // Therefore, if pVecSeqLengths is zero, GetAllCoverages() and
  // GetCoveragesAtMost() will not return intervals of coverage 0;
  // GetAllCoveragesBetween() and GetCoveragesAtLeast() will not
  // produce the expected output if minValue is 0; and
  // GetCoveragesExactly() will produce nothing if value is 0.

  /// returns all the coverage regions, such that no two adjoining regions
  /// will have the same coverage amount. The level of coverage is saved
  /// in the interval_id member of each seq_interval.
  void GetAllCoverages(vec<seq_interval>& vecIntervals) const;

  /// returns one seq_interval per window_size bases on each sequence,
  /// setting the interval_id to the average coverage level in that
  /// window.  since the coverage level in a window is generally not
  /// an integer, the value is scaled by the given amount.
  ///
  /// Only usable if the sequence lengths were passed in with the
  /// constructor or with CreateCoverages().
  void GetAllCoverageInWindows(vec<seq_interval>& vecIntervals,
                               const unsigned int windowSize,
                               const unsigned int scale ) const;

  /// Return the number of bases with coverage 0...n in counts
  void CountAllCoverages(vec<longlong> & counts) const;

  /// returns all coverage regions where the coverage amount is:
  /// minValue <= coverage <= maxValue
  /// if minValue = 0, this includes regions with no ReadLocations
  void GetCoveragesBetween(unsigned int minValue, unsigned int maxValue,
			   vec<seq_interval>& vecIntervals) const;
  
  /// returns all coverage regions where the coverage amount is:
  /// minValue <= coverage
  inline void GetCoveragesAtLeast(unsigned int minValue,
				  vec<seq_interval>& vecIntervals) const
  { GetCoveragesBetween(minValue, 0xFFFFFFFF, vecIntervals); }

  /// returns all coverage regions where the coverage amount is:
  /// coverage = value
  inline void GetCoveragesExactly(unsigned int value,
				  vec<seq_interval>& vecIntervals) const
  { GetCoveragesBetween(value, value, vecIntervals); }

  /// returns all coverage regions where the coverage amount is:
  /// coverage <= maxValue
  inline void GetCoveragesAtMost(unsigned int maxValue,
				 vec<seq_interval>& vecIntervals) const
  { GetCoveragesBetween(0, maxValue, vecIntervals); }


  /// Print statistics of coverage for each 
  /// contig: contig#, mean, Q3/Q1
  void PrintCoverageStats(ostream & os, bool printZeroCoverageContigs=true);


  /// Print coverage for each base in three columns: sequenceId, base,coverage.
  /// Only prints lines for bases with coverage 0 if printZeros=true.
  /// If printLastInContig, then always print last base in contig even if 
  /// it has 0 coverage (so we know where contig ends).
  void PrintTextCoverage(ostream & os, bool printLastInContig=false,
			 bool printZeros = false) const;
    
private:
  class CoverageChange {
  public:
    CoverageChange( int id, int pos, int amt )
      : seqId( id ), coord( pos ), amount( amt ) 
    {}

    int seqId;
    int coord;
    int amount;

    bool IsEmpty() const { return amount == 0; }

    bool operator< ( const CoverageChange& other ) const
    {
      return ( seqId < other.seqId || 
	       seqId == other.seqId && coord < other.coord );
    }

    friend ostream& operator<< ( ostream& out, const CoverageChange& cc )
    {
      return out << cc.seqId 
		 << '@' << cc.coord
		 << ( cc.amount > 0 ? "+" : "" ) << cc.amount;
    }
  };
    

  vec<CoverageChange> m_covChanges;

  vec<int>* m_pSeqLengths;

  ostream* m_pOStream;
};


#endif
