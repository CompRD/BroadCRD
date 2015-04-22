/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef KMERFREQ_KMERFREQUENCYTABLE_H
#define KMERFREQ_KMERFREQUENCYTABLE_H

#include "Basevector.h"
#include "String.h"
#include "Vec.h"

#include "kmer_freq/KmerShortMap.h"

/**
   Class: KmerFrequencyTable
  
   A histogram of kmer frequencies: for each frequency range, how many kmers are
   there in the read set with that frequency?

   Not to be confused with <KmerShortMap>, which _for each kmer_ can represent the
   frequency of that kmer.  A KmerFrequencyTable, by contrast, for each frequency
   represents how many kmers there are with that particular frequency -- without
   representing the actual list of kmers with each possible frequency.
   A KmerFrequencyTable can be constructed from a KmerShortMap associating each
   kmer with its frequency.

   A KmerFrequencyTable _contains_ a KmerShortMap giving the frequency of each kmer,
   and has methods to access this information.  For example, <GetCount()> returns
   the number of occurrences of a given kmer, querying the underlying KmerShortMap
   for this information.  In fact, most of the methods are delegated to the
   underlying KmerShortMap, which does the actual work.
*/
class KmerFrequencyTable {
 public:
  KmerFrequencyTable( const KmerShapeId& KSHAPE, const String& filename );

 private:
  /// Unimplemented: KmerFrequencyTables cannot be copied.
  KmerFrequencyTable( const KmerFrequencyTable& );
  /// Unimplemented: KmerFrequencyTables cannot be copied.
  KmerFrequencyTable& operator= ( const KmerFrequencyTable& );

 public:
  ~KmerFrequencyTable();

  int GetKmerSize() const {
    return m_kmerMap.GetKmerSize();
  }

  kmer_freq_t GetMaxCount() const {
    return m_kmerMap.GetMaxValue();
  }

  // Method: GetCount
  // 
  // Returns a count equal to min(count,GetMaxCount()) if the
  // frequency of the kmer is > 0 and -1 otherwise.  The basevector
  // "kmer" should be Canonicalized.
  kmer_freq_t GetCount( const basevector& kmer ) const {
    return m_kmerMap.GetValue( kmer );
  }

  // Method: GetCounts
  // For each kmer in seq, sets corresponding entry in counts to
  // GetCount() for that kmer.
  void GetCounts( const basevector& seq, vec<int>& counts ) const {
    m_kmerMap.GetValues( seq, counts );
  }

  /// Method: IsStrong
  /// Returns true if every kmer in seq has count >= minFreq.
  bool IsStrong( const basevector& seq, int minFreq ) const {
    return m_kmerMap.IsStrong( seq, minFreq );
  }

  /// Method: IsStrong by gc content
  /// Returns true if every kmer in seq has count >= minFreq[gc],
  /// where gc is the number of GC bases in the kmer.
  bool IsStrong( const basevector& seq, const vec<int>& minFreqByGC ) const {
    return m_kmerMap.IsStrong( seq, minFreqByGC );
  }

  /// Method: GetHistogram
  /// Calculates histogram from counts for all kmers.
  void GetHistogram( vec<kmer_count_t>& hist ) const {
    m_kmerMap.GetHistogram( hist );
  }

  /// Method: GetGCHistograms
  /// Calculates histogram from counts for kmers by GC content.
  void GetGCHistograms( vec< vec<kmer_count_t> >& histByGC ) const {
    m_kmerMap.GetGCHistograms( histByGC );
  }

  /// Calculates two histograms, given a 2nd KmerShortMap that contains 'true' kmers,
  /// one from those kmers that are 'true', and one from those kmers that are 'false'.
  void GetGCHistograms( vec< vec<kmer_count_t> >& trueHist,
			vec< vec<kmer_count_t> >& falseHist,
			KmerShortMap& truth ) const {
  return m_kmerMap.GetGCHistograms( trueHist, falseHist, truth );
}

  /**
     Method: GetPeak
     
     Returns the first kmer frequency (of at least 2) frequency with the
     max kmer count.  This should correspond to the expected frequency
     in the reads of a kmer which occurs exactly once in the genome.
  */
  int GetPeak( const int gc_content = -1 ) const;

  /// Method: GetFirstLocalMin
  /// Returns what appears to be the first local minimum, if any, with
  /// the option to restrict to kmers with particular GC content.  If
  /// there are fewer than min_samples or there appears to be no local
  /// minimum frequency, -1 is returned.
  int GetFirstLocalMin( const double threshold = 2.0, 
                        const int gc_content = -1,
                        const int min_samples = 0 ) const;

 private:
  void GetGCBands( const vec<kmer_count_t>& totalByGC,
                   const int min_samples,
                   const int maxSampleRatio,
                   const bool includeTails,
                   vec< pair<int,int> >& bands ) const;

  int GetFLM( const vec<kmer_count_t>& hist, 
              const double threshold, 
              const int minSamples ) const;

  KmerShortMap m_kmerMap;

  struct cache {
    cache() : peak(-1),localMin(-2) {}
    int peak;
    vec<int> peakGC;
    int localMin;
    vec<int> localMinGC;
  };

  cache* m_pCache;
};

#endif

// Synonyms: Various synonyms
//
//    kmer frequency table - See <KmerFrequencyTable>
//    KmerFrequencyTable.h - See <KmerFrequencyTable>

