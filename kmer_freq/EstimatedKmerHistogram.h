// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ESTIMATEDKMERHISTOGRAM_H
#define ESTIMATEDKMERHISTOGRAM_H

#include "Basevector.h"
#include "Histogram.h"
#include "Vec.h"

// A KmerHistogram produces an approximate count of the kmers in a set
// of basevectors.  It is approximate because the counting is done
// quickly and in a relatively small amount of memory by using a hash
// table with a count for each hash value.  Hence, two kmers that have
// a hash collision will have their counts conflated.  

// Some care has been taken in choosing a hash that minimizes the
// clustering of these collisions, but regardless, this table should
// be sparse to reduce collisions when mapping the hash value to an
// index into the hash table.


class KmerHistogram;

// KmerHistograms cannot be instantiated directly; they should be
// created via this helper function.  See KmerHistogram.cc for
// allowable values of K.
KmerHistogram* GetNewKmerHistogram( const int K, const size_t capacity );


class KmerHistogram {

 public:
  KmerHistogram( const int K )
    : m_K( K )
  {}

  virtual ~KmerHistogram() {}

  int GetKmerSize() { return m_K; }
  
  /// Increment the count for each kmer found in bases.
  void CountKmers( const basevector& bases )
  {
    this->DoCountKmers( bases );
  }

  /// Resize counts to the number of kmers in bases, and fill it with
  /// the count for the corresponding kmer.
  void GetCounts( const basevector& bases, vec<int>& counts ) const
  {
    this->DoGetCounts( bases, counts );
  }

  /// Resize counts to the number of kmers in bases, and fill it with
  /// the count for every fourth kmer (0,4,8...) and zero for all
  /// other kmers.  More than four times faster than GetCounts().
  void GetSomeCounts( const basevector& bases, vec<int>& counts ) const
  {
    this->DoGetSomeCounts( bases, counts );
  }

  /// Get the approximate number of distinct kmers counted so far
  /// (modulo hash collisions).
  size_t GetNumDistinctKmers() const
  {
    return this->DoGetNumDistinctKmers();
  }
  
  /// Get the median kmer frequency of all kmers that appear more than once.
  int GetMedianNonuniqueKmerCount() const
  {
    return this->DoGetMedianNonuniqueKmerCount();
  }

  /// Fill out the given histogram with the observed kmer frequencies.
  void FillKmerCountHistogram( histogram<int>& theHistogram ) const
  {
    return this->DoFillKmerCountHistogram( theHistogram);
  }

 protected:
  int m_K;

  virtual void DoCountKmers( const basevector& bases ) = 0;

  virtual void DoGetCounts( const basevector& bases, vec<int>& counts ) const = 0;

  virtual void DoGetSomeCounts( const basevector& bases, vec<int>& counts ) const = 0;

  virtual size_t DoGetNumDistinctKmers() const = 0;

  virtual int DoGetMedianNonuniqueKmerCount() const = 0;

  virtual void DoFillKmerCountHistogram( histogram<int>& theHistogram ) const = 0;
};


#endif
