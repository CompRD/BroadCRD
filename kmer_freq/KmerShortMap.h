/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef KMERFREQ_KMERSHORTMAP_H
#define KMERFREQ_KMERSHORTMAP_H

#include "String.h"
#include "Vec.h"
#include "Basevector.h"
#include "kmers/KmerShape.h"

// Type: kmer_count_t
//
// Represents the actual count of different kmers sharing some property -- for example,
// the number of kmers with a given frequency of occurrence in the reads.
// Must be signed, to allow for the "NULL" value of -1.
typedef longlong kmer_count_t;

// Type: kmer_freq_t
//
// Represents the number of occurrences of a particular kmer in the genome.
// Must be signed, to allow for the "NULL" value of -1 that differs from any
// possible kmer frequency.
typedef int kmer_freq_t;

// Type: kmer_gc_content_t
//
// Logical type for representing the gc content of a kmer;
// for a kmer of size K, this can vary from 0 to K, so there are K+1
// possible values.
// Must be signed to allow for the "NULL" value of -1 that differs from
// ever possible gc content.
typedef int kmer_gc_content_t;

/*
   Class: KmerShortMap
  
   Class that maps kmers to short int values -- not in the sense of assigning
   these int values as ids, but as just a way to associate a short int value
   with each kmer.  This value can be the kmer's frequency
   (see <KmerFrequencyTable>), a Boolean value indicating whether we think
   the kmer is genomic, or some other value.  The value can even be ignored --
   then we treat the KmerShortMap as just a set of kmers, for example
   <genomic kmers>.

   A KmerShortMap can be constructed from a file on disk.  Currently,
   that file is created by <FindKmerFrequencies>.

   Transformations on KmerShortMaps can be done using <TransformKmerShortMap()>.

   KmerShortMaps are stored in disk in <KmerShortMap files>.  These can be used
   as follows:

   (begin example)
   vec< kmer_with_count<K> > kmers;
   BinaryReader::readFile( filename_in, &kmers );
   BinaryIteratingWriter< vec<kmer_with_count<K> > > writer( filename_out );
   ...   writer.write( kmer_with_count<K>( theKmer, new_value ) );
   writer.close();
   (end example)

   A KmerShortMap can be seen as a wrapper for a vec< kmer_with_count<K> >, which provides
   useful methods for working with this set of (kmer, short int) pairs -- such as building
   kmer frequency histograms, and builds an index of the kmers in the map so that the
   short int value associated with a given kmer can be quickly retrieved.
   The original vec< kmer_with_count<K> > is just a list of kmers, and searching for a
   given kmer in the list can take long.
*/
class KmerShortMap {
 public:
  /// Constructor that loads kmer short map from file
  KmerShortMap( const KmerShapeId& kmerShapeId, const String& filename );

  /// Type: value_t
  /// The type for values that this map associates with each kmer.
  typedef int value_t;

 private:
  /// Unimplemented: KmerShortMaps cannot be copied.
  KmerShortMap( const KmerShortMap& );
  /// Unimplemented: KmerShortMaps cannot be copied.
  KmerShortMap& operator= ( const KmerShortMap& );

 public:
  ~KmerShortMap();

  /// Returns the size of the kmers that are mapped.
  int GetKmerSize() const;

  /// Returns the upper bound on value that can be returned.
  value_t GetMaxValue() const {
    return USHRT_MAX;
  }

  /// Returns a value equal to min(value,GetMaxValue()) if value > 0 and -1
  /// otherwise.  The basevector "kmer" should be Canonicalized.
  value_t GetValue( const basevector& kmer ) const;

  /// For each kmer in seq, sets corresponding entry in values to
  /// GetValue() for that kmer.
  void GetValues( const basevector& seq, vec<value_t>& values ) const;

  /// Returns true if every kmer in seq has value >= minValue.
  /// If minValue = 0 then returns true if all kmers in the seq
  /// are present in the map, regardless of their values.
  bool IsStrong( const basevector& seq, value_t minValue = 0 ) const;

  /// Returns true if every kmer in the seq interval has value >= minValue.
  /// If minValue = 0 then returns true if all kmers in the seq interval
  /// are present in the map, regardless of their values.
  bool IsStrong( const basevector& seq, int startPos, int endPos, value_t minValue = 0 ) const;
  
  /// Returns true if every kmer in seq has value >= minValue[gc], where
  /// gc is the number of GC bases in the kmer.
  bool IsStrong( const basevector& seq, const vec<value_t>& minValueByGC ) const;

  /// Calculates histogram from values for all kmers
  void GetHistogram( vec<kmer_count_t>& hist ) const;

  /// Calculates histogram from values for kmers by GC content
  void GetGCHistograms( vec< vec<kmer_count_t> >& histByGC ) const;

  /// Given a 2nd KmerShortMap that contains 'true' kmers, generates two histograms:
  /// one from those kmers that are 'true', and one from those kmers that are
  /// 'false'. The sum of these two histograms is the histogram produced by
  /// GetGCHistograms( vec< vec<kmer_count_t> >& histByGC) above.
  void GetGCHistograms( vec< vec<kmer_count_t> >& trueHistByGc,
			vec< vec<kmer_count_t> >& falseHistByGc,
			KmerShortMap& truth )  const;


 private:
  class Imp;

  Imp* m_pImp;

  template <class KSHAPE> class ConcreteImp;
};  // class KmerShortMap

/**
   Class: KmerSet

   Represents a set of kmers.  Lets you efficiently test a kmer for membership in the set.

   Implemented as a <KmerShortMap>, but makes it more clear that we're working with a set and not
   a map, via class name and method names.
*/
class KmerSet: private KmerShortMap {
 public:

  /// Constructor: KmerSet constructor
  /// Loads this kmer set from a <KmerShortMap file>. 
  KmerSet( const KmerShapeId& kmerShapeId, const String& filename ):
    KmerShortMap( kmerShapeId, filename ) { }

  /// Method: containsKmer
  /// Tells whether this set contains the given kmer.
  inline bool ContainsKmer( const kmer_t& kmer ) const { return IsStrong( kmer ); }
};  // class KmerSet

#endif
// #ifndef KMERFREQ_KMERSHORTMAP_H
