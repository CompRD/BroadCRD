/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "feudal/BinaryStream.h"
#include "kmer_freq/KmerShortMap.h"
#include "math/Functions.h"
#include "kmers/KmerRecord.h"
#include "Vec.h"

/*
   Class: KmerShortMap::Imp

   Abstract base class for concrete implementations of <KmerShortMap>.
*/
class KmerShortMap::Imp {
 public:
  Imp() {}
  
  virtual ~Imp() {}
  
  virtual int GetKmerSize() const = 0;
  
  virtual value_t GetValue( const basevector& kmer ) = 0;
  
  virtual void GetValues( const basevector& seq, vec<value_t>& values ) = 0;
  
  virtual bool IsStrong( const basevector& seq, value_t minValue ) = 0;

  virtual bool IsStrong( const basevector& seq, int endPos, int startPos, value_t minValue ) = 0;
  
  virtual bool IsStrong( const basevector& seq, const vec<value_t>& minValueByGC ) = 0;

  virtual void GetHistogram( vec<kmer_count_t>& hist ) const = 0;

  virtual void GetGCHistograms( vec< vec<kmer_count_t> >& hist ) const = 0;

  virtual void GetGCHistograms( vec< vec<kmer_count_t> >& trueHist,
				vec< vec<kmer_count_t> >& falseHist,
				KmerShortMap& truth ) const = 0;
};

/*
   Class: KmerShortMap::ConcreteImp

   A concrete implementation of <KmerShortMap::Imp>.
*/
template <class KSHAPE>
class KmerShortMap::ConcreteImp : public KmerShortMap::Imp {
 private:
  static const int K = KSHAPE::KSIZE;

  static const value_t s_maxValue = kmer_with_count<K>::max_count;

  void MakeIndex()
  {
    const vec< kmer_with_count<K> >& kmers = mKmers;
    m_K0 = ( K > 14 ? 14 : 8 );

    kmer_count_t fourToK0 = 1;
    for ( int i = 0; i < m_K0; i++ )
      fourToK0 *= 4;

    m_index.resize( fourToK0 + 1, -1 );

    if ( K >= 16 ) 
      m_indexShift = 2 * (16 - m_K0);
    else
      m_indexShift = 2 * (8 - m_K0);
    m_indexMask = (1 << m_K0*2) - 1;
    
    for ( kmer_count_t i = kmers.size( ) - 1; i >= 0; i-- ) {
      unsigned int x;
      if ( K >= 16 )
        x = *kmers[i].Ints( );
      else
        x = *kmers[i].Shorts( );
      x >>= m_indexShift;
      x &= m_indexMask;
      m_index[ x ] = i;    
    }
    m_index[fourToK0] = kmers.size( );

    for ( kmer_count_t i = fourToK0 - 1; i >= 0; i-- )
      if ( m_index[i] == -1 ) 
        m_index[i] = m_index[i+1];
  }

  bool IsIndexed() const {
    return ( ! m_index.empty() );
  }

 public:
  ConcreteImp( const String& filename )
    : KmerShortMap::Imp()
  {
    BinaryReader::readFile( filename, &mKmers );
  }


  int GetKmerSize() const {
    return K;
  }

  struct kmerOnlyLess
    : public binary_function<kmer_with_count<K>,kmer_with_count<K>,bool>
  {
    bool operator() ( const kmer_with_count<K>& lhs, const kmer_with_count<K>& rhs ) const {
      return lt_kmer( lhs, rhs );
    }
  };

  value_t GetValue( const basevector& kmer ) {
    if ( ! IsIndexed() )
      MakeIndex();

    AssertEq( kmer.size(), static_cast<unsigned>(K) );

    kmer_with_count<K> target( kmer, 0 );
    
    unsigned int indexEntry = ( K >= 16 ? *target.Ints( ) : (unsigned int) *target.Shorts( ) );
    indexEntry >>= m_indexShift;
    indexEntry &= m_indexMask;

    switch ( m_index[indexEntry+1] - m_index[indexEntry] ) {
      case 0:
        return -1;
      case 1:
        if ( eq_kmer( target, mKmers[m_index[indexEntry]] ) )
          return mKmers[m_index[indexEntry]].Count();
        else
          return -1;
      default:
        typename vec< kmer_with_count<K> >::const_iterator iter = 
          lower_bound( mKmers.begin( ) + m_index[indexEntry],
                       mKmers.begin( ) + m_index[indexEntry+1],
                       target,
                       kmerOnlyLess() ); 
    
        if ( iter != mKmers.begin() + m_index[indexEntry+1] &&
             eq_kmer( target, *iter ) )
          return iter->Count();
        else
          return -1;
    }
  }


  void GetValues( const basevector& seq, vec<value_t>& values ) {
    int numKmers = seq.size() - KSHAPE::KSPAN + 1;
    if ( numKmers < 0 )
      numKmers = 0;
    values.resize( numKmers );

    basevector kmer(K);
    for ( int kmerStart = 0; kmerStart < numKmers; ++kmerStart ) {
      KSHAPE::extractKmer( kmer, seq, kmerStart );
      kmer.Canonicalize();
      values[ kmerStart ] = this->GetValue( kmer );
    }
  }


  bool IsStrong( const basevector& seq, value_t minValue ) {
    const int numKmers = seq.size() - KSHAPE::KSPAN + 1;

    if ( numKmers < 0 ) return false;

    basevector kmer(K);
    for ( int kmerStart = numKmers - 1; kmerStart >= 0; --kmerStart ) {
      KSHAPE::extractKmer( kmer, seq, kmerStart );
      kmer.Canonicalize();

      if ( this->GetValue( kmer ) < minValue )
        return false;
    }
    
    return true;
  }


  bool IsStrong( const basevector& seq, int startPos, int endPos, value_t minValue ) {
    const int numKmers = endPos - startPos + 1 - KSHAPE::KSPAN + 1;

    if ( numKmers < 0 ) return false;

    basevector kmer(K);
    for ( int kmerStart = startPos + numKmers - 1; kmerStart >= startPos; --kmerStart ) {
      KSHAPE::extractKmer( kmer, seq, kmerStart );
      kmer.Canonicalize();

      if ( this->GetValue( kmer ) < minValue )
	return false;
 
    }
    
    return true;
  }


  bool IsStrong( const basevector& seq, const vec<value_t>& minValueByGC ) {
    AssertEq( minValueByGC.size(), K+1u );

    const int numKmers = seq.size() - KSHAPE::KSPAN + 1;

    if ( numKmers < 0 ) return false;

    kmer_gc_content_t gc = seq.GcBases( 0, K );

    basevector kmer(K);
    for ( int kmerStart = 0; kmerStart < numKmers; ++kmerStart ) {
      KSHAPE::extractKmer( kmer, seq, kmerStart );
      kmer.Canonicalize();

      if ( this->GetValue( kmer ) < minValueByGC[gc] )
        return false;

      gc -= IsGC( seq[kmerStart] );
      if ( kmerStart+K < (int)seq.size() )
        gc += IsGC( seq[kmerStart+K] );
    }
    
    return true;
  }


  void GetHistogram( vec<kmer_count_t>& hist ) const {
    hist.clear();
    hist.resize(s_maxValue+1,0);

    typename vec< kmer_with_count<K> >::const_iterator iKmer;
    for ( iKmer = mKmers.begin(); iKmer != mKmers.end(); ++iKmer )
      ++hist[ iKmer->Count() ];
  }

  /// Get histograms by kmer GC content
  void GetGCHistograms( vec< vec<kmer_count_t> >& histByGC ) const {
    histByGC.clear();
    histByGC.resize( K+1, vec<kmer_count_t>(s_maxValue+1,0) );
    basevector kmer;
    typename vec< kmer_with_count<K> >::const_iterator iKmer;
    for ( iKmer = mKmers.begin(); iKmer != mKmers.end(); ++iKmer ) {
      iKmer->GetBasevector(kmer);
      kmer_gc_content_t gc = kmer.GcBases();
      histByGC[gc][iKmer->Count()]++;
    }
  }

  /// Get histograms by kmer GC content, split into those kmers
  /// that are true and those which are false - according to the
  /// kmers contained in the KmerShortMap truth.
  void GetGCHistograms( vec< vec<kmer_count_t> >& trueHistByGC,
			vec< vec<kmer_count_t> >& falseHistByGC,
			KmerShortMap& truth ) const {
    trueHistByGC.clear();
    trueHistByGC.resize( K+1, vec<kmer_count_t>(s_maxValue+1,0) );
    falseHistByGC.clear();
    falseHistByGC.resize( K+1, vec<kmer_count_t>(s_maxValue+1,0) );
    basevector kmer;
    typename vec< kmer_with_count<K> >::const_iterator iKmer;
    for ( iKmer = mKmers.begin(); iKmer != mKmers.end(); ++iKmer ) {
      iKmer->GetBasevector(kmer);
      kmer_gc_content_t gc = kmer.GcBases();
      if (truth.IsStrong(kmer))
	trueHistByGC[gc][iKmer->Count()]++;
      else
	falseHistByGC[gc][iKmer->Count()]++;
    }
  }


 private:
  vec< kmer_with_count<K> > mKmers;
  int m_K0;
  vec<kmer_count_t> m_index;
  int m_indexShift;
  int m_indexMask;
};


KmerShortMap::KmerShortMap( const KmerShapeId& kmerShapeId, const String& filename ) {
#define MK_CONCRETE_IMP(_KSHAPE) m_pImp = new ConcreteImp<_KSHAPE>( filename )
  DISPATCH_ON_KSHAPE(kmerShapeId, MK_CONCRETE_IMP);
}

KmerShortMap::~KmerShortMap() {
  delete m_pImp;
}

int KmerShortMap::GetKmerSize() const {
  return m_pImp->GetKmerSize();
}

KmerShortMap::value_t KmerShortMap::GetValue( const basevector& kmer ) const {
  return m_pImp->GetValue( kmer );
}

void KmerShortMap::GetValues( const basevector& seq, vec<value_t>& values ) const {
  m_pImp->GetValues( seq, values );
}

bool KmerShortMap::IsStrong( const basevector& seq, value_t minValue ) const {
  return m_pImp->IsStrong( seq, minValue );
}

bool KmerShortMap::IsStrong( const basevector& seq, int startPos, int endPos, value_t minValue ) const {
  return m_pImp->IsStrong( seq, startPos, endPos, minValue );
}

bool KmerShortMap::IsStrong( const basevector& seq, const vec<value_t>& minValueByGC ) const {
  return m_pImp->IsStrong( seq, minValueByGC );
}

void KmerShortMap::GetHistogram( vec<kmer_count_t>& hist ) const {
  return m_pImp->GetHistogram( hist );
}

void KmerShortMap::GetGCHistograms( vec< vec<kmer_count_t> >& hist  ) const {
  return m_pImp->GetGCHistograms( hist );
}

void KmerShortMap::GetGCHistograms( vec< vec<kmer_count_t> >& trueHist,
				    vec< vec<kmer_count_t> >& falseHist,
				    KmerShortMap& truth ) const {
  return m_pImp->GetGCHistograms( trueHist, falseHist, truth );
}


#define INSTANTIATE(KSHAPE, dummy) template class KmerShortMap::ConcreteImp<KSHAPE>
FOR_ALL_KSHAPES(INSTANTIATE,unused);

