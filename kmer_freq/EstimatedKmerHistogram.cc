// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "kmer_freq/EstimatedKmerHistogram.h"

#if 0
template <int K>
unsigned int HashB( const basevector& bv, const int offset )
{
  STATIC_ASSERT_M( K%4 == 0, bad_kmer_size );

  AssertGe( bv.size() - 4*offset, K );

  // This is a 64-bit FNV-1a hash XOR-folded to 32 bits.
  ulonglong hash = ULLCONST(14695981039346656037);
  const ulonglong prime = ULLCONST(1099511628211);
  const unsigned int numBytes = K/4;
  for ( unsigned int ii = 0; ii < numBytes; ++ii )
  {
    hash ^= bv.DataAsBytes(offset+ii);
    hash *= prime;
  }
  const ulonglong mask32 = ULLCONST(1<<32)-1;
  hash = (hash>>32) ^ (hash & mask32);

  return static_cast<unsigned int>(hash);
}
#endif

template <int K>
class KmerHashHistogram : public KmerHistogram
{
 public:
  KmerHashHistogram( size_t capacity )
    : KmerHistogram( K ),
      m_numFreeEntries( 0 ),
      m_hashMask( 0 )
  {
    this->ResizeHash( capacity );
  }

  virtual ~KmerHashHistogram()
  {
  }

 protected:
  struct Entry {
    unsigned char count;

    Entry() : count(0) {}
  };


  vec<Entry> m_table;
  size_t     m_numFreeEntries;
  size_t     m_hashMask;

  virtual void DoCountKmers( const basevector& bases );
  
  virtual void DoGetCounts( const basevector& bases, vec<int>& counts ) const;

  virtual void DoGetSomeCounts( const basevector& bases, vec<int>& counts ) const;

  virtual int DoGetMedianNonuniqueKmerCount() const;

  virtual void DoFillKmerCountHistogram( histogram<int>& theHistogram ) const;

  virtual size_t DoGetNumDistinctKmers() const;

  const Entry* LookupKmerInHash( const basevector& bv, const int offset ) const;
  Entry*       LookupKmerInHash( const basevector& bv, const int offset );

  size_t GetIndexOfKmerInHash( const basevector& bv, const int offset ) const;

  void ResizeHash( const size_t newSize );

 private:
  KmerHashHistogram<K>* mp_histogram;

  class CountKmersAction;
  friend class CountKmersAction;

  class GetCountsAction;
  friend class GetCountsAction;
};


template <int K>
void KmerHashHistogram<K>::ResizeHash( const size_t newMinSize )
{
  size_t newCapacity = ( m_table.size() < 1024 ? 1024 : m_table.size() );

  while ( newCapacity < newMinSize )
    newCapacity *= 2;

  m_hashMask = newCapacity - 1;

  vec<Entry> newTable( newCapacity );

  m_numFreeEntries += ( newTable.size() - m_table.size() );
  m_table.swap( newTable );
}
 

template <int K>
size_t KmerHashHistogram<K>::GetIndexOfKmerInHash( const basevector& bv, const int offset ) const
{
  AssertEq(K%4,0);
  return bv.hash(offset,K/4) & m_hashMask;
}


template <int K>
const typename KmerHashHistogram<K>::Entry* 
KmerHashHistogram<K>::LookupKmerInHash( const basevector& bv, const int offset ) const
{
  size_t index = this->GetIndexOfKmerInHash( bv, offset );

  return &(m_table[index]);
}


template <int K>
typename KmerHashHistogram<K>::Entry* 
KmerHashHistogram<K>::LookupKmerInHash( const basevector& bv, const int offset ) 
{
  size_t index = this->GetIndexOfKmerInHash( bv, offset );

  return &(m_table[index]);
}


/// This is an algorithm, more or less, that walks with some
/// efficiency across the kmers in bases (but not in order) and
/// applies the specified Action::operator()(int, const basevector&,
/// int) on each of them.  This is so we don't have to duplicate this
/// walking functionality in DoCountKmers(), DoGetCounts(), and
/// DoGetSomeCounts().
///
/// The noShift option determines whether the kmers starting at an
/// index not divisible by four are Acted upon.  If it is true, then no
/// shifting need take place, which saves not only the time spent in
/// processing those kmers, but also the time spent extracting those
/// kmers (whose unclean division by four requires a great deal of
/// shifting about in the basevector).  This means that noShift=true
/// takes less than a fourth of the time of noShift=false.
/// 
/// In general, call WalkAllKmers<K>() or WalkSomeKmers<K>(), which
/// will call WalkKmers<K> with the right noShift parameter value.
template <int K, class Action>
void WalkKmers( const basevector& bases, Action theAction, bool noShift = false )
{
  int numKmers = bases.size() - K + 1;

  if ( numKmers <= 0 )
    return;

  const int basesPerByte = 4;
  basevector shiftedBasesFw;
  basevector shiftedBasesRc;
  for ( int shift = 0; shift < basesPerByte; ++shift )
  {
    // Copy bases into shiftedBasesFw so that every ith kmer, where
    // i%4 == shift, is byte-aligned.
    shiftedBasesFw.SetToSubOf( bases, shift, bases.size() - shift );
    
    // Copy these bases into the rc basevector.
    shiftedBasesRc = shiftedBasesFw;
    // Add enough space to the rc version such that when it is rc'ed,
    // the rc's of the byte-aligned kmers above are also byte aligned.
    const int extraSpace = ( basesPerByte - shiftedBasesRc.size() % basesPerByte );
    shiftedBasesRc.resize( shiftedBasesRc.size() + extraSpace );
    shiftedBasesRc.ReverseComplement();
    const int numBytesInRc = shiftedBasesRc.size()/basesPerByte;

    const int numBytesPerKmer = K/basesPerByte;
    // The number of bytes that start kmers:
    const int numStartBytes = (numKmers+basesPerByte-1)/basesPerByte;
    for ( int fwStartByteIdx = 0; fwStartByteIdx < numStartBytes; ++fwStartByteIdx )
    {
      // Which byte is the first in the rc version of this kmer?
      int rcStartByteIdx = numBytesInRc - fwStartByteIdx - numBytesPerKmer;

      basevector* p_bvWithLesserKmer = &shiftedBasesFw;
      int startByteIdx = fwStartByteIdx;

      int fwByteIdx = fwStartByteIdx;
      int rcByteIdx = rcStartByteIdx;

      // Switch which basevector we're pointing to if the rc turns out to be less than the fw.
      for ( int ii = 0; ii < numBytesPerKmer; ++ii, ++fwByteIdx, ++rcByteIdx )
      {
          unsigned char fwByte = shiftedBasesFw.extractKmer(fwByteIdx*4,4);
          unsigned char rcByte = shiftedBasesRc.extractKmer(rcByteIdx*4,4);
          if ( fwByte < rcByte )
              break;
          if ( rcByte < fwByte )
          {
              p_bvWithLesserKmer = &shiftedBasesRc;
              startByteIdx = rcStartByteIdx;
              break;
          }
      }

      (theAction)( fwStartByteIdx * basesPerByte + shift, *p_bvWithLesserKmer, startByteIdx );
    }

    if ( noShift )
      break;

    --numKmers;
    if ( numKmers <= 0 )
      break;
  }
}

template <int K, class Action>
void WalkAllKmers( const basevector& bases, Action theAction )
{
  WalkKmers<K>( bases, theAction, false );
}

template <int K, class Action>
void WalkSomeKmers( const basevector& bases, Action theAction )
{
  WalkKmers<K>( bases, theAction, true );
}

// This action, specifically for use with WalkKmers above, looks up
// the hash table entry corresponding to the given kmer and increments
// its count.
template <int K>
class KmerHashHistogram<K>::CountKmersAction
{
 public:
  CountKmersAction( KmerHashHistogram<K>* p_histogram )
    : mp_histogram( p_histogram )
  {}
  
  void operator() ( const int kmerIdx, const basevector& bv, const int byteOffset ) const
  {
    Entry* p_entry = mp_histogram->LookupKmerInHash( bv, byteOffset );
    
    if ( p_entry->count == 0 )
      --mp_histogram->m_numFreeEntries;
    
    if ( p_entry->count != 255 )
      ++(p_entry->count);
  }

 private:
  KmerHashHistogram<K>* mp_histogram;
};

template <int K>
void KmerHashHistogram<K>::DoCountKmers( const basevector& bases )
{
  WalkAllKmers<K>( bases, CountKmersAction( this ) );
}

             
// This action, specifically for use with WalkKmers above, looks up
// the hash table entry corresponding to the given kmer and copies its
// count to the given index of the counts vector supplied in the
// constructor.
template <int K>
class KmerHashHistogram<K>::GetCountsAction
{
 public:
  GetCountsAction( const KmerHashHistogram<K>* p_histogram, vec<int>* p_counts )
    : mp_histogram( p_histogram ),
      mp_counts( p_counts )
  {}
  
  void operator() ( const int kmerIdx, const basevector& bv, const int byteOffset ) const
  {
    const Entry* p_entry = mp_histogram->LookupKmerInHash( bv, byteOffset );

    (*mp_counts)[ kmerIdx ] = p_entry->count;
  }

 private:
  const KmerHashHistogram<K>* mp_histogram;
  vec<int>* mp_counts;
};


// Documented in header.
template <int K>
void KmerHashHistogram<K>::DoGetCounts( const basevector& bases, vec<int>& counts ) const
{
  counts.clear();
  
  int numKmers = bases.size() - K + 1;
  
  if ( numKmers <= 0 )
    return;

  counts.resize( numKmers, 0 );

  WalkAllKmers<K>( bases, GetCountsAction( this, &counts ) );
}


// Documented in header.
template <int K>
void KmerHashHistogram<K>::DoGetSomeCounts( const basevector& bases, vec<int>& counts ) const
{
  counts.clear();
  
  int numKmers = bases.size() - K + 1;
  
  if ( numKmers <= 0 )
    return;

  counts.resize( numKmers, 0 );

  WalkSomeKmers<K>( bases, GetCountsAction( this, &counts ) );
}


// Documented in header.
template <int K>
size_t KmerHashHistogram<K>::DoGetNumDistinctKmers() const
{
  return m_table.size() - m_numFreeEntries;
}  


// Documented in header.
template <int K>
int KmerHashHistogram<K>::DoGetMedianNonuniqueKmerCount() const
{
  vec<unsigned char> allCounts;
  
  size_t numDistinctKmers = this->GetNumDistinctKmers();

  if ( numDistinctKmers == 0 )
    return 0;

  allCounts.reserve( numDistinctKmers );

  for ( typename vec<Entry>::const_iterator entryIter = m_table.begin();
        entryIter != m_table.end(); ++entryIter )
    if ( entryIter->count > 1 )
      allCounts.push_back( entryIter->count );
  
  vec<unsigned char>::iterator median = allCounts.begin() + allCounts.size()/2;
  nth_element( allCounts.begin(), median, allCounts.end() );

  return *median;
}  


// Documented in header.
template <int K>
void KmerHashHistogram<K>::DoFillKmerCountHistogram( histogram<int>& theHistogram ) const
{
  for ( typename vec<Entry>::const_iterator entryIter = m_table.begin();
        entryIter != m_table.end(); ++entryIter )
    if ( entryIter->count > 0 )
      theHistogram.AddDatum( entryIter->count );
}  


// Documented in header.
KmerHistogram* GetNewKmerHistogram( const int K, const size_t capacity )
{
  switch ( K ) {
    #define CASE(kmersize) case kmersize: return new KmerHashHistogram<kmersize>(capacity); break
    CASE(16);
    CASE(24);
    CASE(32);
    CASE(48);
    default: 
      cout << "The class KmerHashHistogram<" << K << "> has not been implemented." << endl;
      TracebackThisProcess();
  }
  
  return 0;
}

#define INSTANTIATE(K) \
template void KmerHashHistogram<K>::DoCountKmers( const basevector& bases ); \
template void KmerHashHistogram<K>::DoGetCounts( const basevector& bases, \
                                                 vec<int>& counts ) const; \
template void KmerHashHistogram<K>::DoGetSomeCounts( const basevector& bases, \
                                                     vec<int>& counts ) const; \
template int KmerHashHistogram<K>::DoGetMedianNonuniqueKmerCount() const

INSTANTIATE(16);
INSTANTIATE(24);
INSTANTIATE(32);
INSTANTIATE(48);

