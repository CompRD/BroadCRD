/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "kmer_freq/WriteKmerFrequencies.h"
#include "Vec.h"
#include "kmers/KmerRecord.h"
#include "kmers/SortKmers.h"
#include "feudal/BinaryStream.h"

// Inefficiency: WriteKmereFrequencies inefficiency
// The use of class <kmer> to store intermediate results and class <kmer_with_count>
// to store final answers is somewhat inefficient since both classes occupy 4n
// bytes for some n, even if the actual space requirement is not divisible by 4.


template <class KSHAPE, class KSELECTOR>
void
WriteKmerFrequencies( const vecbasevector& seqs,  const String& filename,
                      const bool include_unique, const KSELECTOR& kselector)
{
  const int K = KSHAPE::KSIZE;
  BinaryIteratingWriter< vec<kmer_with_count<K> > > writer( filename );
  int N = seqs.size( );
  vec<int> rid(N);
  for ( int i = 0; i < N; i++ )
    rid[i] = i;
  unsigned int S = 0;
  for ( size_t l = 0; l < seqs.size( ); l++ )
    S += seqs[l].size( ) - K + 1;
  S += S/4;
  S /= 33;
  typedef typename KSELECTOR::record_t record_t;
  vec< record_t > R(S);
  cout << "pass " << flush;
  for ( int pass = 0; pass < 100; pass++ ) {
    dummy<100> d100;
    SortKmers< KSHAPE, record_t >( d100, seqs, rid, pass, R, S );
    basevector theKmer(K);
    for ( int i = 0; i < (int) S; i++ ) {
      int j;
      for ( j = i+1; j < (int) S; j++ )
        if ( R[i] < R[j] ) break;
      unsigned short count = kselector(R, i, j);
      if  ((include_unique && count == 1) || count > 1) {
        R[i].GetBasevector(theKmer);
        writer.write( kmer_with_count<K>( theKmer, count ) );
      }
      i = j - 1;    
    }
    Dot( cout, pass );
  }
  cout << endl;
  writer.close();

  vec< kmer_with_count<K> > kmers;
  BinaryReader::readFile( filename, &kmers );
  cout << "Sorting... " << flush;
  Sort(kmers);
  cout << "done." << endl;
  BinaryWriter::writeFile( filename, kmers );
  PRINT2( filename, kmers.size() );
}    

/**
   Local class: DefaultKSelector

   Defines a <KSelector> that allows all kmers through.
*/
template <int K>
class DefaultKSelector {
public:
  typedef kmer<K> record_t;  // doesn't really matter which type, since we don't really look 

  unsigned short operator() (const vec< record_t >&, int fromKmerIdx, int toKmerIdx) const {
    return Min( toKmerIdx - fromKmerIdx, 65535 );
  }
};

/**
    Local class: TrustedKSelector

    Defines a <KSelector> that allows through only those kmers that have at least one
    _trusted occurrence_, that is, an occurrence all of whose bases are <trusted bases>.
    It can also be configured to require at least one trusted occurrence of the kmer
    in each direction (that is, to get one trusted read of the kmer from each strand).
    Returns the count of kmer occurances (trusted and untrusted) if there is at least one
    trusted occurance, otherwise returns 0.
*/
template <class KSHAPE>
class TrustedKSelector {
public:
  typedef kmer_record<KSHAPE::KSIZE> record_t;
  
  TrustedKSelector(const vecbitvector& trusted,
		   int threshold = 1,
		   bool requireTrustedBothDirs = false) : trusted_(trusted),
							  threshold_(threshold),
							  requireTrustedBothDirs_(requireTrustedBothDirs) {}
  
  unsigned short operator() (const vec< record_t >& R, int fromKmerIdx, int toKmerIdx) const {
    bool foundTrustedFwd = false, foundTrustedBwd = false;

    int trusted_occurrences = 0;
    for (int i = fromKmerIdx; i < toKmerIdx; ++i) {
      int read_id = R[i].GetId();
      int pos = R[i].GetPos();
      int origPos = pos;
      if (pos < 0) pos = -pos;
      Bool isTrusted = True;
      for (int bit = 0; bit < KSHAPE::KSIZE; ++bit) {
	if (!trusted_[read_id][ (pos-1) + KSHAPE::getShapeOffset(bit) ]) {
	  isTrusted = False;
	  break;
	}
      }
      if (isTrusted) {
	++trusted_occurrences;
	if (origPos > 0)
	  foundTrustedFwd = true;
	else
	  foundTrustedBwd = true;
	
	if (trusted_occurrences >= threshold_  &&
	    (!requireTrustedBothDirs_ ||  (foundTrustedFwd && foundTrustedBwd)))
	    return Min( toKmerIdx - fromKmerIdx, 65535 );
      }
    }
    return 0;
  }
  
private:
  const vecbitvector& trusted_;
  const int threshold_;
  const bool requireTrustedBothDirs_;
  
};  // class TrustedKSelector


/**
    Local class: TrustedCountKSelector

    Defines a <KSelector> that allows through only those kmers that have at least one
    _trusted occurrence_, that is, an occurrence all of whose bases are <trusted bases>.
    Returns the count of trusted occurances.
*/
template <class KSHAPE>
class TrustedCountKSelector {
public:
  typedef kmer_record<KSHAPE::KSIZE> record_t;
  
  TrustedCountKSelector(const vecbitvector& trusted) : trusted_(trusted) {}
  
  unsigned short operator() (const vec< record_t >& R, int fromKmerIdx, int toKmerIdx) const {

    int trusted_occurrences = 0;
    for (int i = fromKmerIdx; i < toKmerIdx; ++i) {
      int read_id = R[i].GetId();
      int pos = R[i].GetPos();
      if (pos < 0) pos = -pos;
      Bool isTrusted = True;
      for (int bit = 0; bit < KSHAPE::KSIZE; ++bit) {
	if (!trusted_[read_id][ (pos-1) + KSHAPE::getShapeOffset(bit) ]) {
	  isTrusted = False;
	  break;
	}
      }
      if (isTrusted)
	++trusted_occurrences;
    }
    return trusted_occurrences;
  }
  
private:
  const vecbitvector& trusted_;
  
};  // class TrustedCountKSelector


     
template <class KSHAPE>
void
WriteKmerFrequencies( const vecbasevector& seqs,  const String& filename, 
                      bool include_unique ) {
  typedef DefaultKSelector< KSHAPE::KSIZE > ksel_t;
  WriteKmerFrequencies<KSHAPE, ksel_t > (seqs, filename, include_unique, ksel_t());
}


template <class KSHAPE>
void
WriteKmerFrequencies( const vecbasevector& seqs, const String& filename,
		      const vecbitvector& trusted, int threshold,
		      bool requireTrustedBothDirs, 
		      bool include_unique ) {
  typedef TrustedKSelector< KSHAPE > ksel_t;
  WriteKmerFrequencies<KSHAPE, ksel_t > (seqs, filename, include_unique,
					 ksel_t(trusted, threshold, requireTrustedBothDirs));
}

template <class KSHAPE>
void
WriteKmerTrustedFrequencies( const vecbasevector& seqs, const String& filename,
			     const vecbitvector& trusted) {
  typedef TrustedCountKSelector< KSHAPE > ksel_t;
  WriteKmerFrequencies<KSHAPE, ksel_t > (seqs, filename, true, ksel_t(trusted));
}


#define INSTANTIATE(KSHAPE, dummy) \
 template void WriteKmerFrequencies<KSHAPE>( const vecbasevector&, const String& kfreqOutputFile, const bool include_unique ); \
 template void WriteKmerFrequencies<KSHAPE>( const vecbasevector&, const String& kfreqOutputFile, \
					     const vecbitvector& trusted, int threshold, bool requiredTrustedBothDirs, const bool include_unique ); \
 template void WriteKmerTrustedFrequencies<KSHAPE>( const vecbasevector&, \
                                                    const String& kfreqOutputFile, \
					            const vecbitvector& trusted ); \
 template void WriteKmerFrequencies<KSHAPE, DefaultKSelector< KSHAPE::KSIZE>  >( const vecbasevector&, const String&, const bool,  const DefaultKSelector< KSHAPE::KSIZE >&); \
 template void WriteKmerFrequencies<KSHAPE, TrustedKSelector< KSHAPE > >( const vecbasevector&, const String&, const bool,  const TrustedKSelector< KSHAPE >&); \
 template void WriteKmerFrequencies<KSHAPE, TrustedCountKSelector< KSHAPE > >( const vecbasevector&, const String&, const bool,  const TrustedCountKSelector< KSHAPE >&)


FOR_ALL_KSHAPES(INSTANTIATE,unused);



