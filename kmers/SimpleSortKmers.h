/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef SIMPLE_SORT_KMERS_H
#define SIMPLE_SORT_KMERS_H

#include "Basevector.h"
#include "Vec.h"
#include "kmers/KmerRecord.h"
#include "AminoAcidConverter.h"


///Sort kmers in one basevector without dealing with reverse complements.
/// Set I to 2 to have 4-byte position values for the kmers.
/// All readIds are set to 0 in the kmer records, since there is only
/// one basevector.
template<int K, int I>
inline void SimpleSortKmers(const basevector & bases,  
			    vec<kmer_record<K,I> > & kmers) {
  const int LAST_KMER = bases.size() - K + 1;
  if (LAST_KMER <= 0) {
    kmers.clear();
    return;
  }
  kmers.resize(LAST_KMER);
  basevector b;
  for (int i = 0; i != LAST_KMER; ++i) {
    b.SetToSubOf(bases, i, K, 0); //get the next basevector
    kmers[i].Set(b, 0, i);
  }
  sort(kmers.begin(), kmers.end());
}


/// Sort kmers in a vecbasevector, ignoring reverse complements.
/// Set I to 2 to have 4-byte position values for the kmers.
/// The kmer_record id values are set to the index of the basevector the 
/// kmer came from.
/// \param selection which basevectors from bases to use. If this is 
/// empty (default), then use all of the basevectors in bases. The idea 
/// is to 
/// provide support for sparse vecbasevectors, such as ones made from
/// a bunch of reads.
template<int K, int I>
void SimpleSortKmers(const vecbasevector & bases, 
		     vec<kmer_record<K,I> > & kmers,		     
                     const vec<int> & selection = vec<int>(),
		     bool bUseAminoAcids = false,
		     bool bCanonicalize = false) {

  if (selection.empty()) { //accept all elements of bases.
    vec<int> & temp=const_cast<vec<int> &>(selection);//const cheat...
    for (size_t i = 0; i != bases.size(); ++i) temp.push_back(i);
  }
  
  const int BASEVECTORS = selection.size();
  unsigned int estimateKmers = 0;
  for (int i = 0; i != BASEVECTORS; ++i) {
    int toadd = bases[selection[i]].size() - K + 1;
    if (toadd > 0) {
      estimateKmers += toadd;
    }
  }

  //cout << "Allocating " << estimateKmers << " kmers @ " << sizeof(kmer_record<K,I>) << " bytes each" << endl;

  kmers.resize(estimateKmers);

  unsigned int currentKmer = 0;
  basevector b;
  //get all the kmers for all the basevectors and load them into kmers.
  for (int j=0; j != BASEVECTORS; ++j) {
    const basevector & bvec = bases[selection[j]];
    if (bvec.size() < K) continue;
    const int LAST_KMER = bvec.size() - K + 1;
    for (int i = 0; i !=  LAST_KMER; ++i) {
      b.SetToSubOf(bvec, i, K, 0); //get the next basevector

      if (bUseAminoAcids) {
	ReduceToAminoAcids(b);
      }

      kmers[currentKmer].Set(b, j, i);
      if (bCanonicalize) {
	kmers[currentKmer].Canonicalize(bvec.size());
      }
      ++currentKmer;
    }
  }
  //Make sure our estimate was correct!
  ForceAssertEq(currentKmer, kmers.size() );
  sort(kmers.begin(), kmers.end());
}                        

#endif
