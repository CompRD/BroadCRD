///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/** Do a quick alignment of a read to a basevector using kmers.

\class KmerAligner

This class contains an internal basevector or vecbasevector set with 
SetBases. You can 
then do a quick alignment of any basevector or read to this internal
vec/basevector using FindPossibleAlignments.

Here is the way the alignment works: we cut up the read into kmers
separated by m_kmerStep bases. If K is big, this is 16, otherwise it is
 (K/2) % 4. It can be seen/changed with G/SetKmerStep.

Then we try to match any of these kmers to any of 
the kmers in the internal basevector. We save all the positions for all
the matches we find, and return those in a vec<int> or a vec<PerfectMatch>.

The class needs an Assembly so that it can appropriately load the bases
from a ReadLocation. Note that if you are going to align more than a 
few reads, I have had fastest results when setting PrefetchAllStrategy
for both reads and contigs in the assembly.

The template parameter K is the number of bases in the kmers we are using
to align things. Note that we only use kmer_records with I=2 because
we want to be sure we have 4 bytes for the kmer position!
    
*/

#ifndef KMER_ALIGNER_H
#define KMER_ALIGNER_H


#include "CoreTools.h"
#include "TaskTimer.h"
#include "assembly/Assembly.h"
#include "assembly/PrefetchStrategy.h"
#include "pairwise_aligners/PerfectMatch.h"
#include "kmers/SimpleSortKmers.h"
#include "ReadLocation.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/LocalAlign.h"
#include "RemoveUnique.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"

#include <algorithm>

template<int K=24>
class KmerAligner {
  friend class KmerAlignerTester;

 public:
  typedef kmer_record<K,2> krec;
  typedef vec<krec> container_type;
  typedef typename container_type::iterator kiter;
  static const int MIN_KMER_STEP=4;
  static const int MIN_K=4;

  KmerAligner(const Assembly & theAssembly = Assembly());
  ~KmerAligner();

  ///Set the offset between kmers tested, make sure it's big enough.
  void SetKmerStep(int step) { m_kmerStep = max(step, MIN_KMER_STEP); }

  int GetKmerStep() { return m_kmerStep; }

  ///This sets the contig's consensus as MyBases
  void SetBases(const Contig & contig) {
    SetBases(contig.GetBases());
  }
 
  ///align against this basevector.
  ///This will sort the basevector's kmers and save the
  ///sorted information.
  void SetBases(const basevector & bases) {
    SimpleSortKmers(bases, m_kmers);
  }

  ///align against this vecbasevector.
  ///This will sort the vecbasevector's kmers and save the
  ///sorted information.
  /// \param selection which basevectors from bases to use. If this is 
  /// empty (default), then use all of the basevectors in bases. 
  /// The idea is to 
  /// provide support for sparse vecbasevectors, such as ones made from
  /// a bunch of reads.
  void SetBases(const vecbasevector & bases, 
                const vec<int> & selection = vec<int>()) {
    SimpleSortKmers(bases, m_kmers, selection);
  }

  void SetBases(const vecbasevector & bases, 
                bool bConvertToAminoAcids,
		const vec<int> & selection = vec<int>(),
		bool bCanonicalize = false) {
		  
    SimpleSortKmers(bases, m_kmers, selection, bConvertToAminoAcids, 
		    bCanonicalize);
  }
 
 
  longlong KmerSize() const {
    longlong ret = m_kmers.size();
    return ret * sizeof(krec);
  }

  ///Canonicalize all our kmers and re-sort.
  void Canonicalize(const vec<unsigned int> & lengths) {
    for (kiter i= m_kmers.begin(); i != m_kmers.end(); ++i) {
      i->Canonicalize(lengths[i->GetId()]);
    }
    sort(m_kmers.begin(), m_kmers.end());
  }
 
  /// Remove all unique kmers from the kmer table.
  void RemoveUnique();

  ///Helper metod, calls FindPossibleAlignments(basevector...).
  void FindPossibleAlignments(const ReadLocation & loc, vec<int> & positions);
     
  ///Helper metod, calls FindPossibleAlignments(basevector...).
  void FindPossibleAlignments(const ReadLocation & loc, vec<int> & positions, 
                              vec<int> & ids);

  /** Returns all possible alignments found with kmers of size K
      staggered 16 bases apart (for efficiency, we do not check every kmer)
      of the basevector to our internal basevector
      in the unsorted vec positions. The vector contains all possible 
      alignments found with kmers, so it may contain multiple copies of
      the same position in the contig if bases aligns there based on
      more than one of its kmers.
  */
  void FindPossibleAlignments(const basevector &, vec<int> & positions);


  /** Returns all possible alignments found with kmers of size K
      staggered 16 bases apart (for efficiency, we do not check every kmer)
      of the basevector to our internal vecbasevector
      in the unsorted vec positions. The vector contains all possible 
      alignments found with kmers, so it may contain multiple copies of
      the same position in the contig if bases aligns there based on
      more than one of its kmers.
      \param basevector the basevector to match
      \param positions list of possible positions within a basevector, 
      \param ids list, parallel to positions, of
             the id of the basevector that matched, i.e. its 
             position in the vecbasevector that was set with SetBases.
      Note that if SetBases was called with a Contig or a simple
      basevector, the call will still work correctly and all ids will be 0. 
  */
  void FindPossibleAlignments(const basevector &, 
                              vec<int>& positions, vec<int> & ids);


  
  /** Returns all possible alignments found with kmers of size K
      staggered 16 bases apart (for efficiency, we do not check every kmer)
      of the basevector to our internal vecbasevector
      in the unsorted vec aligns. The vector contains all possible 
      alignments found with kmers, so it may contain multiple copies of
      the same position in the contig if bases aligns there based on
      more than one of its kmers.
      \param basevector the basevector to match
      \param aligns list of possible positions within a basevector, given
as PerfectMatches (i.e. including ids, start positions, etc...). 
      \param queryId identifier for the basevector being passed in.
      Note that if SetBases was called with a Contig or a simple
      basevector, the call will still work correctly and all target 
      ids will be 0. 
  */
  void FindPossibleAlignments(const basevector & query, 
                              vec<PerfectMatch> & aligns,
                              int queryId = 0);

  /// As above, but will also find alignemnts with reverse complement of query.
  void FindPossibleAlignmentsWithOrientation(const basevector & query, 
                                             vec<PerfectMatch> & aligns,
                                             int queryId = 0) {
    FindPossibleAlignments(query, aligns, queryId);
    vec<PerfectMatch> reverseAligns;
    basevector reversed;
    reversed.ReverseComplement(query);
    FindPossibleAlignments(reversed, reverseAligns, queryId);
    for (int i = 0; i != reverseAligns.isize(); ++i) {
      reverseAligns[i].rc1 = true;
      aligns.push_back(reverseAligns[i]);
    }
  }

  ///As above, but takes a vector of queries.
  void FindPossibleAlignments(const vecbasevector & queries, 
                              vec<PerfectMatch> & aligns) {
    aligns.clear();
    vec<PerfectMatch> temp;
    for (int i = 0; i != queries.size(); ++i) {
      FindPossibleAlignments(queries[i], temp, i);
      copy(temp.begin(), temp.end(), back_inserter(aligns));
    }
  }
      
  ///As above, checks many queries and their reverse complements.
  void FindPossibleAlignmentsWithOrientation(const vecbasevector & queries, 
                                             vec<PerfectMatch> & aligns) {
    aligns.clear();
    vec<PerfectMatch> temp;
    for (int i = 0; i != queries.size(); ++i) {
      FindPossibleAlignmentsWithOrientation(queries[i], temp, i);
      copy(temp.begin(), temp.end(), back_inserter(aligns));
    }
  } 

  /// This method does not work right yet!
  /// Only tries for fw alignments, and the call to SmithWatBandedA
  /// sometimes gives really ugly alignments, so it needs to be tuned.
  /// Given a basevector and a kmeraligner, find seeds, group them, and
  /// extend them using SmithWatBandedA, and save in aligns.
  template<typename AlignVec> 
  void FindExtendedSmithWatAlignments(const basevector & bases,
				      const int query_id,
				      const vecbasevector & ref,
				      AlignVec & aligns,  
				      const int BAND = 2);

  /// Economically find all kmer matches with our own bases.
  ///
  /// Templatized with an input iterator so we can send this into a
  /// vector or into a file to save memory.
  ///
  /// RemoveUnique should be called first for efficiency, but it is not
  /// necessary.
  ///
  /// This will work correctly with kmer_records on which 
  /// this->Canonicalize has been called. Note that this is different
  /// from canonicalizing the basevector and then putting in the kmer_record,
  /// because the position will not be made negative in that case.
  /// 
  /// This method works fine in terms of time, but the resulting vec of 
  /// PerfectMatches is huge, and needs to be brought into memory to be 
  /// sorted. So it's not efficient enough in the end, which is why 
  /// FindAndExtendSelfPerfectMatchesById was written: it does the work
  /// separately for each sequence in our original vecbasevector, and
  /// therefore ends up with much more tractable sizes for the vecs
  /// of PerfectMatches that need to be sorted for extension/compaction.
  template<class InputIter>
  void FindSelfPerfectMatches(InputIter in, const vec<unsigned int> & lengths);

  /// Economically find all kmer matches with our own bases.
  ///
  /// Templatized with an input iterator so we can send this into a
  /// vector or into a file to save memory.
  ///
  /// RemoveUnique should be called first for efficiency, but it is not
  /// strictly necessary.
  ///
  /// This will work correctly with kmer_records on which 
  /// this->Canonicalize has been called. Note that this is different
  /// from canonicalizing the basevector and then putting in the kmer_record,
  /// because the position will not be made negative in that case. 
  ///
  /// This method does the work
  /// separately for each sequence in our original vecbasevector, and
  /// therefore ends up with much more tractable sizes for the vecs
  /// of PerfectMatches that need to be sorted for extension/compaction.
  ///
  /// The trick is to sort the kmers twice: they are sorted by kmer
  /// in m_kmers, and we do an indirect sort on a vector of indices
  /// that ends up sorted by id and position instead.
  template<class InputIter>
  void FindAndExtendSelfPerfectMatchesById(InputIter in, 
					   const vecbasevector & bases,
					   int minKmerDistance = 0);

  ///Write the internal kmer table to file for later reuse.
  void WriteBasesToFile(const String & filename);

  ///Read the internal kmer table from file, use instead of SetBases().
  void ReadBasesFromFile(const String & filename);

  ///Print the bases to an ostream for debugging, one per line.
  void Print(ostream & out) const;

  // Give access to the kmer table.
  const container_type & GetKmerTable() const {return m_kmers;} 

  // How big is the table?
  int GetKmerTableSize() const {return (int)m_kmers.size();} 
 private:

  const Assembly & m_assembly;
  container_type m_kmers;///<We save the sorted kmers for reuse.
  int m_kmerStep;

  ///Private: disallow copy.
  KmerAligner(const KmerAligner &);
  ///Private: disallow assignment.
  KmerAligner & operator=(const KmerAligner &);

  ///Find all perfect matches between m_kmers[pos] and all our other
  ///kmers, and put them into in.
  ///Precondition: our m_kmers are sorted by kmer (as they always should be!).
  template<typename InputIter>
  void AddSelfMatches(int pos, InputIter in,
		      const vec<unsigned int> & lengths) const;

  
  /// Return a perfect match between these two kmers.
  /// Because we know they are our own kmers and belong to the same
  /// sequence, we only need one lengths vector, and we can also switch
  /// the order if it is convenient for us.
  PerfectMatch SelfMatch(const krec & first,
			 const krec & second, 
			 const vec<unsigned int> & lengths) const;
};

template<int K>
KmerAligner<K>::KmerAligner(const Assembly & assembly): 
  m_assembly(assembly),
  m_kmers(),
  m_kmerStep(16) //original fixed value.
{
  STATIC_ASSERT_M(K>=MIN_K, kmer_size_too_small);
  // if small K, set the step to be smaller than K, 
  // and measured in bytes (4bases=1byte).
  if (K<16) {
    SetKmerStep(K/2 - (K/2% 4));
  }
}

template<int K>
KmerAligner<K>::~KmerAligner() {}

template<int K>
const int KmerAligner<K>::MIN_KMER_STEP;

template<int K>
const int KmerAligner<K>::MIN_K;

template<int K>
void KmerAligner<K>::RemoveUnique() {
  CompareForwardKmers< kmer_record<K,2> > eq;
  kiter newend = ::RemoveUnique(m_kmers.begin(), m_kmers.end(), eq);
  m_kmers.erase(newend, m_kmers.end());
}

template<int K>
void KmerAligner<K>::FindPossibleAlignments(const ReadLocation & loc, 
                                         vec<int> & positions) {
  vec<int> ids;
  FindPossibleAlignments(loc, positions, ids);
}

template<int K>
void KmerAligner<K>::FindPossibleAlignments(const ReadLocation & loc, 
                                            vec<int> & positions, 
                                            vec<int> & ids) {

  int id = loc.GetRead().GetId();
  Read read = m_assembly.GetRead(id);
  if (loc.IsForward()) {
    FindPossibleAlignments(read.GetBases(), positions, ids);
  } else {
    basevector reversed = read.GetBases();
    reversed.ReverseComplement();
    FindPossibleAlignments(reversed, positions, ids);
  }
  
}

template<int K>
void KmerAligner<K>::FindPossibleAlignments(const basevector & bases, 
                                            vec<int> & positions) {
  vec<int> ids;
  FindPossibleAlignments( bases, positions, ids);
}

template<int K>
void KmerAligner<K>::FindPossibleAlignments(const basevector & bases, 
                            vec<int >& positions,
                            vec<int> & ids) {

  positions.clear();
  ids.clear();
  //Start off at a reasonable size.
  positions.reserve(1000);
  ids.reserve(1000);
  basevector kmerbases;

  const int LAST_KMER = bases.size() - K;
  kmer_record<K,2> kmer;
  pair<kiter, kiter> range;
  for (int i = 0; i <= LAST_KMER; i += m_kmerStep) {
    //we move forward in groups of m_kmerStep. This may be less
    //efficient than always setting the step to 16 because SetToSubOf() works
    //fastest if it does not have do do any shifting.
    kmerbases.SetToSubOf(bases, i, K);
    kmer.Set(kmerbases, 0, 0);//we don't use the location or readId
    range = equal_range(m_kmers.begin(), m_kmers.end(), kmer);
    int rangeSize = range.second - range.first;
    if (0 == rangeSize  ) { //not found, or duplicated.
      continue;
    } else {
      for (kiter iter = range.first; iter != range.second; ++iter) {
        positions.push_back(iter->GetPos() - i);
        ids.push_back(iter->GetId());
      }
    }
  }
}  

template<int K>
void KmerAligner<K>::FindPossibleAlignments(const basevector & bases, 
                                            vec<PerfectMatch >& positions,
                                            int queryId){

  positions.clear();
  //Start off at a reasonable size.
  positions.reserve(1000);
  basevector kmerbases;

  const int LAST_KMER = bases.size() - K;
  kmer_record<K,2> kmer;
  pair<kiter, kiter> range;
  for (int i = 0; i <= LAST_KMER; i += 16) {
    //we move forward in groups of 16 bases = 1 int.
    //This is the most efficient way to do it because SetToSubOf() works
    //fastest if it does not have do do any shifting.
    kmerbases.SetToSubOf(bases, i, K);
    kmer.Set(kmerbases, 0, 0);// we don't use location or id here
    range = equal_range(m_kmers.begin(), m_kmers.end(), kmer);
    int rangeSize = range.second - range.first;
    if (0 == rangeSize  ) { //not found, or duplicated.
      continue;
    } else {
      for (kiter iter = range.first; iter != range.second; ++iter) {
        positions.push_back
          (PerfectMatch(queryId, i, iter->GetId(), iter->GetPos(), K, false));
      }
    }
  }
}  
 
template<int K>
void KmerAligner<K>::WriteBasesToFile(const String & filename) {
  BinaryWriter::writeFile(filename, m_kmers);
}

template<int K>
void KmerAligner<K>::ReadBasesFromFile(const String & filename) {
  BinaryReader::readFile(filename, &m_kmers);
}
 
template<int K>
void KmerAligner<K>::Print(ostream & os) const {
  const int S = m_kmers.size();
  for (int i = 0; i != S; ++i) {
    os << m_kmers[i];
  }
}

template<int K>
ostream & operator<<(ostream & os, const KmerAligner<K> & aligner) {
  aligner.Print(os);
  return os;
}

/* 
algorithm;
for each source
  for each target
    find all forward matches on that target
   sort matches by offset and then start
    for each offset
         extend the first match and mark as last extended
          for each other match
               if included in last extended mark as bad
               else extend
  remove all bad matches
  add remaining matches to total matches and adjust their indices.
*/

template<int K> template<typename InputIter>
void KmerAligner<K>::FindSelfPerfectMatches(InputIter in,
					    const vec<unsigned int> & lengths) {
  kiter start = m_kmers.begin();
  kiter last = start;
  while ( start != m_kmers.end() ) {
    //find range for all equal kmers.
    while (last != m_kmers.end() && last->EqualKmers(*start)) ++last; 
    //find all N(N+1)/2 pairs and make them into PerfectMatches.
    for (kiter first = start; first != last; ++first) { 
      kiter second = first;
      for (++second; second != last; ++second) {
	int pos1, pos2;
	Bool r;
	*in = SelfMatch(*first, *second, lengths);	
	++in;
      } //end of second loop
    } //end of first loop
    start = last; //go on to the next group of equal kmers
  } //end of start loop
}


template<int K> template<typename OutputIter>
void KmerAligner<K>::FindAndExtendSelfPerfectMatchesById
(OutputIter out, const vecbasevector & bases, int minKmerDistance) {

  vec<unsigned int> lengths;
  bases.ElementSizes(lengths);
  
  vec<int> byId(m_kmers.size());
  for (unsigned int i=0; i != m_kmers.size(); ++i) byId[i]=i;
  indirect_compare<container_type,LessByIdAndPos<krec> > indirById(m_kmers);
  sort(byId.begin(), byId.end(), indirById);
  //PRINT(byId);
  vec<int>::iterator current, start, last;
  vec<PerfectMatch> matches;
  cout << "Analyzing " << bases.size() << " contigs in groups of 1000\n";
  int contig=0;
  for (start = byId.begin(); start != byId.end(); start = current) {
    //find all the matches for one contig;
    matches.clear();
    AddSelfMatches(*start,back_inserter(matches), lengths);
    last = start;
    for (current = start+1; 
	 current != byId.end() && 
	   m_kmers[*current].GetId() == m_kmers[*start].GetId(); ++current) {
      //ignore if too close and same orientation to previous kmer.
      if (0 == minKmerDistance  //check for default first to save time.
	  || ( (m_kmers[*current].GetPos() - m_kmers[*last].GetPos()) 
	       > minKmerDistance )
	  || (m_kmers[*current].IsReversed() != m_kmers[*last].IsReversed() )
	  ) {
	AddSelfMatches(*current,back_inserter(matches), lengths);
	last = current;
      }
    }

    vec<PerfectMatch>::iterator newend = 
      ExtendAndRemoveOverlaps(matches.begin(), matches.end(), bases, bases);
    copy(matches.begin(), newend, out);
    DotMod(cout, ++contig, 1000);
    //don't bother erasing, the clear() at top takes care of that.
  }
}
    
template<int K> template<typename InputIter>
void KmerAligner<K>::AddSelfMatches(int pos, InputIter in,
                    const vec<unsigned int> & lengths) const {
  const krec & k = m_kmers[pos];
  int id = k.GetId();
  for (int i=pos-1; i >= 0 && k.EqualKmers(m_kmers[i]); --i) {
    if (id < m_kmers[i].GetId()) continue;
    *in = SelfMatch(m_kmers[pos], m_kmers[i], lengths);
    ++in;
  }
  for (unsigned int i=pos+1; 
       i < m_kmers.size() && k.EqualKmers(m_kmers[i]); ++i) {
    if (id < m_kmers[i].GetId()) continue;
    *in = SelfMatch(m_kmers[pos], m_kmers[i], lengths);
    ++in;
  }
}
  
/// Return a perfect match between these two kmers.
/// Because we know they are our own kmers and belong to the same
/// sequence, we only need one lengths vector.
/// We make sure to order them by id and start position, so that mirror-image
/// PerfectMatches will be extended and compacted together.
template<int K> 
PerfectMatch KmerAligner<K>::SelfMatch(const krec & first,
				       const krec & second, 
				       const vec<unsigned int> & lengths) const {
  //Make sure PerfectMatches are always created in the same order,
  //so we only get one copy of each after extension.
  const krec * firstp = &first;
  const krec * secondp = &second;
  if (first.GetId() > second.GetId()
      || ( first.GetId() == second.GetId() 
	   && first.GetPos() > second.GetPos() ) ) {
    swap(secondp, firstp);
  }

  //remember, the target of a PerfectMatch must always be forward!
  switch (firstp->IsReversed()) {
  case true: 
    switch (secondp->IsReversed()) {
    case true: //reverse both of them
      return PerfectMatch(firstp->GetId(), 
			  lengths[firstp->GetId()] - K - firstp->TruePos(), 
			  secondp->GetId(), 
			  lengths[secondp->GetId()]-K - secondp->TruePos(), 
			  K, false);
    case false:
      return PerfectMatch(firstp->GetId(), firstp->TruePos(),
			  secondp->GetId(), secondp->GetPos(),
			  K, true);
    }
  case false:
    switch (secondp->IsReversed()) {
    case true: //reverse both of them
      return PerfectMatch(firstp->GetId(), 
			  lengths[firstp->GetId()] - K - firstp->TruePos(),
			  secondp->GetId(), 
			  lengths[secondp->GetId()]- K - secondp->TruePos(), 
			  K, true);
    case false:
      return PerfectMatch(firstp->GetId(), firstp->GetPos(), 
			  secondp->GetId(), secondp->GetPos(), 
			  K, false);
    }
  }
  FatalErr("should never get here");//avoid dumb compiler warning.
}
  
  


/// Find all possible oriented alignments based on kmers of length K.
/// This global method saves you from having to handle a KmerAligner 
/// directly.
/// It also allows for the use of sparse vecbasevectors for both 
/// target and queries. That is, you can pass in large vecbasevectors
/// and only use the subset of them specified by targetIds and queryIds.
///
/// If the xxxIds vectors are empty, then the complete 
/// corresponding vecbasevector is used.
template<int K>
void FindPossibleOrientedAlignments(const vecbasevector & target,
                                    const vecbasevector & queries,
                                    vec<PerfectMatch> & aligns,
                                    const vec<int> & targetIds = vec<int>(),
                                    const vec<int> & queryIds = vec<int>()) {
  //we don't need to deal with the empty targetIds vector because
  //SimpleSortKmers knows how to handle that.
  if (queryIds.empty()) {
    vec<int> & temp=const_cast<vec<int> &>(queryIds);
    for (size_t i=0; i != queries.size(); ++i) temp.push_back(i);
  }
  KmerAligner<K> aligner;
  aligner.SetBases(target, targetIds);
  aligns.clear();
  vec<PerfectMatch> temp;
  vec<PerfectMatch>::iterator lastGood;
  const int QUERIES = queryIds.size();
  for (int i=0; i != QUERIES; ++i) {
    temp.clear();
    //find all kmer matches for this query
    aligner.FindPossibleAlignmentsWithOrientation
      (queries[queryIds[i]], temp, queryIds[i]);
    //compress all kmer matches that go together and extend them.
    lastGood = ExtendAndRemoveOverlaps(temp.begin(), temp.end(), 
	  queries, target);
    //add these to the overall set of matches.
    copy(temp.begin(), lastGood, back_inserter(aligns));
  }
}

///Wrapper around the templated version of the method for one target.
/// Allows KMER_SIZE to be determined at runtime from(8,12,16,24,32,48,64,96).
void FindPossibleOrientedAlignments
(const int KMER_SIZE, 
 const basevector & target, 
 const vecbasevector & queries,
 vec<PerfectMatch> & aligns,
 const vec<int> & queryIds = vec<int>() );

///Wrapper around the templated method for multiple targets.
/// Allows KMER_SIZE to be determined at runtime.
void FindPossibleOrientedAlignments
(const int KMER_SIZE, 
 const vecbasevector & target, 
 const vecbasevector & queries,
 vec<PerfectMatch> & aligns,
 const vec<int> & targetIds = vec<int>(),
 const vec<int> & queryIds = vec<int>() );



#include "lookup/LookAlign.h"
struct BadAlignForCompress : public unary_function<look_align,bool> {
  bool operator()(const look_align & a) { 
    return a.target_id == INT_MAX; 
  }
};

template<class LookAlignVec> 
void CompressAligns(LookAlignVec & aligns) {
  BadAlignForCompress badAlign;

  for (int i=0; i != aligns.isize(); ++i) {
    ForceAssertEq(aligns[i].query_id, aligns[0].query_id);
  }
  sort(aligns.begin(), aligns.end(), order_lookalign_TargetBegin());
  for (int i=0; i < aligns.isize() -1; ++i) {
    //if two alignments overlap, pick the better one.
    if (aligns[i].target_id == aligns[i+1].target_id &&
    	aligns[i].query_id == aligns[i+1].query_id &&
    	aligns[i].rc1 == aligns[i+1].rc1 &&
	0 == Distance(aligns[i].pos2(), aligns[i].Pos2(),
		      aligns[i+1].pos2(), aligns[i+1].Pos2()) &&
	0 == Distance(aligns[i].pos1(), aligns[i].Pos1(),
		      aligns[i+1].pos1(), aligns[i+1].Pos1())) {
      if (aligns[i].ErrorRate() < aligns[i+1].ErrorRate()) {
	swap(aligns[i], aligns[i+1]);
      }
      aligns[i].target_id = INT_MAX; //so we can tell it's bad and remove it.
    }
  }
  aligns.erase(remove_if(aligns.begin(),aligns.end(), badAlign), 
	       aligns.end()); 
}
	   
  

/// Given a basevector and a kmeraligner, find seeds, group them, and
/// extend them using SmithWatBandedA, and save in aligns (which is cleared
/// first).
template<int K> template <typename AlignVec> 
void KmerAligner<K>::FindExtendedSmithWatAlignments(const basevector & bases,
						    const int query_id,
						    const vecbasevector & ref,
						    AlignVec & aligns,  
						    const int BAND) {

  vec<PerfectMatch> matches;
  int ignoreErr;
  aligns.clear();
  vecbasevector source;
  source.resize(query_id);
  source.push_back(bases);
  
  FindPossibleAlignments( bases, matches, query_id );
  matches.erase(ExtendAndRemoveOverlaps(matches.begin(), matches.end(), 
	source, ref), matches.end());
  cout << "after extension:" << matches << endl;
     
  //Extend alignments with LocalAlign.
  for (int i=0; i< matches.isize( ); i++) {
    look_align la;
    basevector shortref;
    int MARGIN=100;
    int leftmargin = max(0,matches[i].startOnTarget - MARGIN);
    int rightmargin = min(int(ref[matches[i].target].size()),
	matches[i].EndOnTarget() + MARGIN);
    shortref.SetToSubOf(ref[matches[i].target], leftmargin, rightmargin - leftmargin);
    LocalAlign(bases, shortref, la.a, 1,3,2);

    //adjust the align to refer to the original reference, not the short form.
    la.a.Setpos2(leftmargin + la.a.pos2());

    vec<int> errs = 
      la.a.MutationsGap1Gap2( bases, ref[matches[i].target]);

    //Reject alignments with >14% errors.
    int totalErrs = errs[0]+errs[1]+errs[2];
    PRINT3(totalErrs, la.a.extent1(), i);
    if (totalErrs > (int(la.a.extent1()) / 7) ) continue;

    //save the alignment.
    la.target_id = matches[i].target;
    la.query_id = query_id;
    la.query_length = bases.size( );
    la.target_length = ref[la.target_id].size( );
    la.rc1 = false;
    la.mutations = errs[0];
    la.indels = errs[1] + errs[2];
    aligns.push_back(la);
  }
  //clean up the alignments.
  for(int i=0; i != aligns.isize(); ++i) aligns[i].PrintReadableBrief(cout);
  CompressAligns(aligns);
  cout << "after compression\n";
  for(int i=0; i != aligns.isize(); ++i) aligns[i].PrintReadableBrief(cout);
}

#endif //KMER_ALIGNER_H
