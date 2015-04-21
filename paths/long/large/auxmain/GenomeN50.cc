///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * GenomeN50.cc
 *
 *  Created on: Jun 24, 2014
 *      Author: tsharpe
 */
#include "MainTools.h"
#include "feudal/HashSet.h"
#include "kmers/BigKMer.h"
#include "kmers/KMerHasher.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"

namespace
{

template <unsigned BIGK>
using BigDict = HashSet<BigKMer<BIGK>,typename BigKMer<BIGK>::hasher>;


template <unsigned BIGK>
class BigKMerizer
{
public:
    typedef BigDict<BIGK> BigKDict;
    BigKMerizer( BigKDict* pDict ) : mDict(*pDict) {}

    void kmerize( bvec const& bv )
    { if ( bv.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      auto itr = bv.begin();
      size_t hash = hasher.hash(itr);
      if ( bv.size() == BIGK ) { add(BigKMer<BIGK>(bv,hash)); return; }
      BigKMer<BIGK> kmer(bv,hash,KMerContext::initialContext(bv[BIGK]));
      add(kmer);
      auto end = bv.end()-BIGK;
      while ( ++itr != end )
      { hash = hasher.stepF(itr);
        kmer.successor(hash,KMerContext(itr[-1],itr[BIGK]));
        add(kmer); }
      hash = hasher.stepF(itr);
      kmer.successor(hash,KMerContext::finalContext(itr[-1]));
      add(kmer); }

private:
    void add( BigKMer<BIGK> const& kmer )
    { canonicalAdd( kmer.isRev() ? kmer.rc() : kmer ); }

    void canonicalAdd( BigKMer<BIGK> const& kmer ) const
    { mDict.apply(kmer,
        [&kmer]( BigKMer<BIGK> const& entry )
        { const_cast<BigKMer<BIGK>&>(entry).addContext(kmer.getContext()); }); }

    BigKDict& mDict;
};

template <unsigned BIGK>
class BigKEdgeBuilder
{
public:
    typedef BigKMer<BIGK> BigKmer;
    typedef BigKMerizer<BIGK> BigKmerizer;
    typedef BigDict<BIGK> BigKDict;

    static void buildEdges( BigKDict const& dict, vecbvec* pEdges )
    {
        BigKEdgeBuilder eb(dict,pEdges);
        dict.parallelForEachHHS(
                [eb]( typename BigKDict::HHS const& hhs ) mutable
                { for ( BigKmer const& entry : hhs )
                    if ( entry.isUnassigned() )
                      eb.buildEdge(entry); });

        size_t nRegularEdges = pEdges->size();
        size_t nTotalLength = pEdges->SizeSum();

        // if a kmer isn't marked as being on an edge, it must be a part of a
        // smooth circle: add those edges, too.  simpleCircle method isn't
        // thread-safe, so this part is single-threaded.
        for ( auto const& hhs : dict )
            for ( auto const& entry : hhs )
                if ( entry.isUnassigned() )
                    eb.simpleCircle(entry);
    }

private:
    BigKEdgeBuilder( BigKDict const& dict, vecbvec* pEdges )
    : mDict(dict), mEdges(*pEdges) {}

    void buildEdge( BigKmer const& entry )
    { if ( isPalindrome(entry) )
        make1KmerEdge(entry);
      else if ( upstreamExtensionPossible(entry) )
      { if ( downstreamExtensionPossible(entry) )
          return;
        BigKmer rc = entry.rc();
        mEdgeSeq.assign(rc.begin(),rc.end());
        mEdgeEntries.push_back(&entry);
        extend(rc.getContext()); }
      else if ( downstreamExtensionPossible(entry) )
      { mEdgeSeq.assign(entry.begin(),entry.end());
        mEdgeEntries.push_back(&entry);
        extend(entry.getContext()); }
      else
        make1KmerEdge(entry); }

    // not thread-safe
    void simpleCircle( BigKmer const& entry )
    { BigKmer const* pFirstEntry = &entry;
      mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(pFirstEntry);
      auto itr = mEdgeSeq.begin();
      KMerHasher<BIGK> hasher;
      BigKmer kmer(mEdgeSeq,hasher.hash(itr));
      KMerContext context = entry.getContext();
      while ( true )
      { ForceAssertEq(context.getPredecessorCount(),1u);
        ForceAssertEq(context.getSuccessorCount(),1u);
        mEdgeSeq.push_back(context.getSingleSuccessor());
        kmer.successor(hasher.stepF(++itr));
        BigKmer const* pEntry = lookup(kmer,&context);
        if ( pEntry == pFirstEntry )
        { mEdgeSeq.pop_back();
          break; }
        if ( !pEntry->isUnassigned() )
        { std::cout << "Failed to close circle.\n";
          for ( auto beg=mEdgeEntries.begin(),end=mEdgeEntries.end();
                  beg != end; ++beg )
            std::cout << *beg << ' ' << **beg;
          std::cout << pEntry << ' ' << *pEntry << std::endl;
          CRD::exit(1); }
        mEdgeEntries.push_back(pEntry); }
      addEdge(); }

    bool isPalindrome( BigKmer const& kmer )
    { auto itr = kmer.begin();
      if ( !(BIGK&1) )
        return CF<BIGK>::isPalindrome(itr);
      if ( CF<BIGK-1>::isPalindrome(itr) )
          return true;
      return CF<BIGK-1>::isPalindrome(++itr); }

    bool upstreamExtensionPossible( BigKmer const& entry )
    { KMerContext context = entry.getContext();
      if ( context.getPredecessorCount() != 1 )
        return false;
      mEdgeSeq.clear().push_back(context.getSinglePredecessor());
      mEdgeSeq.append(entry.begin(),entry.end()-1);
      if ( CF<BIGK>::isPalindrome(mEdgeSeq.begin()) )
        return false;
      BigKmer kmer(mEdgeSeq,KMerHasher<BIGK>()(mEdgeSeq.begin()));
      lookup(kmer,&context);
      return context.getSuccessorCount() == 1; }

    bool downstreamExtensionPossible( BigKmer const& entry )
    { KMerContext context = entry.getContext();
      if ( context.getSuccessorCount() != 1 )
        return false;
      mEdgeSeq.assign(entry.begin()+1,entry.end());
      mEdgeSeq.push_back(context.getSingleSuccessor());
      if ( CF<BIGK>::isPalindrome(mEdgeSeq.begin()) )
        return false;
      BigKmer kmer(mEdgeSeq,KMerHasher<BIGK>()(mEdgeSeq.begin()));
      lookup(kmer,&context);
      return context.getPredecessorCount() == 1; }

    void make1KmerEdge( BigKmer const& entry )
    { mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(&entry);
      addEdge(); }

    void extend( KMerContext context )
    { auto itr = mEdgeSeq.begin();
      KMerHasher<BIGK> hasher;
      BigKmer kmer(mEdgeSeq,hasher.hash(itr));
      while ( context.getSuccessorCount() == 1 )
      { mEdgeSeq.push_back(context.getSingleSuccessor());
        kmer.successor(hasher.stepF(++itr));
        if ( isPalindrome(kmer) )
        { mEdgeSeq.pop_back();
          break; }
        BigKmer const* pEntry = lookup(kmer,&context);
        if ( context.getPredecessorCount() != 1 )
        { mEdgeSeq.pop_back();
          break; }
        mEdgeEntries.push_back(pEntry); }
      switch (mEdgeSeq.getCanonicalForm())
      {case CanonicalForm::PALINDROME:
         ForceAssertEq(mEdgeSeq.size(),BIGK);
        // allow flow-through to FWD case
       case CanonicalForm::FWD:
         addEdge(); break;
       case CanonicalForm::REV:
         mEdgeSeq.clear(); mEdgeEntries.clear(); break; } }

    BigKmer const* lookup( BigKmer const& kmer, KMerContext* pContext )
    { BigKmer const* result;
      if ( kmer.isRev() )
      { result = mDict.lookup(kmer.rc());
        ForceAssert(result);
        *pContext = result->getContext().rc(); }
      else
      { result = mDict.lookup(kmer);
        ForceAssert(result);
        *pContext = result->getContext(); }
      return result; }

    void addEdge()
    { static SpinLockedData gLock;
      size_t edgeId;
      if ( mEdgeSeq.getCanonicalForm() == CanonicalForm::REV )
          mEdgeSeq.ReverseComplement();
      if ( true )
      { SpinLocker lock(gLock);
        edgeId = mEdges.size();
        mEdges.push_back(mEdgeSeq); }
      unsigned offset = 0;
      bool err = false;
      for ( BigKmer const* pEnt : mEdgeEntries )
      { if ( pEnt->isUnassigned() )
          const_cast<BigKmer*>(pEnt)->setAssigned();
        else
        { std::cout << edgeId << ':' << offset++ << ' ' << pEnt
                    << " Already occupied." << std::endl;
          err = true; } }
        if ( err )
            FatalErr("Having trouble with preoccupied kmers.");
        mEdgeSeq.clear();
        mEdgeEntries.clear(); }

    BigKDict const& mDict;
    vecbvec& mEdges;
    std::vector<BigKmer const*> mEdgeEntries;
    bvec mEdgeSeq;
};

template <unsigned BIGK>
size_t N50( vecbvec const& reads )
{
    size_t nKmers = reads.getKmerCount(BIGK);
    BigDict<BIGK> bigDict(nKmers);
    BigKMerizer<BIGK> kmerizer(&bigDict);
    parallelFor(0ul,reads.size(),
        [kmerizer,&reads]( size_t readId ) mutable
        { kmerizer.kmerize(reads[readId]); });

    vecbvec edges;
    edges.reserve(bigDict.size()/100);
    BigKEdgeBuilder<BIGK>::buildEdges(bigDict,&edges);
    ForceAssertEq(edges.getKmerCount(BIGK),bigDict.size());

    vec<size_t> lengths;
    size_t nnn = edges.size();
    lengths.reserve(nnn);
    for ( bvec const& bv : edges )
        lengths.push_back(bv.size());
    std::sort(lengths.begin(),lengths.end());
    return N50(lengths);
}

}

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArgumentsAcceptEmptyArgListNoHeader;
    CommandArgument_String(READS);
    EndCommandArguments;

    vecbvec reads(READS);
    std::cout << "100\t" << N50<100ul>(reads) << std::endl;
    std::cout << "200\t" << N50<200ul>(reads) << std::endl;
    std::cout << "300\t" << N50<300ul>(reads) << std::endl;
    std::cout << "400\t" << N50<400ul>(reads) << std::endl;
    std::cout << "500\t" << N50<500ul>(reads) << std::endl;
    std::cout << "600\t" << N50<600ul>(reads) << std::endl;
    std::cout << "700\t" << N50<700ul>(reads) << std::endl;
    std::cout << "800\t" << N50<800ul>(reads) << std::endl;
    std::cout << "900\t" << N50<900ul>(reads) << std::endl;
    std::cout << "1000\t" << N50<1000ul>(reads) << std::endl;
}
