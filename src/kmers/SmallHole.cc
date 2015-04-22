///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SmallHole.cc
 * \author tsharpe
 * \date Apr 5, 2012
 *
 * \brief
 */

#include "MainTools.h"

#include "kmers/ReadPatherDefs.h"
#include "system/Assert.h"
#include "system/SpinLockedData.h"
#include "system/SysConf.h"
#include "system/System.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

namespace
{

unsigned const K = 25;

void dumpFASTA( ostream& os, HugeBVec::const_iterator itr,
                    HugeBVec::const_iterator end )
{
    std::ostream_iterator<char> outItr(os);
    while ( end-itr >= 80 )
    {
        std::transform(itr,itr+80,outItr,&Base::val2Char);
        os << '\n';
        itr += 80;
    }
    if ( itr != end )
    {
        std::transform(itr,end,outItr,&Base::val2Char);
        os << '\n';
    }
}

void dumpFASTG( ostream& os, UnipathGraph<K> const& ug )
{
    UnipathEdgeVec const& edges = ug.getAllEdges();
    os << "#begin FASTG version0.62\n";
    size_t connId = 0;
    size_t nnn = edges.size();
    for ( size_t idx = 0; idx != nnn; ++idx )
    {
        ComponentID const& cID = edges[idx].getComponentID();
        size_t end = idx + 1;
        while ( end != nnn && edges[end].getComponentID() == cID )
            end += 1;
        EdgeID edgeID(idx);
        if ( end-idx == 1 )
        {
            os << "\n>" << "edge" << idx << '\n';
            dumpFASTA(os,ug.getBases(edgeID),ug.getBasesEnd(edgeID));
        }
        else
        {
            os << '\n';
            for ( size_t idx2 = idx; idx2 != end; ++idx2 )
            {
                EdgeID edgeID2(idx2);
                size_t node = 2*idx2;
                os << '>' << "edge" << idx2 << ':'
                        << node << '-' << node+1 << ";\n";
                dumpFASTA(os,ug.getBases(edgeID2),ug.getBasesEnd(edgeID2));
            }
            for ( size_t idx2 = idx; idx2 != end; ++idx2 )
            {
                UnipathEdge const& edge = edges[idx2];
                for ( unsigned base = 0; base < 4u; ++base )
                {
                    if ( edge.hasPredecessor(base) )
                    {
                        size_t edgeId = edge.getPredecessor(base).val();
                        if ( edgeId >= idx2 )
                        {
                            os << ">c" << connId++ << ':';
                            if ( edge.isPredecessorRC(base) )
                                os << 2*idx2 << '-' << 2*edgeId;
                            else
                                os << 2*edgeId+1 << '-' << 2*idx2;
                            os << ":connection;\n";
                        }
                    }
                    if ( edge.hasSuccessor(base) )
                    {
                        size_t edgeId = edge.getSuccessor(base).val();
                        if ( edgeId >= idx2 )
                        {
                            os << ">c" << connId++ << ':';
                            if ( edge.isSuccessorRC(base) )
                                os << 2*idx2+1 << '-' << 2*edgeId+1;
                            else
                                os << 2*idx2+1 << '-' << 2*edgeId;
                            os << ":connection;\n";
                        }
                    }
                }
            }
        }
        idx = end - 1;
    }
    os << "\n#end FASTG\n";
}

class SeedPair
{
public:
    SeedPair( EdgeID const& edgeID1, EdgeID const& edgeID2,
                bool edge2Upstream, bool edge2NC,
                size_t readId, unsigned off1, unsigned off2 )
    : mEdgeID1(edgeID1), mEdgeID2(edgeID2),
      mEdge2Upstream(edge2Upstream), mEdge2NC(edge2NC),
      mReadID(readId), mOff1(off1), mOff2(off2) {}

    // compiler-supplied copying and destructor are OK

    EdgeID const& getEdgeID1() const { return mEdgeID1; }
    EdgeID const& getEdgeID2() const { return mEdgeID2; }
    bool isEdge2Upstream() const { return mEdge2Upstream; }
    bool isEdge2NonCanonical() const { return mEdge2NC; }
    size_t getGapSize() const
    { return mOff1 > mOff2 ? mOff1-mOff2 : mOff2-mOff1; }

    friend bool samePair( SeedPair const& sp1, SeedPair const& sp2 )
    { return sp1.mEdgeID1 == sp2.mEdgeID1 &&
                sp1.mEdgeID2 == sp2.mEdgeID2 &&
                sp1.mEdge2Upstream == sp2.mEdge2Upstream &&
                sp1.mEdge2NC == sp2.mEdge2NC; }

    friend int compare( SeedPair const& sp1, SeedPair const& sp2 )
    { int result = compare(sp1.mEdgeID1,sp2.mEdgeID1);
      if ( !result )
      { result = compare(sp1.mEdgeID2,sp2.mEdgeID2);
        if ( !result )
        { result = compare(sp1.mEdge2Upstream,sp2.mEdge2Upstream);
          if ( !result )
          { result = compare(sp1.mEdge2NC,sp2.mEdge2NC);
            if ( !result )
              result = compare(sp1.getGapSize(),sp2.getGapSize()); } } }
      return result; }

    friend bool operator<( SeedPair const& sp1, SeedPair const& sp2 )
    { return compare(sp1,sp2) < 0; }

    friend std::ostream& operator<<( std::ostream& os, SeedPair const& sp )
    { os << sp.mEdgeID1 << '\t' << sp.mEdgeID2 << '\t'
         << (sp.mEdge2Upstream?'U':'D') << (sp.mEdge2NC?'N':'C') << '\t'
         << sp.mReadID << '\t' << sp.getGapSize()
         << '\t' << sp.mOff1 << '\t' << sp.mOff2;
      return os; }

private:
    EdgeID mEdgeID1;
    EdgeID mEdgeID2;

    // Edge-2 may be either upstream (in the 5' direction) or downstream
    // of a canonical (meaning, as it happened to be recorded in the
    // HugeBVec) edge-1.
    bool mEdge2Upstream;

    // Edge-2 may be either canonical or non-canonical w.r.t. a canonical
    // edge-1.
    bool mEdge2NC;

    // The part of the read that spans the gap between edge-1 and edge-2.
    ID<6> mReadID;
    unsigned mOff1;
    unsigned mOff2;
};

typedef std::vector<SeedPair> VecSeedPair;
typedef VecSeedPair::const_iterator SPItr;

//TODO: make this statistically robust
class EvidenceEvaluator
{
public:
    // default constructor must return an object where isAdequate is false,
    // and the gap size is enormous.
    EvidenceEvaluator() : mN(1), mSum(~0ul), mSum2(0) {}

    EvidenceEvaluator( SPItr itr, SPItr const& end )
    : mItr(itr), mN(0), mSum(0), mSum2(0)
    { using std::distance;
      mN = distance(itr,end);
      while ( itr != end )
      { size_t gapSize = itr->getGapSize();
        mSum += gapSize;
        mSum2 += gapSize*gapSize;
        ++itr; } }

    // compiler-supplied copying and destructor are OK

    bool isNull() const { return mN == 1ul; }
    bool isUpstream() const { return mItr->isEdge2Upstream(); }
    EdgeID const& getDownstreamEdgeID() const
    { return isUpstream() ? mItr->getEdgeID1() : mItr->getEdgeID2(); }
    EdgeID const& getUpstreamEdgeID() const
    { return isUpstream() ? mItr->getEdgeID2() : mItr->getEdgeID1(); }

    bool isAdequate() const
    { if ( mN < 3 ) return false;
      double mean = 1.*mSum/mN;
      double sDev = (1.*mSum2 - 1.*mSum*mSum/mN)/mN;
      return sDev < sqrt(mean); }

    size_t getGapSize() const { return (mSum+mN/2)/mN; }

    friend std::ostream& operator<<( std::ostream& os,
                                        EvidenceEvaluator const& ee )
    { for ( SPItr itr(ee.mItr), end(ee.mItr+ee.mN); itr != end; ++itr )
        os << *itr << '\n';
      return os; }

private:
    SPItr mItr;
    size_t mN;
    size_t mSum;
    size_t mSum2;
};

template <unsigned K>
class SeedPairKmerEater
{
public:
    SeedPairKmerEater( UnipathGraph<K> const& ug, VecSeedPair* pairs )
    : mpUG(&ug), mpPairs(pairs), mOffset(0) {}

    // compiler-supplied copying and destructor are OK

    void operator()( KMer<K> const& kmer, KMerContext, size_t rId, size_t rOff )
    { KDef const* pKDef = mpUG->getDict().lookup(kmer);
      if ( pKDef )
      { EdgeID const& nextID = pKDef->getEdgeID();
        bool nextNC = mpUG->getBases(*pKDef)[K/2] != kmer[K/2];
        size_t edgeOffset = pKDef->getEdgeOffset();
        if ( nextID != mLastID )
        { if ( !mLastID.isNull() )
          { size_t offset = rOff -
                          (nextNC ?
                            mpUG->getEdge(nextID).getLength() - edgeOffset :
                            edgeOffset + 1);
            if ( mLastID < nextID )
               mpPairs->push_back(SeedPair(mLastID,nextID,
                                           mLastNC,mLastNC!=nextNC,
                                           rId,mOffset,offset+1));
            else
               mpPairs->push_back(SeedPair(nextID,mLastID,
                                           !nextNC,mLastNC!=nextNC,
                                           rId,offset,mOffset-1)); }
          mLastID = nextID; }
        mLastNC = nextNC;
        mOffset = rOff +
                    (nextNC ?
                      edgeOffset + 1 :
                      mpUG->getEdge(nextID).getLength() - edgeOffset); } }

private:
    UnipathGraph<K> const* mpUG;
    VecSeedPair* mpPairs;
    EdgeID mLastID;
    bool mLastNC;
    size_t mOffset;
};

template <unsigned K>
class SeedPairProcessor
{
public:
    SeedPairProcessor( VirtualMasterVec<bvec> const& reads,
                     UnipathGraph<K> const& ug,
                     size_t batchSize,
                     VecSeedPair* pOut )
    : mReads(reads), mpUG(&ug), mBatchSize(batchSize), mpOut(pOut)
    { mOut.reserve(1000000); }

    // compiler-supplied copying and destructor are OK

    void operator()( size_t batchNo );

private:
    VirtualMasterVec<bvec> mReads;
    UnipathGraph<K> const* mpUG;
    size_t mBatchSize;
    VecSeedPair* mpOut;
    VecSeedPair mOut;
    static SpinLockedData gLock;
};

template <unsigned K>
void SeedPairProcessor<K>::operator()( size_t batchNo )
{
    size_t idx = batchNo * mBatchSize;
    if ( idx > mReads.size() ) idx = mReads.size();
    size_t stop = idx + mBatchSize;
    if ( stop > mReads.size() ) stop = mReads.size();
    typedef VirtualMasterVec<bvec>::const_iterator Itr;
    Itr end(mReads.begin(stop));
    for ( Itr itr(mReads.begin(idx)); itr != end; ++itr )
    {
        SeedPairKmerEater<K> eater(*mpUG,&mOut);
        KMer<K>::kmerizeIntoEater(itr->begin(),itr->end(),eater,idx++);
    }
    if ( mOut.size() )
    {
        SpinLocker locker(gLock);
        mpOut->insert(mpOut->end(),mOut.begin(),mOut.end());
    }
    mOut.clear();
}

template <unsigned K>
SpinLockedData SeedPairProcessor<K>::gLock;

void bestEvidence( VecSeedPair const& pairs )
{
    SPItr beg(pairs.begin());
    SPItr end(pairs.end());
    EdgeID curEdge = beg->getEdgeID1();
    EvidenceEvaluator bestEE[2];
    while ( beg != end )
    {
        if ( beg->getEdgeID1() != curEdge )
        {
            if ( bestEE[1].isAdequate() )
            {
                std::cout << bestEE[1];
                bestEE[1] = EvidenceEvaluator();
            }
            if ( bestEE[0].isAdequate() )
            {
                std::cout << bestEE[0];
                bestEE[0] = EvidenceEvaluator();
            }
            curEdge = beg->getEdgeID1();
        }
        SPItr itr(beg);
        while ( ++itr != end )
            if ( !samePair(*itr,*beg) )
                break;
        EvidenceEvaluator ee(beg,itr);
        if ( ee.isAdequate() )
        {
            EvidenceEvaluator& best = bestEE[ee.isUpstream()];
            if ( ee.getGapSize() < best.getGapSize() )
                best = ee;
        }
        beg = itr;
    }
    if ( bestEE[1].isAdequate() )
        std::cout << bestEE[1];
    if ( bestEE[0].isAdequate() )
        std::cout << bestEE[0];
}

void siftEvidence( VecSeedPair const& pairs,
                        bool* pHasPred, EvidenceEvaluator* pSucc )
{
    SPItr beg(pairs.begin());
    SPItr end(pairs.end());
    while ( beg != end )
    {
        SPItr itr(beg);
        while ( ++itr != end )
            if ( !samePair(*itr,*beg) )
                break;
        EvidenceEvaluator ee(beg,itr);
        if ( ee.isAdequate() )
        {
            EvidenceEvaluator* pEE = pSucc + ee.getUpstreamEdgeID().val();
            if ( ee.getGapSize() < pEE->getGapSize() )
            {
                *pEE = ee;
                pHasPred[ee.getDownstreamEdgeID().val()] = true;
            }
        }
        beg = itr;
    }
}

template <unsigned K>
void dumpContigs( UnipathGraph<K> const& ug,
                        bool const* pHasPred,
                        EvidenceEvaluator* pSucc )
{
    size_t contigID = 0;
    size_t nnn = ug.getNEdges();
    for ( size_t idx = 0; idx < nnn; ++idx )
    {
        if ( !pHasPred[idx] )
        {
            std::cout << ">contig" << contigID++ << '\n';
            dumpFASTA( std::cout, ug.getBases(EdgeID(idx)),
                        ug.getBasesEnd(EdgeID(idx)) );
            EvidenceEvaluator* pEE = pSucc + idx;
            while ( !pEE->isNull() )
            {
                size_t nnn = pEE->getGapSize();
                while ( nnn-- )
                    std::cout << 'N';
                std::cout << '\n';
                EdgeID const& edgeID = pEE->getDownstreamEdgeID();
                dumpFASTA( std::cout, ug.getBases(edgeID),
                            ug.getBasesEnd(edgeID) );
                *pEE = EvidenceEvaluator();
                pEE = pSucc + edgeID.val();
            }
        }
    }
}

} // end of anonymous namespace

int main( int argc, char** argv )
{
    String const EMPTY;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(READS, "Reads.");
    CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE,30u,"Coverage.");
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0u,"Threads.");
    CommandArgument_String_OrDefault_Doc(GENOME, EMPTY, "Genome.");
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);

    size_t const BATCH_SIZE = 100;

    VirtualMasterVec<bvec> reads(READS);
    size_t dictSize = 4*reads.sizeSum();
    std::cout << Date() << ": Creating dictionary." << std::endl;
    KmerDict<K> dict( 5*dictSize/4 );
    dict.process(reads,NUM_THREADS,BATCH_SIZE);
    std::cout << Date( ) << ": there are " << dict.size() << " kmers (expected ~"
              << dictSize << ")" << std::endl;

    size_t const MAX_COUNT = 200;
    Spectrum<K> spectrum(dict,MAX_COUNT);
    std::cout << "Kmer spectrum:" << std::endl;
    std::cout << spectrum << std::endl;
    size_t troughIdx = spectrum.getFirstTroughIndex();
    size_t peakIdx = spectrum.getFirstPeakIndex();
    std::cout << "Trough has " << spectrum[troughIdx] << " kmers at "
              << troughIdx << "x\nPeak has " << spectrum[peakIdx]
              << " kmers at " << peakIdx << 'x' << std::endl;

    if ( !spectrum.isFirstPeakValid() )
        FatalErr("Unable to find CN=1 peak.");
    if ( !spectrum.isPeakDistinct() )
        FatalErr("Insufficient peak/trough separation.");

    double delta = sqrt(peakIdx);
    size_t minIdx = peakIdx - .5*delta;
    size_t maxIdx = peakIdx + .5*delta + 1;
    std::cout << Date() << ": Extracting CN=1 (" << minIdx << ':' << maxIdx
                << ") kmers." << std::endl;

    dictSize = spectrum.sum(minIdx,maxIdx);
    dict.clean(KmerDict<K>::BadKmerCountFunctor(minIdx,maxIdx));
    AssertEq(dictSize,dict.size());
    std::cout << Date( ) << ": there are " << dictSize << " CN=1 kmers" << std::endl;

    if ( GENOME != EMPTY )
    {
        std::cout << Date() << ": Checking CN=1 kmers against genome."
                << std::endl;
        VirtualMasterVec<bvec> genome(GENOME);
        dictSize = 4*genome.sizeSum();
        KmerDict<K> genomeDict( 5*dictSize/4 );
        genomeDict.process(genome,NUM_THREADS,1);
        std::cout << "There are " << genomeDict.size()
                  << " genomic kmers (expected ~" << dictSize << ")."
                  << std::endl;

        size_t const GENOME_MAX_COUNT = 10;
        std::vector<size_t> kCounts(GENOME_MAX_COUNT+1);
        typedef KmerDict<K>::OCItr OCItr;
        typedef KmerDict<K>::ICItr ICItr;
        for ( OCItr oItr(dict.cbegin()),oEnd(dict.cend()); oItr!=oEnd; ++oItr )
        {
            for ( ICItr itr(oItr->begin()), end(oItr->end()); itr!=end; ++itr )
            {
                size_t count = 0;
                KDef const* pKDef = genomeDict.lookup(*itr);
                if ( pKDef )
                {
                    count = pKDef->getCount();
                    if ( count > GENOME_MAX_COUNT )
                        count = GENOME_MAX_COUNT;
                }
                kCounts[count] += 1;
            }
        }
        std::cout << "True copy numbers of the ostensibly CN=1 kmers:"
                  << std::endl;
        for ( size_t idx = 0; idx <= GENOME_MAX_COUNT; ++idx )
            std::cout << idx << '\t' << kCounts[idx] << '\n';
    }

    UnipathGraph<K> ug(dict);

    std::cout << Date() << ": Dumping FASTG." << std::endl;
    ofstream osFASTG(READS.ReplaceExtension(".fastb",".fastg").c_str());
    dumpFASTG(osFASTG,ug);
    osFASTG.close();

    VecSeedPair pairs;
    pairs.reserve(20000000ul);
    if ( true )
    {
        SeedPairProcessor<K> sp(reads,ug,BATCH_SIZE,&pairs);
        size_t nnn = (reads.size()+BATCH_SIZE-1)/BATCH_SIZE;
        std::cout << Date() << ": Finding seed pairs in " << nnn
                << " batches of " << BATCH_SIZE << '.' << std::endl;
        parallelFor(0ul,nnn,sp,NUM_THREADS);
    }

    std::cout << Date() << ": Sorting seed pairs." << std::endl;
    sortInPlaceParallel(pairs.begin(),pairs.end(),NUM_THREADS);

    std::cout << Date() << ": Sifting evidence." << std::endl;
#if 1
    bestEvidence(pairs);
#else
    bool* pHasPred = new bool[ug.getNEdges()];
    EvidenceEvaluator* pSucc = new EvidenceEvaluator[ug.getNEdges()];
    siftEvidence(pairs,pHasPred,pSucc);

    std::cout << Date() << ": Writing contigs." << std::endl;
    dumpContigs(ug,pHasPred,pSucc);
    delete [] pSucc;
    delete [] pHasPred;
#endif
    std::cout << Date() << ": Done." << std::endl;
}
