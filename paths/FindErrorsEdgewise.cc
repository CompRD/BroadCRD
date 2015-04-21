///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FindErrorsNew.cc
 * \author tsharpe
 * \date Jan 11, 2011
 *
 * \brief
 */
#include "MainTools.h"
#include "Basevector.h"
#include "Charvector.h"
#include "FeudalMimic.h"
#include "PairsManager.h"
#include "dna/Bases.h"
#include "feudal/BinaryStreamTraits.h"
#include "feudal/BitVec.h"
#include "feudal/QualNibbleVec.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/ReadPatherDefs.h"
#include "system/ExternalSorter.h"
#include "system/LockedData.h"
#include "system/SortInPlace.h"
#include "system/SysConf.h"
#include "system/System.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <cstring>
#include <vector>

#define ULOG(x) std::cout << x << std::endl
#define DLOG(x) ULOG(Date() << " (FE): " << x)

namespace
{

// worklist procedure for sorting kmers
template <unsigned K>
class KmerSortProc
{
public:
    typedef std::vector<KMer<K> > VecKmer;
    class Common
    {
    public:
        Common( KmerDict<K>& dict, vecbvec const& reads, unsigned minKmerFreq,
                unsigned sortNThreads, size_t nBatches, size_t batchSize,
                size_t kmersPerBatch, String const& tmpDir )
        : mSorter(tmpDir+"/FEKmers.tmp."), mDict(dict), mReads(reads),
          mMinKmerFreq(minKmerFreq), mSortNThreads(sortNThreads),
          mDotter(nBatches), mBatchSize(batchSize),
          mKmersPerBatch(kmersPerBatch)
        {}

        size_t getBatchSize() const { return mBatchSize; }
        size_t getKmersPerBatch() const { return mKmersPerBatch; }
        ExternalSorter<KMer<K> >& getSorter() { return mSorter; }
        size_t getNReads() const { return mReads.size(); }
        bvec const& getRead( size_t idx ) const { return mReads[idx]; }

        void dump( VecKmer& vk )
        { typedef typename VecKmer::iterator Itr;
          using std::swap;
          Itr begin(vk.begin()); Itr end(vk.end());
          if ( mSortNThreads > 1 )
              sortInPlaceParallel(begin,end,mSortNThreads);
          else
              sortInPlace(begin,end);
          Itr out(begin);
          while ( begin != end )
          { KMer<K> const& kmer = *begin;
            Itr itr(begin);
            while ( ++itr != end )
              if ( *itr != kmer ) break;
            if ( itr - begin >= mMinKmerFreq ) mDict.insertCanonical(kmer);
            else if ( begin != out )
            { while ( begin != itr ) { swap(*begin,*out); ++begin; ++out; } }
            begin = itr; }
          begin = vk.begin();
          if ( begin != out )
              mSorter.insert(begin,out);
          vk.clear();
          mDotter.batchDone(); }

    private:
        Common( Common const& ); // unimplemented -- no copying
        Common& operator=( Common const& ); // unimplemented -- no copying

        ExternalSorter<KMer<K> > mSorter;
        KmerDict<K>& mDict;
        vecbvec const& mReads;
        unsigned mMinKmerFreq;
        unsigned mSortNThreads;
        Dotter mDotter;
        size_t mBatchSize;
        size_t mKmersPerBatch;
    };

    KmerSortProc( Common& common )
    : mCommon(common)
    {}

    KmerSortProc( KmerSortProc const& that )
    : mCommon(that.mCommon)
    { mKmers.reserve(mCommon.getKmersPerBatch()); }

    // compiler-supplied destructor OK, operator= cannot be implemented

    void operator()( size_t batchNo );

private:
    Common& mCommon;
    VecKmer mKmers;
};

template <unsigned K>
void KmerSortProc<K>::operator()( size_t batchNo )
{
    using std::min;
    size_t nReads = mCommon.getNReads();
    size_t batchSize = mCommon.getBatchSize();
    size_t readId = min(nReads,batchNo*batchSize);
    size_t endId = min(nReads,readId+batchSize);
    std::back_insert_iterator<VecKmer> oItr = std::back_inserter(mKmers);
    for ( ; readId != endId; ++readId )
    {
        bvec const& bv = mCommon.getRead(readId);
        KMer<K>::kmerize(bv.begin(),bv.end(),oItr);
    }

    mCommon.dump(mKmers);
}

// create a KmerDict on the heap from a fastb file of reads.
// only put into the dictionary those that appear minKmerFreq times or more.
template <unsigned K>
KmerDict<K>* createDict( vecbvec const& reads, unsigned nThreads,
                        unsigned minKmerFreq, unsigned estimatedCoverage,
                        String const& tmpDir )
{
    size_t nReads = reads.size();
    size_t avgReadLen = (reads.SizeSum() + nReads - 1) / nReads;
    size_t kmersPerRead = avgReadLen > K ? avgReadLen-K+1 : 1;
    size_t nKmers = kmersPerRead*nReads / estimatedCoverage;

    ULOG("Reserving kmer dictionary space for " << nKmers << " kmers.");
    KmerDict<K>* pDict = new KmerDict<K>(nKmers);

    size_t bytesPerRead = kmersPerRead*sizeof(KMer<K>);
    WorklistParameterizer wp(nReads,bytesPerRead,
                             nReads*minKmerFreq/estimatedCoverage,nThreads);
    size_t batchSize = wp.getBatchSize();
    size_t nBatches = wp.getNBatches();
    unsigned sortNThreads = nThreads/wp.getNThreads();
    nThreads = wp.getNThreads();

    typedef typename KmerSortProc<K>::Common Common;
    Common common(*pDict, reads, minKmerFreq, sortNThreads, nBatches, batchSize,
                    batchSize*kmersPerRead, tmpDir);
    ULOG("Kmerizing " << nReads << " reads in " << nBatches << " batches of "
            << batchSize << '.');

    parallelFor(0ul,nBatches,KmerSortProc<K>(common),nThreads);

    if ( nThreads > 1 )
    {
        KMer<K> lastKmer;
        KMer<K> thisKmer;
        size_t kCount = 0;
        typedef typename ExternalSorter<KMer<K> >::template Merger<> Merger;
        Merger merger(common.getSorter());
        while ( merger.getNext(thisKmer) )
        {
            if ( thisKmer == lastKmer )
                kCount += 1;
            else
            {
                if ( kCount >= minKmerFreq )
                    pDict->insertCanonical(lastKmer);
                lastKmer = thisKmer;
                kCount = 1;
            }
        }
    }

    ULOG("\nThere are " << pDict->size() << " kmers in the dictionary.");
    return pDict;
}

// An identified read's alignment onto an unidentified edge.
// The offset is described in a tricky way to make it a non-negative value:
// It describes where to place the last base of the read relative to the first
// base of the unipath.  Thus, the minimum offset of K indicates that the final
// kmer of the read matches the first kmer of the unipath.  The maximum offset
// of edgeLen+K-1 shows that the first kmer in the read matches final kmer of
// the edge. For RC alignments, the offset in interpreted the same way, but on
// the RC'd read.
class ReadLoc
{
public:
    ReadLoc() {}

    ReadLoc( ReadID const& readID, size_t offset, bool isRC )
    : mOffset(2*offset+isRC), mReadID(readID)
    { ForceAssertLt(offset,1ul<<31); }

    ReadLoc( ReadLoc const& that )
    { memcpy(this,&that,sizeof(ReadLoc)); }

    // compiler-supplied destructor is OK

    ReadLoc& operator=( ReadLoc const& that )
    { memcpy(this,&that,sizeof(ReadLoc)); return *this; }

    ReadID const& getReadID() const { return mReadID; }
    size_t getOffset() const { return mOffset >> 1; }
    bool isRC() const { return mOffset&1; }

    friend int compare( ReadLoc const& rl1, ReadLoc const& rl2 )
    { int result = compare(rl1.mOffset,rl2.mOffset);
      if ( !result ) result = compare(rl1.mReadID,rl2.mReadID);
      return result; }

private:
    unsigned mOffset;
    ReadID mReadID;
};

// Describes a unipath (edge) and a read that share at least one kmer.
class ReadLineup
{
public:
    ReadLineup() {}
    ReadLineup( EdgeID const& edgeID, ReadID const& readID, size_t offset,
                bool isRC )
    : mEdgeID(edgeID), mReadLoc(readID,offset,isRC) {}
    ReadLineup( ReadLineup const& that )
    : mEdgeID(that.mEdgeID), mReadLoc(that.mReadLoc) {}

    // compiler-supplied destructor is OK

    ReadLineup& operator=( ReadLineup const& that )
    { memcpy(this,&that,sizeof(ReadLineup)); return *this; }

    EdgeID const& getEdgeID() const { return mEdgeID; }
    ReadID const& getReadID() const { return mReadLoc.getReadID(); }
    size_t getOffset() const { return mReadLoc.getOffset(); }
    bool isRC() const { return mReadLoc.isRC(); }
    ReadLoc const& getReadLoc() const { return mReadLoc; }

    friend void swap( ReadLineup& rl1, ReadLineup& rl2 )
    { char buf[sizeof(ReadLineup)];
      memcpy(buf,&rl1,sizeof(ReadLineup));
      memcpy(&rl1,&rl2,sizeof(ReadLineup));
      memcpy(&rl2,buf,sizeof(ReadLineup)); }

    friend int compare( ReadLineup const& rl1, ReadLineup const& rl2 )
    { int result = compare(rl1.mEdgeID,rl2.mEdgeID);
      if ( !result ) result = compare(rl1.mReadLoc,rl2.mReadLoc);
      return result; }

    friend bool operator<( ReadLineup const& rl1, ReadLineup const& rl2 )
    { return compare(rl1,rl2) < 0; }
    friend bool operator<=( ReadLineup const& rl1, ReadLineup const& rl2 )
    { return compare(rl1,rl2) <= 0; }
    friend bool operator>( ReadLineup const& rl1, ReadLineup const& rl2 )
    { return compare(rl1,rl2) > 0; }
    friend bool operator>=( ReadLineup const& rl1, ReadLineup const& rl2 )
    { return compare(rl1,rl2) >= 0; }
    friend bool operator==( ReadLineup const& rl1, ReadLineup const& rl2 )
    { return !memcmp(&rl1,&rl2,sizeof(ReadLineup)); }
    friend bool operator!=( ReadLineup const& rl1, ReadLineup const& rl2 )
    { return memcmp(&rl1,&rl2,sizeof(ReadLineup)); }

    friend ostream& operator<<( ostream& os, ReadLineup const& rl )
    { os << rl.getEdgeID() << ':' << rl.getOffset() <<
        (rl.isRC() ? '-' : '+') << rl.getReadID();
      return os; }

private:
    EdgeID mEdgeID;
    ReadLoc mReadLoc;
};


// a worklist procedure to make and externally sort ReadLineups
template <unsigned K>
class LineupProc
{
public:
    typedef std::vector<ReadLineup> VecReadLineup;

    class Common
    {
    public:
        Common( UnipathGraph<K> const& graph, vecbvec const& reads,
                size_t nBatches, size_t batchSize, bool doPalindromes,
                ExternalSorter<ReadLineup>* pXS )
        : mGraph(graph), mReads(reads), mBatchSize(batchSize),
          mDoPalindromes(doPalindromes), mXS(*pXS), mDotter(nBatches)
        {}

        size_t getBatchSize() const { return mBatchSize; }
        size_t getNReads() const { return mReads.size(); }
        bvec const& getRead( size_t idx ) const { return mReads[idx]; }

        KDef const* lookup( KMer<K> const& kmer ) const
        { return mGraph.getDict().lookup(kmer); }

        bool doPalindromes() const { return mDoPalindromes; }
        bool isPalindrome( EdgeID const& edgeID ) const
        { return mGraph.getEdge(edgeID).isPalindrome(); }

        size_t getEdgeLen( EdgeID const& edgeID ) const
        { return mGraph.getEdge(edgeID).getLength(); }

        HugeBVec::const_iterator getBases( KDef const& kDef ) const
        { return mGraph.getBases(kDef); }

        void dump( VecReadLineup& vrl )
        { typedef VecReadLineup::iterator Itr;
          Itr start(vrl.begin()); Itr stop(vrl.end());
          sortInPlace(start,stop);
          using std::unique; mXS.insert(start,unique(start,stop));
          vrl.clear();
          mDotter.batchDone(); }

    private:
        Common( Common const& ); // unimplemented -- no copying
        Common& operator=( Common const& ); // unimplemented -- no copying

        UnipathGraph<K> const& mGraph;
        vecbvec const& mReads;
        size_t mBatchSize;
        bool mDoPalindromes;
        ExternalSorter<ReadLineup>& mXS;
        Dotter mDotter;
    };

    LineupProc( Common& common ) : mCommon(common) {}

    LineupProc( LineupProc const& that ) : mCommon(that.mCommon)
    { mVRL.reserve(EDGES_PER_READ*mCommon.getBatchSize()); }

    // copying prohibited by reference, compiler-supplied destructor is OK

    void operator()( size_t batchNo );

    // assume that there are 3 unipaths per read, on average
    static size_t const EDGES_PER_READ = 3;

private:
    typedef bvec::const_iterator RItr;
    typedef HugeBVec::const_iterator UItr;

    void alignRead( ReadID const& readID, bvec const& bv );

    size_t getFwdMatchCount( RItr rdItr, UItr edgeItr, size_t len )
    { size_t matchCount = 0;
      while ( matchCount < len && *rdItr == *edgeItr )
      { ++matchCount; ++rdItr; ++edgeItr; }
      return matchCount; }

    size_t getRevMatchCount( RItr rdItr, UItr edgeItr, size_t len )
    { size_t matchCount = 0;
      while ( matchCount < len && *rdItr == GetComplementaryBase(*edgeItr) )
      { ++matchCount; ++rdItr; --edgeItr; }
      return matchCount; }

    Common& mCommon;
    VecReadLineup mVRL;
};

template <unsigned K>
void LineupProc<K>::operator()( size_t batchNo )
{
    using std::min;
    size_t nReads = mCommon.getNReads();
    size_t batchSize = mCommon.getBatchSize();
    size_t readId = min(nReads,batchNo*batchSize);
    size_t endId = min(nReads,readId+batchSize);
    for ( ; readId != endId; ++readId )
        alignRead(ReadID(readId),mCommon.getRead(readId));

    mCommon.dump(mVRL);
}

template <unsigned K>
void LineupProc<K>::alignRead( ReadID const& readID, bvec const& bv )
{
    size_t rdLen = bv.size();
    if ( rdLen < K )
        return; // EARLY RETURN!

    RItr begin(bv.begin());
    RItr eEnd(bv.end() - K + 1);
    while ( true )
    {
        RItr itr(begin);
        KDef const* pDef = 0;
        for ( ; itr < eEnd; ++itr )
        {
            KMer<K> kkk(itr);
            if ( (pDef = mCommon.lookup(kkk)) )
            {
                if ( mCommon.doPalindromes()
                        || !mCommon.isPalindrome(pDef->getEdgeID()) )
                    break;
                pDef = 0;
            }
        }

        if ( !pDef )
            break;

        size_t rdLenRemaining = rdLen - itr.pos();
        size_t edgeOff = pDef->getEdgeOffset();
        UItr edgeItr = mCommon.getBases(*pDef);
        size_t edgeLen = mCommon.getEdgeLen(pDef->getEdgeID()) + K - 1;
        using std::min;
        size_t maxAlgnLen = min(rdLenRemaining, edgeLen - edgeOff);
        size_t matchCount = getFwdMatchCount(itr, edgeItr, maxAlgnLen);
        if ( matchCount >= K )
        {
            size_t offset = edgeOff + rdLenRemaining;
            ReadLineup rl(pDef->getEdgeID(), readID, offset, false);
            if ( mVRL.empty() || rl != mVRL.back() )
                mVRL.push_back(rl);
        }
        else
        {
            maxAlgnLen = min(rdLenRemaining, edgeOff + K);
            matchCount = getRevMatchCount(itr, edgeItr + K - 1, maxAlgnLen);
            if ( matchCount < K )
                FatalErr("Retrieved kmer matches neither fwd nor rev.");
            size_t offset = edgeOff + K + itr.pos();
            ReadLineup rl(pDef->getEdgeID(), readID, offset, true);
            if ( mVRL.empty() || rl != mVRL.back() )
                mVRL.push_back(rl);
        }

        begin = itr + matchCount - K + 1;
        // if we ended mid-edge, we know the next kmer isn't in the dictionary
        if ( matchCount < edgeLen )
            begin += 1;
    }
}

template <unsigned K>
void createReadLineups( UnipathGraph<K> const& graph, vecbvec const& reads,
                        unsigned nThreads, bool doPalindromes,
                        ExternalSorter<ReadLineup>* pXS )
{
    DLOG("Creating read to unipath alignments.");

    size_t nReads = reads.size();
    size_t bytesPerRead = LineupProc<K>::EDGES_PER_READ*sizeof(ReadLineup);
    WorklistParameterizer wp(nReads,bytesPerRead,1000,nThreads);
    size_t nBatches = wp.getNBatches();
    size_t batchSize = wp.getBatchSize();
    nThreads = wp.getNThreads();

    typedef typename LineupProc<K>::Common Common;
    Common common(graph, reads, nBatches, batchSize, doPalindromes, pXS);
    ULOG("Aligning " << nReads << " reads in " << nBatches << " batches of "
                << batchSize << '.');
    LineupProc<K> proc(common);
    parallelFor(0ul,nBatches,proc,nThreads);
}

template <unsigned K>
class StackProc
{
public:
    typedef std::vector<ReadLoc>::const_iterator LocItr;

    class Params
    {
    public:
        // QSS means quality score sum.  It's the sum of the quality scores
        // for a particular base call within a stack column.

        // Parameters that determine whether there is support for a call.
        // If a call is supported, it will not be corrected in the current
        // column.
        unsigned mSupportingScore; // a quality score deemed high-quality
        unsigned mSupportingCount; // no. of high-q calls needed to support
        unsigned mSupportingSum; // high-q sum needed to support (0 means ignore this criterion)

        // Parameters that determine whether there is a winning call in the
        // current column.
        unsigned mWinningSum; // the winner must have a QSS at least this big
        double mPlaceToWinRatio; // (2nd-largest-QSS/largest-QSS) must not exceed this

        // Parameters that determine whether there are loser (i.e., correctable)
        // calls in the current column.  In order to be corrected, there must be
        // a (different) winning call, the call must not be supported, and the
        // QSS must be sufficiently small relative to the winner.
        double mLoserToWinRatio; // loser-QSS/largest-QSS must not exceed this

        // Parameters that determine whether a call is confirmed.  Calls
        // confirmed by a column will not be corrected, even if they are losers
        // in that, or in some other, column.  There are two ways to be
        // confirmed: 1) the QSS must exceed an auto-confirm parameter, or
        unsigned mAutoConfirm; // auto-confirm QSS (0 means ignore this criterion)
        // 2) a call must be the winner, must be supported if current column is
        // within the edge boundaries, no other call can be uncorrectable, and...
        double mPlaceToWinRatioToConfirm; // (2nd-largest-QSS/largest-QSS) must not exceed this
        unsigned mConfirmingSum; // the QSS must be at least this big (0 means ignore this criterion)
        unsigned mConfirmingCount; // no. of high-q calls needed to confirm

        unsigned mReadStackDepth; // minimum no. of reads on an edge to bother with
        unsigned mBaseStackDepth; // minimum no. of calls in a column to bother with

        unsigned mQualCeilRadius; // use the minimum qual score within this radius for QSSs

        bool mDoBranches; // correct beyond edge boundaries even when branch is supported

        Params( unsigned supportingScore, unsigned supportingCount,
                double loserToWinRatio, double placeToWinRatio,
                double placeToWinRatioToConfirm, unsigned confirmingSum,
                unsigned winningSum, unsigned readStackDepth,
                unsigned baseStackDepth, unsigned autoConfirm,
                unsigned supportingSum, unsigned confirmingCount,
                unsigned qualCeilRadius, bool doBranches )
        : mSupportingScore(supportingScore),
          mSupportingCount(supportingCount),
          mSupportingSum(supportingSum),
          mWinningSum(winningSum),
          mPlaceToWinRatio(placeToWinRatio),
          mLoserToWinRatio(loserToWinRatio),
          mAutoConfirm(autoConfirm),
          mPlaceToWinRatioToConfirm(placeToWinRatioToConfirm),
          mConfirmingSum(confirmingSum),
          mConfirmingCount(confirmingCount),
          mReadStackDepth(readStackDepth),
          mBaseStackDepth(baseStackDepth),
          mQualCeilRadius(qualCeilRadius),
          mDoBranches(doBranches)
        {}
    };

    class WorkItem
    {
    public:
        WorkItem( EdgeID const& edgeID, size_t estLocs )
        : mEdgeID(edgeID)
        { mLocs.reserve(estLocs); }

        // compiler-supplied copying and destructor are OK

        EdgeID const& getEdgeID() const { return mEdgeID; }
        std::vector<ReadLoc> const& getLocs() const { return mLocs; }
        size_t getNLocs() const { return mLocs.size(); }

        void addLoc( ReadLoc const& rl ) { mLocs.push_back(rl); }

        friend void swap( WorkItem& item1, WorkItem& item2 )
        { using std::swap;
          swap(item1.mEdgeID,item2.mEdgeID);
          swap(item1.mLocs,item2.mLocs); }

    private:
        EdgeID mEdgeID;
        std::vector<ReadLoc> mLocs;
    };

    enum Action { NO_ACTION, CORRECT, CONFIRM, CONFLICT };

    class Scorer
    {
    public:
        Scorer( Params const& params ) : mParams(params)
        { clearScores(); }

        // compiler-supplied copying and destructor are OK

        void addScore( unsigned char base, unsigned qual )
        { BaseScore& bs = mBaseScores[base];
          bs.mQualSum += qual; bs.mQualCount += 1;
          if ( qual >= mParams.mSupportingScore ) bs.mSupportingCount += 1; }

        class Results
        {
        public:

            Results( bool isBranch )
            : mIsBranch(isBranch), mWinner(NO_WINNER), mAction(0) {}

            // compiler-supplied copying and destructor are OK

            bool isBranch() const { return mIsBranch; }

            Action getAction( unsigned char base ) const
            { return static_cast<Action>((mAction >> (2*base)) & 3); }

            void setAction( unsigned char base, Action action )
            { mAction ^= (getAction(base) ^ action) << (2*base); }

            bool isActionNecessary() const { return mAction; }

            unsigned char getWinner() const { return mWinner; }

            static unsigned char const NO_WINNER =
                                                static_cast<unsigned char>(-1);

            friend ostream& operator<<( ostream& os, Results const& results )
            { static char const* gActionCode = "0-!";
              os << gActionCode[results.getAction(0)]
                 << gActionCode[results.getAction(1)]
                 << gActionCode[results.getAction(2)]
                 << gActionCode[results.getAction(3)]
                 << (results.getWinner() == NO_WINNER ? 'X' :
                        Base::val2Char(results.getWinner()));
              return os; }

        private:
            bool mIsBranch;
            unsigned char mWinner;
            unsigned char mAction;
            friend class Scorer;
        };

        Results const getResults( bool onEdge );

        friend ostream& operator<<( ostream& os, Scorer const& scorer )
        { os << "A(" << scorer.mBaseScores[0]
             << ") C(" << scorer.mBaseScores[1]
             << ") G(" << scorer.mBaseScores[2]
             << ") T(" << scorer.mBaseScores[3] << ')';
          return os; }

    private:
        struct BaseScore
        {
            BaseScore() : mQualSum(0), mQualCount(0), mSupportingCount(0) {}

            unsigned mQualSum;
            unsigned mQualCount;
            unsigned mSupportingCount;

            friend ostream& operator<<( ostream& os, BaseScore const& bs )
            { os << bs.mQualSum << ':'
                << bs.mQualCount << ':'
                << bs.mSupportingCount;
              return os; }
        };

        bool isSupported( BaseScore const& bs ) const
        { return mParams.mSupportingSum ?
                    bs.mQualSum >= mParams.mSupportingSum :
                    bs.mSupportingCount >= mParams.mSupportingCount; }

        bool isCorrectable( unsigned char base, unsigned winningQSS ) const
        { BaseScore const& bs = mBaseScores[base];
          double loserRatio = 1.*bs.mQualSum/winningQSS;
          return loserRatio <= mParams.mLoserToWinRatio && !isSupported(bs); }

        void clearScores() { memset(mBaseScores,0,sizeof(mBaseScores)); }

        Params const& mParams;
        BaseScore mBaseScores[4];
    };

    class Common
    {
    public:
        Common( UnipathGraph<K> const& graph, Params const& params,
                vecbvec& reads, VecQualNibbleVec& quals,
                std::vector<WorkItem> const& workItems )
        : mGraph(graph), mReads(reads), mQuals(quals), mParams(params),
          mWorkItems(workItems), mDotter(workItems.size())
        { Mimic(reads,mActions); }

        // compiler-supplied destructor is OK

        void dot() { mDotter.batchDone(); }

        UnipathGraph<K> const& getGraph() const { return mGraph; }
        bvec const& getRead( size_t idx ) const { return mReads[idx]; }
        QualNibbleVec const& getQual( size_t idx ) const { return mQuals[idx]; }
        Params const& getParams() const { return mParams; }
        WorkItem const& getWorkItem( size_t idx ) { return mWorkItems[idx]; }

        void confirm( size_t readId, unsigned rdOffset )
        { mActions[readId][rdOffset] = CONFIRM; }

        void correct( size_t readId, unsigned rdOffset, unsigned char winner )
        { char volatile* pOldVal = &mActions[readId][rdOffset];
          while ( true )
          {
              char oldVal = *pOldVal;
              Action oldAction = static_cast<Action>(oldVal&3);
              if ( oldAction == CONFIRM || oldAction == CONFLICT ) break;
              char newVal = CORRECT | (winner << 2);
              if ( oldAction == CORRECT )
              {
                  if ( oldVal == newVal ) break;
                  newVal = CONFLICT;
              }
              if ( __sync_bool_compare_and_swap(pOldVal,oldVal,newVal) )
                  break; } }

        void applyActions();

    private:
        Common( Common const& ); // unimplemented -- no copying
        Common& operator=( Common const& ); // unimplemented -- no copying

        UnipathGraph<K> const& mGraph;
        vecbvec& mReads;
        VecQualNibbleVec& mQuals;
        Params mParams;
        std::vector<WorkItem> const& mWorkItems;
        VecCharVec mActions;
        Dotter mDotter;
    };

    StackProc( Common& common )
    : mCommon(common), mScorer(common.getParams()) {}
    StackProc( StackProc const& that )
    : mCommon(that.mCommon), mScorer(that.mScorer) {}

    // compiler-supplied destructor is OK

    void operator()( size_t batchNo );

    void dumpStack( WorkItem const& item, long minOffset );

private:
    StackProc& operator=( StackProc const& ); // unimplemented -- no assignment


    long findMinOffset( std::vector<ReadLoc> const& locs )
    { long minOffset = 0;
      typedef std::vector<ReadLoc>::const_iterator Itr;
      for ( Itr itr(locs.begin()), end(locs.end()); itr != end; ++itr )
      { long offset = itr->getOffset();
        offset -= mCommon.getRead(itr->getReadID().val()).size();
        if ( offset < minOffset ) minOffset = offset; }
        return minOffset; }

    // returns true if there are multiple supported calls
    bool processColumn( long colOffset, bool onEdge,
                            LocItr begin, LocItr* pEnd );

    Common& mCommon;
    Scorer mScorer;
};

template <unsigned K>
typename StackProc<K>::Scorer::Results const
    StackProc<K>::Scorer::getResults( bool onEdge )
{
    unsigned supportCount = 0;
    if ( !mParams.mDoBranches )
    {
        for ( unsigned idx = 0; idx < 4; ++idx )
            if ( isSupported(mBaseScores[idx]) )
                supportCount += 1;
    }
    Results results(supportCount > 1);

    unsigned nCalls = mBaseScores[0].mQualCount+mBaseScores[1].mQualCount+
                        mBaseScores[2].mQualCount+mBaseScores[3].mQualCount;
    if ( nCalls >= mParams.mBaseStackDepth )
    {
        unsigned winner = Results::NO_WINNER;
        unsigned winQualSum = 0;
        unsigned placeQualSum = 0;
        for ( unsigned idx = 0; idx < 4; ++idx )
        {
            unsigned qualSum = mBaseScores[idx].mQualSum;
            if ( qualSum > winQualSum )
            {
                placeQualSum = winQualSum;
                winQualSum = qualSum;
                winner = idx;
            }
            else if ( qualSum > placeQualSum )
                placeQualSum = qualSum;
        }
        if ( winQualSum >= mParams.mWinningSum )
        {
            double placeWinRatio = 1. * placeQualSum / winQualSum;
            if ( placeWinRatio <= mParams.mPlaceToWinRatio )
            {
                results.mWinner = winner;
                bool foundUncorrectable = false;
                for ( unsigned idx = 0; idx < 4; ++idx )
                {
                    if ( idx != winner )
                    {
                        if ( !isCorrectable(idx,winQualSum) )
                            foundUncorrectable = true;
                        else if ( mBaseScores[idx].mQualCount )
                            // only mark if correctable and needs correction
                            results.setAction(idx,CORRECT);
                    }
                }
                BaseScore const& bsWinner = mBaseScores[winner];
                if ( !foundUncorrectable &&
                     placeWinRatio <= mParams.mPlaceToWinRatioToConfirm &&
                     (!onEdge || isSupported(bsWinner)) )
                {
                    if ( mParams.mConfirmingSum )
                    {
                        if ( winQualSum >= mParams.mConfirmingSum )
                            results.setAction(winner,CONFIRM);
                    }
                    else if ( bsWinner.mSupportingCount >=
                                                    mParams.mConfirmingCount )
                        results.setAction(winner,CONFIRM);
                }
            }

            if ( mParams.mAutoConfirm )
                for ( unsigned idx = 0; idx < 4; ++idx )
                    if ( mBaseScores[idx].mQualSum >= mParams.mAutoConfirm )
                        results.setAction(idx,CONFIRM);
        }
    }

    clearScores();
    return results;
}

template <unsigned K>
void StackProc<K>::Common::applyActions()
{
    size_t actionCounts[4];
    memset(actionCounts,0,sizeof(actionCounts));
    vecbvec::iterator oBItr(mReads.begin());
    VecQualNibbleVec::iterator oQItr(mQuals.begin());
    for ( VecCharVec::iterator oItr(mActions.begin()), oEnd(mActions.end());
            oItr != oEnd; ++oItr, ++oBItr, ++oQItr )
    {
        bvec& bv = *oBItr;
        QualNibbleVec& qv = *oQItr;
        CharVec& av = *oItr;
        CharVec::iterator beg(av.begin());
        for ( CharVec::iterator itr(beg), end(av.end()); itr != end; ++itr )
        {
            char val = *itr;
            actionCounts[val&3] += 1;
            if ( static_cast<Action>(val&3) == CORRECT )
            {
                size_t idx = itr - beg;
                bv.set(idx,val >> 2);
                qv.set(idx,0);
            }
        }
    }
    ULOG("\nNo action:  " << actionCounts[0]);
    ULOG("Corrected:  " << actionCounts[1]);
    ULOG("Confirmed:  " << actionCounts[2]);
    ULOG("Conflicted: " << actionCounts[3]);
}

template <unsigned K>
void StackProc<K>::operator()( size_t batchNo )
{
    WorkItem const& item = mCommon.getWorkItem(batchNo);
    std::vector<ReadLoc> const& locs = item.getLocs();
    AssertGt(locs.size(),0ul);
    long minOffset = findMinOffset(locs);
    long maxOffset = locs.back().getOffset();
    AssertLt(minOffset,maxOffset);

    //std::cout << "Processing edge " << item.getEdgeID() << std::endl;
    //dumpStack(item,minOffset);

    LocItr begin(locs.begin());
    LocItr end(locs.end());
    UnipathGraph<K> const& graph = mCommon.getGraph();
    UnipathEdge const& edge = graph.getEdge(item.getEdgeID());
    long edgeLen = edge.getLength()+K-1;

    if ( maxOffset > edgeLen && !mCommon.getParams().mDoBranches )
    {
        if ( processColumn(edgeLen,false,begin,&end) )
            maxOffset = edgeLen;
        else
            end = locs.end();
    }
    long offset = maxOffset;
    while ( --offset >= minOffset && begin != end )
    {
        bool onEdge = offset >= 0 && offset < edgeLen;
        if ( processColumn(offset,onEdge,begin,&end) && offset < 0 )
            break;
    }

    mCommon.dot();
}

template <unsigned K>
void StackProc<K>::dumpStack( WorkItem const& item, long minOffset )
{
    std::vector<ReadLoc> const& locs = item.getLocs();
    bvec scratch;
    for ( LocItr itr(locs.begin()), end(locs.end()); itr != end; ++itr )
    {
        size_t readId = itr->getReadID().val();
        bvec const& read = mCommon.getRead(readId);
        size_t skip = itr->getOffset() - read.size() - minOffset;
        bool isRC = itr->isRC();
        std::cout << std::string(skip,' ')
                    << (isRC ? '-' : '+')
                    << (isRC ? scratch.ReverseComplement(read) : read)
                    << ' ' << readId << std::endl;
    }

    UnipathGraph<K> const& graph = mCommon.getGraph();
    std::cout << std::string(-minOffset,' ') << '=';
    std::transform(graph.getBases(item.getEdgeID()),
                   graph.getBasesEnd(item.getEdgeID()),
                   std::ostream_iterator<char>(std::cout),BaseToCharMapper());
    std::cout << std::endl;
}

template <unsigned K>
bool StackProc<K>::processColumn( long colOffset, bool onEdge,
                                    LocItr begin, LocItr* pEnd )
{
    LocItr newEnd(begin);
    LocItr itr(*pEnd);
    do
    {
        ReadLoc const& loc = *--itr;
        long locOffset = loc.getOffset();
        if ( locOffset <= colOffset )
            break;

        size_t readId = loc.getReadID().val();
        bvec const& read = mCommon.getRead(readId);
        long rdOffset = read.size() - (locOffset - colOffset);
        if ( rdOffset >= 0 )
        {
            if ( newEnd == begin )
                newEnd = itr + 1;
            QualNibbleVec const& qual = mCommon.getQual(readId);
            if ( !loc.isRC() )
                mScorer.addScore( read[rdOffset], qual[rdOffset] );
            else
            {
                rdOffset = read.size() - rdOffset - 1;
                unsigned char base = GetComplementaryBase(read[rdOffset]);
                mScorer.addScore( base, qual[rdOffset] );
            }
        }
    }
    while ( itr != begin );

    //std::cout << colOffset << ": " << mScorer;
    typename Scorer::Results const results = mScorer.getResults(onEdge);
    //std::cout << ' ' << results << std::endl;

    if ( results.isActionNecessary() )
    {
        unsigned char winner = results.getWinner();
        unsigned char cWinner = GetComplementaryBase(winner);
        LocItr itr(newEnd);
        do
        {
            ReadLoc const& loc = *--itr;
            long locOffset = loc.getOffset();
            if ( locOffset <= colOffset )
                break;

            size_t readId = loc.getReadID().val();
            bvec const& read = mCommon.getRead(readId);
            long rdOffset = read.size() - (locOffset - colOffset);
            if ( rdOffset >= 0 )
            {
                unsigned char base;
                if ( !loc.isRC() )
                {
                    base = read[rdOffset];
                    Action action = results.getAction(base);
                    if ( action == CONFIRM )
                        mCommon.confirm(readId,rdOffset);
                    else if ( action == CORRECT )
                        mCommon.correct(readId,rdOffset,winner);
                }
                else
                {
                    unsigned rdOffsetRC = read.size() - rdOffset - 1;
                    base = GetComplementaryBase(read[rdOffsetRC]);
                    Action action = results.getAction(base);
                    if ( action == CONFIRM )
                        mCommon.confirm(readId,rdOffsetRC);
                    else if ( action == CORRECT )
                        mCommon.correct(readId,rdOffsetRC,cWinner);
                }
            }
        }
        while ( itr != begin );
    }

    *pEnd = newEnd;
    return results.isBranch();
}

class SquashProc
{
public:
    class Common
    {
    public:
        Common( VecQualNibbleVec& quals, size_t batchSize, size_t nBatches,
                unsigned radius )
        : mQuals(quals), mBatchSize(batchSize), mRadius(radius),
          mDotter(nBatches) {}

        // compiler-supplied copying and destructor are OK

        VecQualNibbleVec& getQuals() { return mQuals; }
        size_t getBatchSize() const { return mBatchSize; }
        unsigned getRadius() const { return mRadius; }

        void dot() { mDotter.batchDone(); }

    private:
        VecQualNibbleVec& mQuals;
        size_t mBatchSize;
        unsigned mRadius;
        Dotter mDotter;
    };

    SquashProc( Common& common ) : mCommon(common) {}

    // compiler-supplied copying and destructor are OK

    void operator()( size_t batchNo );

private:
    Common& mCommon;
};

void SquashProc::operator()( size_t batchNo )
{
    using std::min;
    VecQualNibbleVec& quals = mCommon.getQuals();
    size_t nReads = quals.size();
    size_t batchSize = mCommon.getBatchSize();
    size_t readId = min(nReads,batchNo*batchSize);
    size_t endId = min(nReads,readId+batchSize);
    unsigned radius = mCommon.getRadius();

    typedef VecQualNibbleVec::iterator Itr;
    Itr end(quals.begin(endId));
    for ( Itr itr(quals.begin(readId)); itr != end; ++itr )
        itr->squash(radius);
    mCommon.dot();
}

void squashQuals( VecQualNibbleVec& quals, unsigned nThreads, unsigned radius )
{
    DLOG("Squashing quality scores.");
    size_t batchSize = (quals.size()+nThreads-1)/nThreads;
    SquashProc::Common common(quals,batchSize,nThreads,radius);
    SquashProc proc(common);
    parallelFor(0u,nThreads,proc,nThreads);
}

template <unsigned K>
void buildStacks( ExternalSorter<ReadLineup>::Merger<>& merger,
                  unsigned minStackDepth, size_t estLocs,
                  std::vector<typename StackProc<K>::WorkItem>& workItems )
{
    typedef typename StackProc<K>::WorkItem WorkItem;
    DLOG("Evaluating stacks for inclusion.");
    ReadLineup rl;
    merger.getNext(rl);
    WorkItem curItem(rl.getEdgeID(),estLocs);
    do
    {
        if ( rl.getEdgeID() != curItem.getEdgeID() )
        {
            WorkItem nextItem(rl.getEdgeID(),estLocs);
            if ( curItem.getNLocs() < minStackDepth )
                swap(curItem,nextItem);
            else
            {
                workItems.push_back(nextItem);
                swap(curItem,workItems.back());
            }
        }
        curItem.addLoc(rl.getReadLoc());
    }
    while ( merger.getNext(rl) );
    if ( curItem.getNLocs() >= minStackDepth )
    {
        workItems.push_back(WorkItem(curItem.getEdgeID(),0));
        swap(curItem,workItems.back());
    }
}

template <unsigned K>
void processStacks( UnipathGraph<K> const& graph, unsigned nThreads,
                    typename StackProc<K>::Params const& params, vecbvec& reads,
                    VecQualNibbleVec& quals, ExternalSorter<ReadLineup>& xs )
{
    if ( params.mQualCeilRadius )
        squashQuals(quals,nThreads,params.mQualCeilRadius);

    size_t estLocs = xs.size()/graph.getNEdges();
    ExternalSorter<ReadLineup>::Merger<> merger(xs);
    if ( !merger.size() )
    {
        ULOG("There are no stacks to process");
        return; // EARLY RETURN!
    }

    std::vector<typename StackProc<K>::WorkItem> workItems;
    workItems.reserve(graph.getNEdges());
    buildStacks<K>(merger,params.mReadStackDepth,estLocs,workItems);

    typedef typename StackProc<K>::Common Common;
    Common common(graph,params,reads,quals,workItems);
    if ( true )
    {
        DLOG("Processing " << workItems.size() << " of " << graph.getNEdges()
                << " stacks.");
        StackProc<K> proc(common);
        if ( !nThreads || nThreads > processorsOnline() )
            nThreads = processorsOnline();
        parallelFor(0ul,workItems.size(),proc,nThreads);
    }
    common.applyActions();
}

void useKmerSpectrum( String const& SPECTRUM_SUMMARY_FILE,
                      unsigned& ESTIMATED_COVERAGE,
                      unsigned& MIN_READSTACK_DEPTH,
                      unsigned& MIN_BASESTACK_DEPTH,
                      unsigned& MIN_QSS_TO_CONFIRM,
                      unsigned& MIN_QSS_TO_SUPPORT,
                      unsigned& MIN_QSS_TO_WIN,
                      unsigned& MAX_KMER_FREQ_TO_MARK )
{
    // KmerSpectrum produces a file with summary of the kmer spectrum curve.
    // From KmerSpectra.h:
    // ---- Separate kmer spectrum in 4 regions based on the kf:
    //    1      ... kf_lo    : bad kmers with low frequency
    //    kf_lo  ... kf_hi    : good kmers CN = 1
    //    kf_hi  ... kf_max   : good kmers CN > 1 (repetitive)
    //    kf_max ... inf      : bad kmers with high frequency
    // in addition, we also capture the mode of the spectrum (highest point)

    size_t kf_lo = 0, kf_mode = 0;

    if ( ESTIMATED_COVERAGE )
    {
        kf_lo = round(ESTIMATED_COVERAGE / 7.0);
        kf_mode = round(4.0 * ESTIMATED_COVERAGE / 7.0);
        DLOG("Setting heuristic parameters based on estimated coverage "
             << "(min=" << kf_lo << ", mode=" << kf_mode << ")");
    }
    else if ( SPECTRUM_SUMMARY_FILE != "" &&
              IsRegularFile(SPECTRUM_SUMMARY_FILE) )
    {
        double genome_size = 0;
        size_t kf_hi = 0, kf_max = 0;
        ifstream is(SPECTRUM_SUMMARY_FILE.c_str());
        is >> genome_size;
        is >> ESTIMATED_COVERAGE;
        is >> kf_lo;
        is >> kf_mode;
        is >> kf_hi;
        is >> kf_max;
        is.close();
        DLOG("Setting heuristic parameters based on kmer spectrum in "
              << SPECTRUM_SUMMARY_FILE
              << " (min=" << kf_lo << ", mode=" << kf_mode << ")");
    }
    else
    {
        DLOG("Can't determine scaling parameters.");
        return;
    }

    // Set min read stack 1/4 of the way off the minimum toward the mode (log scale);
    // reasoning is that we only want to error correct off stacks which we believe to
    // consist of valid kmers, and it's about 50/50 at the minimum.
    //double good_stack = pow((double)kf_lo * kf_lo * kf_lo * kf_mode, .25);
    // Scale base determined by looking at examples at our nominal default of 50x coverage

    double good_stack = kf_lo;
    double scale = good_stack / 7.0;

    MIN_QSS_TO_CONFIRM = round(scale * MIN_QSS_TO_CONFIRM);
    ULOG("Setting MIN_QSS_TO_CONFIRM to " << MIN_QSS_TO_CONFIRM);

    MIN_QSS_TO_SUPPORT = round(scale * MIN_QSS_TO_SUPPORT);
    ULOG("Setting MIN_QSS_TO_SUPPORT to " << MIN_QSS_TO_SUPPORT);

    MIN_QSS_TO_WIN = round(scale * MIN_QSS_TO_WIN);
    ULOG("Setting MIN_QSS_TO_WIN to " << MIN_QSS_TO_WIN);

    MIN_READSTACK_DEPTH = round(good_stack);
    ULOG("Setting MIN_READSTACK_DEPTH to " << MIN_READSTACK_DEPTH);

    // Don't go too far into the fringes off the base kmer stack; since we're scaling MIN_QSS...
    // parameters linearly, this should probalby be linear with the stack size as well.
    MIN_BASESTACK_DEPTH = round(good_stack / 4.0);
    ULOG("Setting MIN_BASESTACK_DEPTH to " << MIN_BASESTACK_DEPTH);

    // Set kmer cutoff halfway between 1 and the spectrum minimum (log scale)
    MAX_KMER_FREQ_TO_MARK = round(sqrt(kf_lo));
    ULOG("Setting MAX_KMER_FREQ_TO_MARK to " << MAX_KMER_FREQ_TO_MARK);
}

template <unsigned K>
class KeepProc
{
public:
    class Common
    {
    public:
        Common( UnipathGraph<K> const& graph, vecbvec const& reads,
                size_t batchSize, BitVec& keepers )
        : mGraph(graph), mReads(reads), mBatchSize(batchSize), mKeepers(keepers)
        {}

        UnipathGraph<K> const& getGraph() const { return mGraph; }
        vecbvec const& getReads() const { return mReads; }
        size_t getBatchSize() const { return mBatchSize; }
        BitVec& getKeepers() { return mKeepers; }

    private:
        UnipathGraph<K> const& mGraph;
        vecbvec const& mReads;
        size_t mBatchSize;
        BitVec& mKeepers;
    };

    KeepProc( Common& common ) : mCommon(common) {}

    // compiler-supplied copying and destructor are OK

    void operator()( size_t batchNo );

private:
    bool isKeeper( bvec const& read );

    Common& mCommon;
};

template <unsigned K>
void KeepProc<K>::operator()( size_t batchNo )
{
    using std::min;
    size_t batchSize = mCommon.getBatchSize();
    vecbvec const& reads = mCommon.getReads();
    size_t readId = min(reads.size(),batchNo*batchSize);
    size_t finalReadId = min(reads.size(),readId+batchSize);
    BitVec& keepers = mCommon.getKeepers();
    for ( ; readId < finalReadId; ++readId )
    {
        if ( !isKeeper(reads[readId]) )
            keepers.set(readId,false);
    }
}

template <unsigned K>
bool KeepProc<K>::isKeeper( bvec const& bv )
{
    typedef bvec::const_iterator RItr;
    typedef HugeBVec::const_iterator UItr;
    typedef HugeBVec::const_rc_iterator URCItr;

    UnipathGraph<K> const& graph = mCommon.getGraph();
    size_t rdLen = bv.size();
    if ( rdLen < K )
        return false; // EARLY RETURN!

    KMer<K> kkk(bv.begin());
    KDef const* pDef = graph.getDict().lookup(kkk);
    if ( !pDef )
        return false; // ANOTHER EARLY RETURN!

    RItr rBeg(bv.begin());
    RItr rK(bv.begin(K));
    RItr rEnd(bv.end());
    UItr uBeg(graph.getBases(*pDef));
    URCItr uRCBeg;

    using std::equal;
    using std::min;
    bool isRC = !equal(rBeg,rK,uBeg);

    EdgeID edgeID = pDef->getEdgeID();
    if ( !isRC )
        uBeg += K;
    else
    {
        size_t off = graph.getEdge(edgeID).getLength() -
                        pDef->getEdgeOffset() - 1;
        uRCBeg = graph.getRCBases(edgeID) + off;
        Assert(equal(rBeg,rK,uRCBeg));
        uRCBeg += K;
    }
    rBeg += K;

    unsigned kLess1 = K - 1;

    while ( rBeg != rEnd )
    {
        if ( isRC )
        {
            URCItr uRCEnd(graph.getRCBasesEnd(edgeID));
            RItr itr(rBeg+min(rEnd-rBeg,uRCEnd-uRCBeg));
            if ( !equal(rBeg,itr,uRCBeg) )
                return false; // YET ANOTHER!
            rBeg = itr;
            if ( rBeg != rEnd )
            {
                unsigned char base = *itr;
                unsigned char cBase = GetComplementaryBase(base);
                UnipathEdge const& edge = graph.getEdge(edgeID);
                edgeID = edge.getPredecessor(cBase);
                if ( edgeID.isNull() )
                    return false; // AND YET ANOTHER!
                isRC = !edge.isPredecessorRC(cBase);
                if ( isRC )
                {
                    uRCBeg = graph.getRCBases(edgeID)+kLess1;
                    AssertEq(*uRCBeg,base);
                }
                else
                {
                    uBeg = graph.getBases(edgeID)+kLess1;
                    AssertEq(*uBeg,base);
                }
            }
        }
        else
        {
            UItr uEnd(graph.getBasesEnd(edgeID));
            RItr itr(rBeg+min(rEnd-rBeg,uEnd-uBeg));
            if ( !equal(rBeg,itr,uBeg) )
                return false; // EARLY RETURNS ALL OVER THE PLACE!
            rBeg = itr;
            if ( rBeg != rEnd )
            {
                unsigned char base = *itr;
                UnipathEdge const& edge = graph.getEdge(edgeID);
                edgeID = edge.getSuccessor(base);
                if ( edgeID.isNull() )
                    return false; // IT NEVER ENDS!
                isRC = edge.isSuccessorRC(base);
                if ( isRC )
                {
                    uRCBeg = graph.getRCBases(edgeID)+kLess1;
                    AssertEq(*uRCBeg,base);
                }
                else
                {
                    uBeg = graph.getBases(edgeID)+kLess1;
                    AssertEq(*uBeg,base);
                }
            }
        }
    }
    return true;
}

template <unsigned K>
size_t findKeepers( UnipathGraph<K> const& graph, vecbvec const& reads,
                    unsigned nThreads, BitVec& keepers )
{
    DLOG("Re-pathing to identify reads with low-frequency kmers.");

    size_t nReads = reads.size();
    keepers.clear().resize(nReads,true);
    size_t batchSize = (nReads+nThreads-1)/nThreads;

    if ( true )
    {
        typename KeepProc<K>::Common common(graph,reads,batchSize,keepers);
        KeepProc<K> proc(common);
        parallelFor(0u,nThreads,proc,nThreads);
    }

    size_t nGood = keepers.Sum();
    size_t nBad = nReads - nGood;
    double pctBad = 100.*nBad/nReads;
    ULOG("There are " << nBad << " (" << pctBad << "%) uncorrectable reads.");
    return nGood;
}

void writeCorrectedFiles( vecbvec const& reads,
                         QualNibbleVecVec const& quals,
                         BitVec const& keepers,
                         size_t nKeepers,
                         String const& inHead,
                         String const& editHead )
{
    vec<Bool> losers;
    losers.reserve(reads.size());
    IncrementalWriter<bvec> readWriter((editHead+".fastb").c_str(),nKeepers);
    IncrementalWriter<qvec> qualWriter((editHead+".qualb").c_str(),nKeepers);
    typedef BitVec::const_iterator Itr;
    for ( Itr itr(keepers.cbegin()), end(keepers.cend()); itr != end; ++itr )
    {
        bool keep = *itr;
        losers.push_back(!keep);
        if ( keep )
        {
            unsigned idx = itr.pos();
            readWriter.add(reads[idx]);
            qualWriter.add(quals[idx].GetQualvector());
        }
    }
    readWriter.close();
    qualWriter.close();

    PairsManager pairs(inHead+".pairs");
    size_t nPairsThen = pairs.nPairs();
    size_t nUnpairedThen = pairs.nReads() - 2*nPairsThen;
    pairs.removeReads( losers, true );
    pairs.Write(editHead+".pairs");
    size_t nPairsNow = pairs.nPairs();
    size_t nUnpairedNow = pairs.nReads() - 2*nPairsNow;
    ULOG("After correction there are " << nPairsNow << " of " << nPairsThen
            << " pairs that survived.");
    ULOG("There are " << nUnpairedNow-nUnpairedThen << " additional orphans.");
}

}

TRIVIALLY_SERIALIZABLE(ReadLineup);

int main(int argc, char *argv[])
{
    String empty;
    RunTime();
    BeginCommandArguments;

    CommandArgument_String_OrDefault_Doc(IN_HEAD,"frag_reads_prec",
        "Looks for IN_HEAD.{fastb,qualb,pairs}");
    CommandArgument_String_OrDefault_Doc(OUT_EDIT_HEAD,"frag_reads_edit",
        "Generates OUT_EDIT_HEAD.{fastb,qualb,keep}");
    CommandArgument_String_OrDefault_Doc(OUT_CORR_HEAD,"frag_reads_corr",
        "If DELETE is TRUE, generates OUT_CORR_HEAD.{fastb,qualb,pairs}");
    CommandArgument_String_OrDefault_Doc(WORKDIR,".",
        "Place to write scratch files.");

    CommandArgument_String_OrDefault_Doc(USE_KMER_SPECTRUM, "",
        "If set to name of a spectrum_summary file created by KmerSpectrum, "
        "try to automatically adjust hueristic parameters based on spectrum "
        "features. Modifies MIN_READSTACK_DEPTH, MIN_BASESTACK_DEPTH, "
        "MAX_KMER_FREQUENCY_TO_MARK, MIN_QSS_TO_WIN, MIN_QSS_TO_CONFIRM and "
        "MIN_QSS_TO_SUPPORT by scaling factors based on kmer spectrum.");

    CommandArgument_UnsignedInt_OrDefault_Doc(AUTO_CONFIRM, 0,
        "Automatically confirm bases having this quality score sum. "
        "(Defaults to no automatic confirmation.)");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_READSTACK_DEPTH, 5,
        "Minimum number of reads in a valid stack.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_BASESTACK_DEPTH, 1,
        "Minimum number of bases in a stack column used for inferences.");
    CommandArgument_UnsignedInt_OrDefault_Doc(QUAL_CEIL_RADIUS, 2,
        "Replace each quality score by the minimum of the quality scores "
        "over a window of this radius.");
    CommandArgument_Bool_OrDefault_Doc(QCR_ALWAYS, True,
        "Apply QUAL_CEIL_RADIUS every cycle. Doesn't really make sense to do "
        "it more than once but may improve results.");
    CommandArgument_Bool_OrDefault_Doc(DO_BRANCHES, False,
       "Correct errors beyond unipath boundaries despite finding supported "
       "branching.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_QSS_TO_WIN, 40,
        "Minimum quality score sum for best call to allow a stack column "
        "to be used for inferences.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_QSS_TO_SUPPORT, 60,
        "Minimum quality score sum to be considered supported.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_QSS_TO_CONFIRM, 90,
        "Minimum quality score sum to lock the call.");
    CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CORRECT, 0.25,
        "Maximum ratio of correctable quality score sum and best QSS.");
    CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CORRECT2, 0.25,
        "Maximum ratio of second best to best quality score sum allowing "
        "any corrections.");
    CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CONFIRM, 1.0,
        "Maximum ratio of second best to best quality score sum to confirm "
        "any base.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_KMER_FREQ_TO_MARK, 2,
        "Reads with minimum kmer frequencies at and below this will be marked "
        "and discarded in _corr.");

    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_CYCLES, 2,
        "Number of correction cycles.");
    CommandArgument_Bool_OrDefault_Doc(DELETE, False,
       "actually delete reads, rather than mark them for deletion");

    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
            "Number of parallel threads, defaults to number of processors.");
    CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);

    EndCommandArguments;

    unsigned const K = 28;

    SetMaxMemory(MAX_MEMORY_GB<<30);
    if ( NUM_THREADS < 1 || NUM_THREADS > processorsOnline() )
        NUM_THREADS = processorsOnline();

    // some things that could be parameters, but which RunAllPathsLG doesn't
    // offer an opportunity to adjust, so they're not parameters right now
    unsigned ESTIMATED_COVERAGE = 0;
    unsigned MIN_Q_TO_SUPPORT = 20;
    unsigned MIN_N_BASES_TO_SUPPORT = 2;
    unsigned MIN_N_BASES_TO_CONFIRM = 3;
    bool DO_PALINDROMES = false;

    if ( USE_KMER_SPECTRUM != empty || ESTIMATED_COVERAGE )
        useKmerSpectrum(USE_KMER_SPECTRUM, ESTIMATED_COVERAGE,
                        MIN_READSTACK_DEPTH, MIN_BASESTACK_DEPTH,
                        MIN_QSS_TO_CONFIRM, MIN_QSS_TO_SUPPORT,
                        MIN_QSS_TO_WIN, MAX_KMER_FREQ_TO_MARK);

    String fastbIn = IN_HEAD + ".fastb";
    DLOG("Reading reads.");
    vecbvec reads(fastbIn);

    DLOG("Building kmer dictionary, excluding low-frequency kmers.");
    KmerDict<K>* pDict = createDict<K>(reads,NUM_THREADS,MIN_READSTACK_DEPTH,
                                       ESTIMATED_COVERAGE?ESTIMATED_COVERAGE:35,
                                       WORKDIR);
    UnipathGraph<K>* pGraph = new UnipathGraph<K>(*pDict);

    DLOG("Reading quals.");
    VecQualNibbleVec quals;
    LoadQualNibbleVec(fastbIn.ReplaceExtension(".fastb",".qualb"),&quals);

    StackProc<K>::Params params(
            MIN_Q_TO_SUPPORT, MIN_N_BASES_TO_SUPPORT, MAX_QSS_RATIO_TO_CORRECT,
            MAX_QSS_RATIO_TO_CORRECT2, MAX_QSS_RATIO_TO_CONFIRM,
            MIN_QSS_TO_CONFIRM, MIN_QSS_TO_WIN, MIN_READSTACK_DEPTH,
            MIN_BASESTACK_DEPTH, AUTO_CONFIRM, MIN_QSS_TO_SUPPORT,
            MIN_N_BASES_TO_CONFIRM,QUAL_CEIL_RADIUS,DO_BRANCHES);

    for ( unsigned cycle = 0; cycle < NUM_CYCLES; ++cycle )
    {
        DLOG("Starting cycle " << cycle+1);
        ExternalSorter<ReadLineup> xs(WORKDIR+"/FELineups.tmp.");
        createReadLineups(*pGraph,reads,NUM_THREADS,DO_PALINDROMES,&xs);

        processStacks(*pGraph,NUM_THREADS,params,reads,quals,xs);
        if ( !QCR_ALWAYS )
            params.mQualCeilRadius = 0;
    }

    BitVec keepers;
    size_t nKeepers = findKeepers(*pGraph,reads,NUM_THREADS,keepers);

    delete pGraph;
    delete pDict;

    DLOG("Writing output files.");
    reads.WriteAll(OUT_EDIT_HEAD+".fastb");
    WriteAll(quals,OUT_EDIT_HEAD+".qualb");
    String keepFile(OUT_EDIT_HEAD+".keep");
    BinaryWriter::writeFile(keepFile.c_str(),keepers);
    if ( DELETE )
        writeCorrectedFiles(reads,quals,keepers,nKeepers,IN_HEAD,OUT_CORR_HEAD);

    DLOG("Done.");
}
