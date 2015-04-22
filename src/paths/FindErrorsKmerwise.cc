///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FindErrorsKmerwise.cc
 * \author tsharpe
 * \date Jul 7, 2011
 *
 * \brief
 */
#include "MainTools.h"
#include "Basevector.h"
#include "Charvector.h"
#include "FeudalMimic.h"
#include "PairsManager.h"
#include "feudal/BitVec.h"
#include "feudal/QualNibbleVec.h"
#include "kmers/ReadPatherDefs.h"
#include "reporting/PerfStat.h"
#include "system/ExternalSorter.h"
#include "system/Worklist.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <utility>
#include <vector>

#define ULOG(x) std::cout << x << std::endl
#define DLOG(x) ULOG(Date() << " (FE): " << x)

namespace
{
unsigned const K = 24;

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
    else if ( SPECTRUM_SUMMARY_FILE != "" && IsRegularFile(
            SPECTRUM_SUMMARY_FILE) )
    {
        DLOG("Reading spectrum summary from " << SPECTRUM_SUMMARY_FILE);
        ifstream is(SPECTRUM_SUMMARY_FILE.c_str());
        double genome_size;
        is >> genome_size;

        is >> ESTIMATED_COVERAGE;
        is >> kf_lo;
        is >> kf_mode;

        //size_t kf_hi;
        //is >> kf_hi;
        //size_t kf_max;
        //is >> kf_max;

        is.close();
        DLOG("Setting heuristic parameters based on kmer spectrum "
                << "(min=" << kf_lo << ", mode=" << kf_mode << ")");
    }
    else
    {
        DLOG("Can't determine scaling parameters.");
        return;
    }

    // Set min read stack 1/4 (log scale) of the way off the minimum toward the
    // mode.
    // Reasoning is that we only want to error correct off stacks which we
    // believe to consist of valid kmers, and it's about 50/50 at the minimum.
    // double good_stack = pow((double)kf_lo * kf_lo * kf_lo * kf_mode, .25);
    // Scale base determined by looking at examples at our nominal default of
    // 50x coverage

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

    // Don't go too far into the fringes off the base kmer stack.
    // Since we're scaling MIN_QSS parameters linearly, this should probably be
    // linear with the stack size as well.
    MIN_BASESTACK_DEPTH = round(good_stack / 4.0);
    ULOG("Setting MIN_BASESTACK_DEPTH to " << MIN_BASESTACK_DEPTH);

    // Set kmer cutoff halfway between 1 and the spectrum minimum (log scale)
    MAX_KMER_FREQ_TO_MARK = round(sqrt(kf_lo));
    ULOG("Setting MAX_KMER_FREQ_TO_MARK to " << MAX_KMER_FREQ_TO_MARK);
}

class ReadLoc
{
public:
    ReadLoc() : mOffset(0ul), mIsRC(false), mID(0ul) {}
    ReadLoc( size_t id, size_t offset, bool isRC )
    : mOffset(offset), mIsRC(isRC), mID(id) { ForceAssertLt(offset,256ul); }
    ReadLoc( ReadLoc const& that )
    { memcpy(this,&that,sizeof(ReadLoc)); }

    // compiler-supplied destructor OK

    ReadLoc& operator=( ReadLoc const& that )
    { memcpy(this,&that,sizeof(ReadLoc)); return *this; }

    bool isRC() const { return mIsRC; }
    void setRC( bool isRC ) { mIsRC = isRC; }

    size_t getID() const { return mID.val(); }
    void setID( size_t id ) { mID.setVal(id); }

    size_t getOffset() const { return mOffset; }
    void setOffset( size_t offset )
    { ForceAssertLt(offset,256ul); mOffset = offset; }

    friend std::ostream& operator<<( std::ostream& os, ReadLoc const& readLoc )
    { os << '@' << readLoc.getID() << ':' << readLoc.getOffset()
         << (readLoc.isRC() ? '-' : '+');
      return os; }

    friend int compare( ReadLoc const& rl1, ReadLoc const& rl2 )
    { return compare(rl1.mOffset,rl2.mOffset); }

    friend bool operator<( ReadLoc const& rl1, ReadLoc const& rl2 )
    { return rl1.mOffset < rl2.mOffset; }

    friend bool operator<=( ReadLoc const& rl1, ReadLoc const& rl2 )
    { return rl1.mOffset <= rl2.mOffset; }

    friend bool operator>( ReadLoc const& rl1, ReadLoc const& rl2 )
    { return rl1.mOffset > rl2.mOffset; }

    friend bool operator>=( ReadLoc const& rl1, ReadLoc const& rl2 )
    { return rl1.mOffset >= rl2.mOffset; }

    friend bool operator==( ReadLoc const& rl1, ReadLoc const& rl2 )
    { return !memcmp(&rl1,&rl2,sizeof(ReadLoc)); }

    friend bool operator!=( ReadLoc const& rl1, ReadLoc const& rl2 )
    { return memcmp(&rl1,&rl2,sizeof(ReadLoc)); }

private:
    unsigned char mOffset;
    bool  mIsRC;
    ID<5> mID;
};

class ReadKmer : public KMer<K>
{
public:
    ReadKmer() {}
    ReadKmer( KMer<K> const& kmer, ReadLoc const& readLoc )
    : KMer<K>(kmer), mReadLoc(readLoc) {}

    // compiler-supplied copying and destructor are OK

    ReadLoc const& getReadLoc() const { return mReadLoc; }
    ReadLoc& getReadLoc() { return mReadLoc; }

private:
    ReadLoc mReadLoc;
};


typedef std::vector<ReadKmer> VecReadKmer;

// Creates a sorted list of kmers and their locations within reads.
class SortProc
{
public:
    class Common
    {
    public:
        Common( vecbvec const& reads, size_t nBatches, size_t batchSize,
                size_t kmersPerRead, bool doPalindromes,
                ExternalSorter<ReadKmer>* pXS )
        : mReads(reads), mDotter(nBatches), mBatchSize(batchSize),
          mKmersPerBatch(batchSize*kmersPerRead), mDoPalindromes(doPalindromes),
          mXS(*pXS) {}

        // compiler-supplied destructor is OK

        vecbvec const& getReads() const { return mReads; }
        size_t getNBatches() const { return mDotter.getNBatches(); }
        size_t getBatchSize() const { return mBatchSize; }
        size_t getKmersPerBatch() const { return mKmersPerBatch; }
        bool getDoPalindromes() const { return mDoPalindromes; }

        void dump( VecReadKmer& readKmers )
        { if ( !readKmers.empty() )
          { ReadKmer* start = &*readKmers.begin();
            ReadKmer* end = &*readKmers.end();
            sortInPlace(start,end,CompareFunctor<KMer<K> >());
            mXS.insert(start,end); }
          mDotter.batchDone(); }

    private:
        Common( Common const& ); // unimplemented -- no copying
        Common& operator=( Common const& ); // unimplemented -- no copying

        vecbvec const& mReads;
        Dotter mDotter;
        size_t mBatchSize;
        size_t mKmersPerBatch;
        bool mDoPalindromes;
        ExternalSorter<ReadKmer>& mXS;
    };

    class OItr
    {
    public:
        OItr( VecReadKmer& readKmers, size_t readId, size_t readLen,
                bool doPalindromes )
        : mReadKmers(readKmers), mReadLoc(readId,0,false),
          mRevOff(readLen-K), mOffset(0),
          mDoPalindromes(doPalindromes) {}

        // compiler-supplied copying and destructor are OK

        OItr& operator*() { return *this; }
        OItr& operator++() { return *this; }
        OItr& operator++(int) { return *this; }

        KMer<K> const& operator=( KMer<K> const& kmer )
        { CanonicalForm form = kmer.getCanonicalForm();
          if ( mDoPalindromes || form != CanonicalForm::PALINDROME )
          { if ( form != CanonicalForm::REV )
            { mReadLoc.setOffset(mOffset++);
              mReadLoc.setRC(false);
              mReadKmers.push_back(ReadKmer(kmer,mReadLoc)); }
            else
            { mReadLoc.setOffset(mRevOff-mOffset++);
              mReadLoc.setRC(true);
              ReadKmer rk(kmer,mReadLoc); rk.rc();
              mReadKmers.push_back(rk); } }
          return kmer; }

    private:
        VecReadKmer& mReadKmers;
        ReadLoc mReadLoc;
        size_t mRevOff;
        size_t mOffset;
        bool mDoPalindromes;
    };

    SortProc( Common& common, unsigned K )
    : mCommon(common), mReadKmers(K)
    {}

    SortProc( SortProc const& that )
    : mCommon(that.mCommon), mReadKmers(that.mReadKmers)
    { mReadKmers.reserve(mCommon.getKmersPerBatch()); }

    // compiler-supplied destructor is OK, assignment prohibited by references

    void operator()( size_t batchNo );

private:
    Common& mCommon;
    VecReadKmer mReadKmers;
};

void SortProc::operator()( size_t batchNo )
{
    using std::min;
    vecbvec const& reads = mCommon.getReads();
    size_t nReads = reads.size();
    size_t batchSize = mCommon.getBatchSize();
    size_t readId = min(nReads,batchNo*batchSize);
    size_t endId = min(nReads,readId+batchSize);
    for ( ; readId != endId; ++ readId )
    {
        bvec const& bv = reads[readId];
        OItr oItr(mReadKmers,readId,bv.size(),mCommon.getDoPalindromes());
        KMer<K>::kmerizeNonCanonically(bv.begin(),bv.end(),oItr);
    }
    mCommon.dump(mReadKmers);
}

void kmerizeReads( vecbvec const& reads, unsigned nThreads,
                    bool doPalindromes, ExternalSorter<ReadKmer>& xs )
{
    size_t nReads = reads.size();
    size_t kmersPerRead = (reads.SizeSum()+nReads-1)/nReads - K + 1;
    size_t bytesPerUnit = kmersPerRead*sizeof(ReadKmer);
    WorklistParameterizer wlp(reads.size(),bytesPerUnit,10000,nThreads);
    DLOG("Kmerizing " << nReads << " reads in " << wlp.getNBatches()
            << " batches of " << wlp.getBatchSize() << '.');

    SortProc::Common common(reads,wlp.getNBatches(),wlp.getBatchSize(),
                                kmersPerRead,doPalindromes,&xs);
    SortProc proc(common,K);

    parallelFor(0ul,wlp.getNBatches(),proc,wlp.getNThreads());
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

class StackProc
{
public:
    typedef SmallVec< ReadLoc, MempoolAllocator<ReadLoc> > ReadStack;
    typedef ReadStack::iterator LocItr;
    typedef OuterVec<ReadStack> Workitem;

    enum Action { NO_ACTION, CORRECT, CONFIRM, CONFLICT };

    class Scorer
    {
    public:
        Scorer( unsigned supportingScore, unsigned supportingCount,
                double loserToWinRatio, double placeToWinRatio,
                double placeToWinRatioToConfirm, unsigned confirmingSum,
                unsigned winningSum, unsigned readStackDepth,
                unsigned baseStackDepth, unsigned autoConfirm,
                unsigned supportingSum, unsigned confirmingCount )
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
          mBaseStackDepth(baseStackDepth)
        { clearScores(); }

        // compiler-supplied copying and destructor are OK

        unsigned getMinReadStackDepth() const { return mReadStackDepth; }
        unsigned getMinBaseStackDepth() const { return mBaseStackDepth; }

        void addScore( unsigned char base, unsigned qual )
        { BaseScore& bs = mBaseScores[base];
          bs.mQualSum += qual; bs.mQualCount += 1;
          if ( qual >= mSupportingScore ) bs.mSupportingCount += 1; }

        class Results
        {
        public:

            Results()
            : mWinner(NO_WINNER), mAction(0), mIsBranchPoint(false) {}

            // compiler-supplied copying and destructor are OK

            bool isBranchPoint() const { return mIsBranchPoint; }
            void setBranchPoint() { mIsBranchPoint = true; }

            Action getAction( unsigned char base ) const
            { return static_cast<Action>((mAction >> (2*base)) & 3); }

            void setAction( unsigned char base, Action action )
            { mAction ^= (getAction(base) ^ action) << (2*base); }

            bool isActionNecessary() const { return mAction; }

            unsigned char getWinner() const { return mWinner; }
            void setWinner( unsigned char winner ) { mWinner = winner; }

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
            unsigned char mWinner;
            unsigned char mAction;
            bool mIsBranchPoint;
        };

        Results const getResults( bool inKmer );

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
        { return mSupportingSum ?
                    bs.mQualSum >= mSupportingSum :
                    bs.mSupportingCount >= mSupportingCount; }

        bool isConfirmable( BaseScore const& bs ) const
        { return bs.mSupportingCount >= mConfirmingCount ||
                    (mConfirmingSum && bs.mQualSum >= mConfirmingSum); }

        void clearScores() { memset(mBaseScores,0,sizeof(mBaseScores)); }

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
        double mPlaceToWinRatio; // (2nd_largest_QSS/largest_QSS) must not exceed this

        // Parameters that determine whether there are loser (i.e., correctable)
        // calls in the current column.  In order to be corrected, there must be
        // a (different) winning call, the call must not be supported, and the
        // QSS must be sufficiently small relative to the winner.
        double mLoserToWinRatio; // loser_QSS/largest_QSS must not exceed this

        // Parameters that determine whether a call is confirmed.  Calls
        // confirmed by a column will not be corrected, even if they are losers
        // in that, or in some other, column.  There are two ways to be
        // confirmed: 1) the QSS must exceed an auto-confirm parameter, or
        unsigned mAutoConfirm; // auto-confirm QSS (0 means ignore this criterion)
        // 2) a call must be the winner, must be supported if current column is
        // within the kmer boundaries, no other call can be uncorrectable, and...
        double mPlaceToWinRatioToConfirm; // (2nd_largest_QSS/largest_QSS) must not exceed this
        unsigned mConfirmingSum; // the QSS must be at least this big (0 means ignore this criterion)
        unsigned mConfirmingCount; // no. of high-q calls needed to confirm

        unsigned mReadStackDepth; // minimum no. of reads in a stack to bother with
        unsigned mBaseStackDepth; // minimum no. of calls in a column to bother with

        BaseScore mBaseScores[4];
    };

    class Common
    {
    public:
        Common( size_t cycle, unsigned K, bool doBranches, bool skipUnder,
                vecbvec& reads, VecQualNibbleVec& quals )
        : mCycle(cycle), mK(K), mNoBranches(!doBranches), mSkipUnder(skipUnder),
          mReads(reads), mQuals(quals), mMaxReadLen(reads.MaxSize())
        { Mimic(reads,mActions); }

        // compiler-supplied destructor is OK

        unsigned getK() const { return mK; }
        bool noBranches() const { return mNoBranches; }
        bool skipUnder() const { return mSkipUnder; }
        bvec const& getRead( size_t idx ) const { return mReads[idx]; }
        size_t getMaxReadLen() const { return mMaxReadLen; }
        QualNibbleVec const& getQual( size_t idx ) const { return mQuals[idx]; }

        void confirm( size_t readId, unsigned rdOffset )
        { mActions[readId][rdOffset] = CONFIRM; }

        void correct( size_t readId, unsigned rdOffset, unsigned char winner )
        { char volatile* pOldVal = &mActions[readId][rdOffset];
          while ( true )
          { char oldVal = *pOldVal;
            Action oldAction = static_cast<Action>(oldVal&3);
            if ( oldAction == CONFIRM || oldAction == CONFLICT ) break;
            char newVal = CORRECT | (winner << 2);
            if ( oldAction == CORRECT )
            { if ( oldVal == newVal ) break;
              newVal = CONFLICT; }
            if ( __sync_bool_compare_and_swap(pOldVal,oldVal,newVal) )
              break; } }

        void applyActions();

    private:
        Common( Common const& ); // unimplemented -- no copying
        Common& operator=( Common const& ); // unimplemented -- no copying

        size_t mCycle;
        unsigned mK;
        bool mNoBranches;
        bool mSkipUnder;
        vecbvec& mReads;
        VecQualNibbleVec& mQuals;
        size_t mMaxReadLen;
        VecCharVec mActions;
    };

    StackProc( Common& common, Scorer const& scorer )
    : mCommon(common), mScorer(scorer) {}
    StackProc( StackProc const& that )
    : mCommon(that.mCommon), mScorer(that.mScorer) {}

    // compiler-supplied destructor is OK

    void operator()( Workitem& );
    void dumpStack( ReadStack const& locs );

private:
    StackProc& operator=( StackProc const& ); // unimplemented -- no assignment

    // returns true if there are multiple supported calls (i.e., a branch point)
    // and we're looking for branch points
    bool processColumn( long column, LocItr begin, LocItr end );

    Common& mCommon;
    Scorer mScorer;

    // set by processColumn to indicate productive parts of the ReadLocs vector
    // on the last trip through.  depending on whether column is ascending or
    // descending, you can use one of these to trim for the next trip.
    LocItr mStart;
    LocItr mStop;
};

StackProc::Scorer::Results const StackProc::Scorer::getResults( bool inKmer )
{
    Results results;

    unsigned nCalls = mBaseScores[0].mQualCount+mBaseScores[1].mQualCount+
                        mBaseScores[2].mQualCount+mBaseScores[3].mQualCount;
    if ( nCalls >= mBaseStackDepth )
    {
        unsigned winner = Results::NO_WINNER;
        unsigned winQualSum = 0;
        unsigned placeQualSum = 0;
        unsigned support = 0;
        unsigned supportCount = 0;
        for ( unsigned idx = 0; idx < 4; ++idx )
        {
            BaseScore const& bs = mBaseScores[idx];
            if ( isSupported(bs) )
            {
                support |= 1<<idx;
                supportCount += 1;
            }

            unsigned qualSum = bs.mQualSum;
            if ( qualSum > winQualSum )
            {
                placeQualSum = winQualSum;
                winQualSum = qualSum;
                winner = idx;
            }
            else if ( qualSum > placeQualSum )
                placeQualSum = qualSum;
        }
        if ( supportCount > 1 )
            results.setBranchPoint();

        if ( winQualSum >= mWinningSum )
        {
            double placeWinRatio = 1. * placeQualSum / winQualSum;
            double biggestLoserRatio = 0.;
            if ( placeWinRatio <= mPlaceToWinRatio )
            {
                results.setWinner(winner);
                for ( unsigned idx = 0; idx < 4; ++idx )
                {
                    if ( idx != winner )
                    {
                        BaseScore const& bs = mBaseScores[idx];
                        if ( bs.mQualCount && !isSupported(bs) )
                        {
                            double loserRatio = 1.*bs.mQualSum/winQualSum;
                            if ( loserRatio > biggestLoserRatio )
                                biggestLoserRatio = loserRatio;
                            if ( loserRatio <= mLoserToWinRatio )
                                results.setAction(idx,CORRECT);
                        }
                    }
                }
            }

            unsigned winMask = 1u << winner;
            if ( !inKmer )
                support |= winMask;
            if ( support == winMask &&
                    placeWinRatio <= mPlaceToWinRatioToConfirm &&
                    biggestLoserRatio <= mLoserToWinRatio &&
                    isConfirmable(mBaseScores[winner]) )
                results.setAction(winner,CONFIRM);

            if ( mAutoConfirm )
                for ( unsigned idx = 0; idx < 4; ++idx )
                    if ( mBaseScores[idx].mQualSum >= mAutoConfirm )
                        results.setAction(idx,CONFIRM);
        }
    }

    clearScores();
    return results;
}

void StackProc::Common::applyActions()
{
    size_t actionCounts[4];
    memset(actionCounts,0,sizeof(actionCounts));
    vecbvec::iterator rItr(mReads.begin());
    VecQualNibbleVec::iterator qItr(mQuals.begin());
    typedef VecCharVec::iterator Itr;
    Itr aEnd(mActions.end());
    for ( Itr aItr(mActions.begin()); aItr != aEnd; ++aItr, ++rItr, ++qItr )
    {
        bvec& bv = *rItr;
        QualNibbleVec& qv = *qItr;
        CharVec& av = *aItr;
        CharVec::iterator beg(av.begin());
        for ( CharVec::iterator itr(beg), end(av.end()); itr != end; ++itr )
        {
            Action action = static_cast<Action>(*itr & 3);
            actionCounts[action] += 1;
            if ( action == CORRECT )
            {
                size_t idx = itr - beg;
                bv.set(idx,*itr >> 2);
                qv.set(idx,0);
            }
        }
    }
    ULOG("\nNo action:  " << actionCounts[0]);
    ULOG("Confirmed:  " << actionCounts[2]);
    ULOG("Corrected:  " << actionCounts[1]);
    ULOG("Conflicted: " << actionCounts[3]);
    size_t nBases = actionCounts[0] + actionCounts[1] +
                        actionCounts[2] + actionCounts[3];
    String cycle = ToString(mCycle);
    PerfStat::log() << std::fixed << std::setprecision(1)
                    << PerfStat("frac_bases_confirmed_" + cycle,
                                "% of bases confirmed in cycle " + cycle,
                                100.0 * actionCounts[2] / nBases);
    PerfStat::log() << std::fixed << std::setprecision(2)
                    << PerfStat("frac_bases_corrected_" + cycle,
                                "% of bases corrected in cycle " + cycle,
                                100.0 * actionCounts[1] / nBases);
    PerfStat::log() << std::fixed << std::setprecision(2)
                    << PerfStat("frac_conflicts_" + cycle,
                    "% of bases with conflicting corrections in cycle " + cycle,
                                100.0 * actionCounts[3] / nBases);
}

void StackProc::operator()( Workitem& wi )
{
    typedef Workitem::iterator Itr;
    for ( Itr itr(wi.begin()), end(wi.end()); itr != end; ++itr )
    {
        LocItr start(itr->begin());
        LocItr stop(itr->end());
        if ( stop-start < mScorer.getMinReadStackDepth() )
            continue; // NON-STRUCTURED!

        using std::sort;
        sort( start, stop );

        long idx;
        if ( mCommon.skipUnder() )
            idx = 0; // just do the columns that precede the kmer
        else
            idx = mCommon.getK(); // do the kmer, too
        long minColIdx = 0 - stop[-1].getOffset();
        long minBaseStack = mScorer.getMinBaseStackDepth();
        mStart = start;
        while ( --idx >= minColIdx && stop-mStart >= minBaseStack )
            if ( processColumn(idx,mStart,stop) )
                break;

        // do the part after the kmer
        idx = mCommon.getK();
        long endColIdx = mCommon.getMaxReadLen();
        mStop = stop;
        while ( idx < endColIdx && mStop-start >= minBaseStack )
            if ( processColumn(idx++,start,mStop) )
                break;
    }
}

void StackProc::dumpStack( ReadStack const& locs )
{
    if ( locs.empty() )
        return; // EARLY RETURN!

    bvec scratch;
    ReadStack::const_iterator end(locs.end());
    size_t maxOff = end[-1].getOffset();
    std::string kmerXs(std::string(maxOff,' ')+std::string(mCommon.getK(),'x'));
    std::cout << '\n' << kmerXs << '\n';
    for ( ReadStack::const_iterator itr(locs.begin()); itr != end; ++itr )
    {
        size_t readId = itr->getID();
        bvec const& read = mCommon.getRead(readId);
        size_t skip = maxOff - itr->getOffset();
        bool isRC = itr->isRC();
        std::cout << std::string(skip,' ')
                    << (isRC ? scratch.ReverseComplement(read) : read)
                    << ' ' << (isRC ? '-' : '+') << readId << ':'
                    << itr->getOffset() << '\n';
    }
    std::cout << kmerXs << '\n' << std::endl;
}

bool StackProc::processColumn( long colIdx, LocItr begin, LocItr end )
{
    mStart = end;
    mStop = begin;

    LocItr itr(begin);
    while ( itr != end )
    {
        ReadLoc const& rl = *itr;
        long readIdx = colIdx + rl.getOffset();
        if ( readIdx < 0 )
        {
            ++itr;
            continue; // Non-structured!
        }
        if ( readIdx > 0 && mStart == end )
            mStart = itr; // note that we haven't yet incremented itr
        ++itr;

        size_t readId = rl.getID();
        bvec const& read = mCommon.getRead(readId);
        long readLen = read.size();
        if ( readIdx < readLen )
        {
            QualNibbleVec const& qual = mCommon.getQual(readId);
            if ( rl.isRC() )
                mScorer.addScore(*read.rcbegin(readIdx),*qual.rbegin(readIdx));
            else
                mScorer.addScore(read[readIdx],qual[readIdx]);
            mStop = itr; // note that we've already incremented itr
        }
    }

    bool inKmer = static_cast<size_t>(colIdx) < mCommon.getK();
    Scorer::Results const results = mScorer.getResults(inKmer);

    if ( results.isActionNecessary() )
    {
        unsigned char winner = results.getWinner();
        unsigned char cWinner = Scorer::Results::NO_WINNER;
        if ( winner != Scorer::Results::NO_WINNER )
            cWinner = GetComplementaryBase(winner);

        itr = begin;
        while ( itr != end )
        {
            ReadLoc const& rl = *itr;
            ++itr;
            long readIdx = colIdx + rl.getOffset();
            if ( readIdx < 0 )
                continue; // Non-structured!
            size_t readId = rl.getID();
            bvec const& read = mCommon.getRead(readId);
            size_t readLen = read.size();
            if ( readIdx < static_cast<long>(readLen) )
            {
                if ( rl.isRC() )
                {
                    readIdx = readLen - readIdx - 1;
                    unsigned char base = GetComplementaryBase(read[readIdx]);
                    Action action = results.getAction(base);
                    if ( action == CONFIRM )
                        mCommon.confirm(readId,readIdx);
                    else if ( action == CORRECT )
                        mCommon.correct(readId,readIdx,cWinner);
                }
                else
                {
                    Action action = results.getAction(read[readIdx]);
                    if ( action == CONFIRM )
                        mCommon.confirm(readId,readIdx);
                    else if ( action == CORRECT )
                        mCommon.correct(readId,readIdx,winner);
                }

            }
        }
    }
    return mCommon.noBranches() && results.isBranchPoint();
}

typedef std::less<KMer<K> > MergeComp;
typedef ExternalSorter<ReadKmer>::Merger< VecReadKmer, MergeComp> Merger;

void processStacks( unsigned nThreads, unsigned coverage,
                    bool doBranches, bool skipUnder, unsigned cycle,
                    StackProc::Scorer const& scorer, vecbvec& reads,
                    VecQualNibbleVec& quals, Merger& merger )
{
    if ( !merger.size() )
    {
        ULOG("There are no stacks to process");
        return; // EARLY RETURN!
    }

    DLOG("Processing read stacks.");
    StackProc::Common common(cycle,K,doBranches,skipUnder,reads,quals);
    if ( true )
    {
        StackProc proc(common,scorer);
        Worklist<StackProc::Workitem,StackProc> wl(proc,nThreads);

        ReadKmer next;
        ReadLoc const& readLoc = next.getReadLoc();
        merger.getNext(next);
        KMer<K> last(next);

        StackProc::Workitem wi;
        size_t const VECS_PER_WORKITEM = 10000;
        wi.reserve(VECS_PER_WORKITEM);

        StackProc::ReadStack stack;
        stack.reserve(10*coverage);
        stack.push_back(readLoc);

        size_t depth = 0;
        size_t maxDepth = 2*(nThreads-1);
        size_t line = 0;
        while ( merger.getNext(next) )
        {
            if ( last != next )
            {
                if ( stack.size() >= scorer.getMinReadStackDepth() )
                {
                    wi.push_back(stack);
                    if ( wi.full() )
                    {
                        if ( depth < maxDepth )
                        {
                            depth = wl.add(wi);
                            std::cout << '.' << std::flush;
                        }
                        else
                        {
                            proc(wi);
                            depth = 0;
                            std::cout << "*" << std::flush;
                        }
                        if ( ++line == 80 )
                        { std::cout << std::endl; line = 0; }
                        wi.clear();
                    }
                }
                stack.clear();
                last = next;
            }
            stack.push_back(readLoc);
        }
        if ( stack.size() >= scorer.getMinReadStackDepth() )
            wi.push_back(stack);
        stack.clear();
        if ( wi.size() )
            wl.add(wi);
        std::cout << std::endl;
        wi.clear();
    }

    common.applyActions();
}

void findKeepers( unsigned tooFewKmers, Merger& merger, BitVec* pKeepers )
{
    if ( !merger.size() )
    {
        ULOG("There are no read kmers to process");
        return; // EARLY RETURN!
    }

    DLOG("Marking reads with low-frequency kmers.");
    std::vector<size_t> rlv;
    rlv.reserve(tooFewKmers+1);

    ReadKmer next;
    ReadLoc const& readLoc = next.getReadLoc();
    merger.getNext(next);
    rlv.push_back(readLoc.getID());
    KMer<K> last(next);

    typedef std::vector<size_t>::iterator Itr;
    while ( merger.getNext(next) )
    {
        if ( last != next )
        {
            if ( rlv.size() <= tooFewKmers )
                for ( Itr itr(rlv.begin()), end(rlv.end()); itr != end; ++itr )
                    pKeepers->set(*itr,false);

            rlv.clear();
            last = next;
        }
        if ( rlv.size() != rlv.capacity() )
            rlv.push_back(readLoc.getID());
    }
    if ( rlv.size() <= tooFewKmers )
        for ( Itr itr(rlv.begin()), end(rlv.end()); itr != end; ++itr )
            pKeepers->set(*itr,false);

    size_t nBad = pKeepers->size() - pKeepers->Sum();
    double pctBad = 100.*nBad/pKeepers->size();
    ULOG("There are " << nBad << " (" << pctBad << "%) uncorrectable reads.");
    PerfStat ps("frac_reads_removed",
                "% of reads removed because of low frequency kmers", pctBad);
    PerfStat::log() << std::fixed << std::setprecision(1) << ps;
}

void writeCorrectedFiles( vecbvec const& reads,
                          QualNibbleVecVec const& quals,
                          BitVec const& keepers,
                          String const& pairsFile,
                          String const& editHead )
{
    size_t nReads = reads.size();
    vec<Bool> losers;
    losers.reserve(nReads);
    IncrementalWriter<bvec> readWriter((editHead+".fastb").c_str(),nReads);
    IncrementalWriter<qvec> qualWriter((editHead+".qualb").c_str(),nReads);
    typedef BitVec::const_iterator Itr;
    for ( Itr itr(keepers.begin()), end(keepers.end()); itr != end; ++itr )
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

    PairsManager pairs(pairsFile);
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

} // end of anonymous namespace

TRIVIALLY_SERIALIZABLE(ReadLoc);
TRIVIALLY_SERIALIZABLE(ReadKmer);
extern template class SmallVec<ReadLoc, MempoolAllocator<ReadLoc> >;
extern template class OuterVec< SmallVec< ReadLoc, MempoolAllocator<ReadLoc> > >;


int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;

  // File locations:

  CommandArgument_String_OrDefault(ROOT, "");
  CommandArgument_String_OrDefault_Doc(IN_HEAD, "frag_reads_prec",
     "looks for ROOT/IN_HEAD.{fastb,qualb}");
  CommandArgument_String_OrDefault_Doc(OUT_EDIT_HEAD, "frag_reads_edit",
     "generates ROOT/OUT_EDIT_HEAD.{fastb,keep,log}");
  CommandArgument_String_OrDefault_Doc(OUT_CORR_HEAD, "frag_reads_corr",
     "generates ROOT/OUT_CORR_HEAD.{fastb,keep,log}");
  CommandArgument_String_OrDefault_Doc(WORKDIR,".",
      "Place to write scratch files.");


  // Computational performance:

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
     "Number of parallel threads.");
  CommandArgument_LongLong_OrDefault_Doc(MAX_MEMORY_GB, 0,
     "Maximum amount of memory to use in Gb");

  // Heuristics:

  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_READSTACK_DEPTH, 5,
     "Minimum number of same-kmer containing reads to consider building a "
     "ReadStack");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_BASESTACK_DEPTH, 1,
     "Minimum number of bases in a column to make inferences from it; if zero, "
     "set to same value as MIN_READSTACK_DEPTH");
  CommandArgument_UnsignedInt_OrDefault_Doc(AUTO_CONFIRM, 0,
     "if set, automatically confirm bases whose quality score sum is at least "
     "this large");
  CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CORRECT, 0.25,
     "Maximum ratio of QSS (Quality score sum) of existing to best base to "
     "mark it for correction");
  CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CORRECT2, 0.25,
     "Maximum ratio of QSS (Quality score sum) of second best to best base to "
     "mark any base for correction");
  CommandArgument_Double_OrDefault_Doc(MAX_QSS_RATIO_TO_CONFIRM, 1.0,
     "Maximum ratio of QSS (Quality score sum) of second best to best base to "
     "confirm any base");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_Q_TO_SUPPORT, 20,
     "Minimum base quality score (Q) to mark it a 'supporting' or 'confirming' "
     "base");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_N_BASES_TO_SUPPORT, 2,
     "Minimum number supporting bases of a nucleotide to mark it 'supported'");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_QSS_TO_SUPPORT, 60,
     "If positive, minimum quality score sum of supporting bases of a "
     "nucleotide to mark it 'supported'; overrides MIN_N_BASES_TO_SUPPORT and "
     "MIN_Q_TO_SUPPORT");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_N_BASES_TO_CONFIRM, 3,
     "Minimum number supporting bases of a nucleotide to mark it 'confirmed'");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_QSS_TO_CONFIRM, 90,
     "Minimum QSS (Quality score sum) of bases of a nucleotide to mark it "
     "'base_locked'; overrides MIN_N_BASES_TO_CONFIRM & MIN_Q_TO_SUPPORT");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_QSS_TO_WIN, 40,
     "Ignore columns for which the quality score sum for the winning base is "
     "not at least this value.");

  CommandArgument_UnsignedInt_OrDefault_Doc(QUAL_CEIL_RADIUS, 2,
     "Replace each quality score by the minimum of the quality scores over a "
     "window of this radius.");
  CommandArgument_Bool_OrDefault_Doc(QCR_ALWAYS, True,
     "Apply QUAL_CEIL_RADIUS every cycle. Doesn't really make sense to do it "
     "more than once but may improve results.");

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_CYCLES, 2,
     "Number of correction cycles.");

  CommandArgument_UnsignedInt_OrDefault_Doc(MAX_KMER_FREQ_TO_MARK, 2,
     "Reads with minimum kmer frequencies at and below this will be marked and "
     "discarded in _corr.");

  CommandArgument_Bool_OrDefault_Doc(DO_PALINDROMES, False,
     "Includes read stacks built from palidromic kmers.");
  CommandArgument_Bool_OrDefault_Doc(DO_BRANCHES, False,
     "Continues error correction after finding a suspected branching.");
  CommandArgument_Bool_OrDefault_Doc(SKIP_UNDER, True,
     "Don't do anything with bases under kmers");

  CommandArgument_String_OrDefault_Doc(USE_KMER_SPECTRUM, "",
     "If set to name of a spectrum_summary file created by KmerSpectrum, "
     "try to automatically adjust hueristic parameters based on spectrum "
     "features.  Modifies MIN_READSTACK_DEPTH, MIN_BASESTACK_DEPTH, "
     "MAX_KMER_FREQUENCY_TO_MARK, MIN_QSS_TO_WIN, MIN_QSS_TO_CONFIRM and "
     "MIN_QSS_TO_SUPPORT by scaling factors based on kmer spectrum.");
  CommandArgument_UnsignedInt_OrDefault_Doc(ESTIMATED_COVERAGE, 0,
     "If positive, use this as estimated coverage for scaling instead of the "
     "estimate from the kmer spectrum.");

  CommandArgument_Bool_OrDefault_Doc(DELETE, False,
     "actually delete reads, rather than mark them for deletion");

  EndCommandArguments;

  // Threads and memory limits

  SetMaxMemory(MAX_MEMORY_GB << 30);
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Parameter fiddling

  if (MIN_BASESTACK_DEPTH == 0)
    MIN_BASESTACK_DEPTH = MIN_READSTACK_DEPTH;
  ForceAssertLe(MIN_BASESTACK_DEPTH, MIN_READSTACK_DEPTH);

  if (USE_KMER_SPECTRUM != "" || ESTIMATED_COVERAGE > 0)
    useKmerSpectrum(USE_KMER_SPECTRUM,
                    ESTIMATED_COVERAGE,
                    MIN_READSTACK_DEPTH,
                    MIN_BASESTACK_DEPTH,
                    MIN_QSS_TO_CONFIRM,
                    MIN_QSS_TO_SUPPORT,
                    MIN_QSS_TO_WIN,
                    MAX_KMER_FREQ_TO_MARK);

  if ( !ESTIMATED_COVERAGE )
      ESTIMATED_COVERAGE = 75;

  StackProc::Scorer scorer(
          MIN_Q_TO_SUPPORT, MIN_N_BASES_TO_SUPPORT, MAX_QSS_RATIO_TO_CORRECT,
          MAX_QSS_RATIO_TO_CORRECT2, MAX_QSS_RATIO_TO_CONFIRM,
          MIN_QSS_TO_CONFIRM, MIN_QSS_TO_WIN, MIN_READSTACK_DEPTH,
          MIN_BASESTACK_DEPTH, AUTO_CONFIRM, MIN_QSS_TO_SUPPORT,
          MIN_N_BASES_TO_CONFIRM);

  // Path twiddling

  if ( ROOT.size() )
  {
      if ( ROOT.back() != '/' )
          ROOT.append(1,'/');
      IN_HEAD = ROOT + IN_HEAD;
      OUT_EDIT_HEAD = ROOT + OUT_EDIT_HEAD;
      OUT_CORR_HEAD = ROOT + OUT_CORR_HEAD;
  }

  DLOG("Reading reads.");
  vecbvec reads(IN_HEAD+".fastb");

  DLOG("Reading quals.");
  VecQualNibbleVec quals;
  LoadQualNibbleVec(IN_HEAD+".qualb",&quals);

  ExternalSorter<ReadKmer> xs(WORKDIR+"/FE.ReadKmers.tmp.");
  VecReadKmer scratch(K);

  for ( unsigned cycle = 0; cycle < NUM_CYCLES; ++cycle )
  {
      ULOG("");
      DLOG("Starting cycle " << cycle+1);

      if ( QUAL_CEIL_RADIUS )
          squashQuals(quals,NUM_THREADS,QUAL_CEIL_RADIUS);
      if ( !QCR_ALWAYS )
          QUAL_CEIL_RADIUS = 0;

      kmerizeReads(reads,NUM_THREADS,DO_PALINDROMES,xs);
      Merger merger(xs,scratch);
      processStacks(NUM_THREADS,ESTIMATED_COVERAGE,DO_BRANCHES,SKIP_UNDER,
                      cycle,scorer,reads,quals,merger);
  }

  ULOG("");
  DLOG("Finding keepers.");
  BitVec keepers;
  keepers.resize(reads.size(),true);
  if ( true )
  {
      kmerizeReads(reads,NUM_THREADS,true,xs);
      Merger merger(xs,scratch);
      findKeepers(MAX_KMER_FREQ_TO_MARK,merger,&keepers);
  }

  ULOG("");
  DLOG("Writing edited output files.");
  reads.WriteAll(OUT_EDIT_HEAD+".fastb");
  WriteAll(quals,OUT_EDIT_HEAD+".qualb");
  String keepFile(OUT_EDIT_HEAD+".keep");
  BinaryWriter::writeFile(keepFile.c_str(),keepers);

  if ( DELETE )
  {
      DLOG("Writing corrected output files.");
      writeCorrectedFiles(reads,quals,keepers,IN_HEAD+".pairs",OUT_CORR_HEAD);
  }

  DLOG("Done.");
}

#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

template class SmallVec<ReadLoc, MempoolAllocator<ReadLoc> >;
template class OuterVec< SmallVec< ReadLoc, MempoolAllocator<ReadLoc> > >;
