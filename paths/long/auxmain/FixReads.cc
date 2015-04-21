///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FixReads.cc
 * \author tsharpe
 * \date Sep 25, 2012
 *
 * \brief Given reads and quals, pull read stack from high-x kmers flanking
 * low-x region, and substitute the stack's consensus for the low-x region.
 */

#include "MainTools.h"
#include "Basevector.h"
#include "ParallelVecUtilities.h"
#include "ReadError.h"
#include "kmers/ReadPatherDefs.h"
#include "math/Combinatorics.h"
#include "pairwise_aligners/EditDistance.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/Logging.h"
#include "paths/long/ultra/MultipleAligner.h"
#include "simulation/BamReadGenerator.h"
#include "simulation/ReadSimulatorSimpleGenerator.h"
#include "simulation/ReadErrorAnalyzer.h"
#include "system/WorklistN.h"
#include "paths/HyperBasevector.h"
#include "paths/long/RefTrace.h"
#include <algorithm>

namespace
{

class RLoc
{
public:
    static unsigned const MAX_READ_LEN = 255;

    RLoc() : mOffset(0), mRC(false) {}
    RLoc( ReadID const& readID, unsigned offset, bool rc = false )
    : mReadID(readID), mOffset(offset), mRC(rc)
    { ForceAssertLe(offset,MAX_READ_LEN); }

    // compiler-supplied copying and destructor are OK

    ReadID const& getReadID() const { return mReadID; }
    unsigned getOffset() const { return mOffset; }
    bool isRC() const { return mRC; }

    friend bool operator<( RLoc const& rloc1, RLoc const& rloc2 )
    { bool result = false;
      int cmp = compare(rloc1.mReadID,rloc2.mReadID);
      if ( cmp < 0 ) result = true;
      else if ( !cmp ) result = rloc1.mOffset < rloc2.mOffset;
      return result; }

    friend ostream& operator<<( ostream& os, RLoc const& rloc )
    { os << "Read=" << rloc.mReadID << " Off=" << (unsigned)rloc.mOffset
        << " RC=" << (rloc.mRC?'Y':'N');
      return os; }

private:
    ReadID mReadID;
    unsigned char mOffset;
    bool mRC;
};
unsigned const RLoc::MAX_READ_LEN;
typedef vec<RLoc> RLocVec;

class KLoc
{
public:
    KLoc() : mLen(0) {}
    KLoc( KmerID const& kmerID, unsigned len )
    : mKmerID(kmerID), mLen(len) {}
    KLoc( KmerID const& kmerID, size_t readId, unsigned offset, bool rc )
    : mKmerID(kmerID), mLen(1), mRLoc(ReadID(readId),offset,rc) {}

    // compiler-supplied copying and destructor are OK

    bool inBounds( size_t kmerId )
    { size_t val = mKmerID.val();
      return val <= kmerId && kmerId < val+mLen; }

    KmerID const& getKmerID() const { return mKmerID; }
    unsigned getLen() const { return mLen; }
    RLoc const& getRLoc() const { return mRLoc; }
    RLoc getRLoc( size_t kmerId )
    { unsigned diff = kmerId - mKmerID.val();
      unsigned off = mRLoc.isRC() ? mRLoc.getOffset()+mLen-1-diff :
                                    mRLoc.getOffset()+diff;
      return RLoc(mRLoc.getReadID(),off,mRLoc.isRC()); }

    bool extend( KLoc const& kLoc )
    { bool result = false;
      if ( mRLoc.getReadID() == kLoc.mRLoc.getReadID() &&
              mRLoc.getOffset()+1 == kLoc.mRLoc.getOffset() &&
              mRLoc.isRC() == kLoc.mRLoc.isRC() &&
              mKmerID.val() == kLoc.mKmerID.val()+(mRLoc.isRC()?1:-1) )
      { ForceAssertLt(mLen,RLoc::MAX_READ_LEN); mLen += 1; result = true;
        if ( mRLoc.isRC() ) mKmerID = kLoc.mKmerID; }
      return result; }

    friend bool operator<( KLoc const& kloc1, KLoc const& kloc2 )
    { bool result = false;
      int cmp = compare(kloc1.mKmerID,kloc2.mKmerID);
      if ( cmp < 0 ) result = true;
      else if ( !cmp ) result = kloc1.mLen < kloc2.mLen;
      return result; }

private:
    KmerID mKmerID;
    unsigned char mLen;
    RLoc mRLoc;
};
typedef vec<KLoc> KLocVec;

template <unsigned K>
class ReadFixer
{
public:
    ReadFixer( Generator const& generator, unsigned nThreads, double errRate )
    : mReads(generator.getReads()), mpDict(buildDict(nThreads)), mGraph(*mpDict),
      mErrRate(errRate), mpCorrectedReads(0)
    { buildKLocVec(nThreads); mErrors.resize(mReads.size()); }

    ~ReadFixer() { delete mpDict; }

    // finds errors, and corrects reads
    void correctReads( unsigned nThreads, vecbvec* pVBV )
    { pVBV->assign(mReads.begin(),mReads.end());
      mpCorrectedReads = pVBV;
      findErrors(nThreads); }

    // just finds errors.  call getErrs to learn what we found.
    void findErrors( unsigned nThreads )
    { size_t nReads = mReads.size();
      Dotter dotter(nReads);
      FixerProc proc(this,&dotter);
      parallelFor(0ul,nReads,proc,nThreads); }

    void findReadErrors( size_t readId );

    ReadErrorVecVec const& getErrs() const { return mErrors; }

private:
    ReadFixer( ReadFixer const& ); // copying prohibited
    ReadFixer& operator=( ReadFixer const& ); // copying prohibited
    static size_t const BATCH_SIZE = 1000;

    KmerDict<K>* buildDict( unsigned nThreads );
    void buildKLocVec( unsigned nThreads );
    void fixReadSegment( size_t readId, unsigned start, unsigned end, bool rc,
                            size_t kStart, size_t kEnd, ReadErrorVec& errs );

    static bool sepsNearlyEqual( int sep1, int sep2 )
    { return static_cast<unsigned>(std::abs(sep1-sep2)) <= SEP_TOLERANCE; }

    // TODO: hard-coded value isn't good
    static unsigned const SEP_TOLERANCE = 2;

    vecbvec const& mReads;
    KmerDict<K>* mpDict;
    UnipathGraph<K> mGraph;
    double mErrRate;
    vecbvec* mpCorrectedReads;

    // TODO: ought to use a more sophisticated data structure like an
    // Interval Tree instead of a vector
    KLocVec mKLocVec;
    Scorer mScorer; // edit-distance scorer
    ReadErrorVecVec mErrors;

    // KmerEater that emits KLoc's.
    class KLocKmerEater
    {
    public:
        KLocKmerEater( UnipathGraph<K>* pGraph, KLocVec* pKLocVec )
        : mpGraph(pGraph), mpKLocVec(pKLocVec) {}

        ~KLocKmerEater()
        { if ( mLastLoc.getLen() ) mKLocVec.push_back(mLastLoc);
          if ( !mKLocVec.empty() )
          { static SpinLockedData gLock;
            SpinLocker locker(gLock);
            mpKLocVec->insert(mpKLocVec->end(),
                                mKLocVec.begin(),mKLocVec.end()); } }

        // compiler-supplied copying is OK

        void operator()(KMer<K> const& kmer, KMerContext,
                                    size_t rId, size_t rOff);

    private:
        UnipathGraph<K>* mpGraph;
        KLocVec mKLocVec;
        KLocVec* mpKLocVec;
        KLoc mLastLoc;
    };

    // worklist proc for parallelizing read fixing
    class FixerProc
    {
    public:
        FixerProc( ReadFixer* pFixer, Dotter* pDotter )
        : mpFixer(pFixer), mpDotter(pDotter) {}

        // compiler-supplied copying and destructor are OK

        void operator()( size_t readId )
        { mpFixer->findReadErrors(readId); mpDotter->batchDone(); }

    private:
        ReadFixer* mpFixer;
        Dotter* mpDotter;
    };

};

// class to apply consensus findings to patch a read
class Patcher
{
public:
    typedef MultipleAligner::ConsScore Score;

    Patcher( double errRate )
    : mErrRate(errRate) {}

    // compiler-supplied copy constructor and destructor are OK
    // assignment prohibited by reference member

    void patch( MultipleAligner& ma, bvec const& seg, unsigned off,
                    ReadErrorVec& errs )
    { gvec alignment(ma.alignRead(seg.begin(),seg.end()));
      typedef MultipleAligner::VecConsScore VCS;
      typedef VCS::const_iterator Itr;
      VCS const& vcs = ma.getScores();
      ForceAssertEq(alignment.size(),vcs.size());
      gvec::const_iterator rdItr(alignment.begin());
      Itr end(vcs.end());
      for ( Itr itr(vcs.begin()); itr != end; ++itr, ++rdItr, ++off )
      { Score const& score = *itr;
        gvec::value_type consBase = score.bestVal();
        gvec::value_type readBase = *rdItr;
        if ( readBase != consBase && isCorrectable(score) )
        { ReadError::ErrType type = ReadError::SUBSTITUTION;
          if ( consBase == gvec::GAP )
            type = ReadError::INSERTION;
          else if ( readBase == gvec::GAP )
            type = ReadError::DELETION;
          errs.push_back(ReadError(off,type,readBase,consBase)); } } }

private:
    bool isCorrectable( Score const& score )
    { double p = 1. - BinomialProbCumFailure(mErrRate,
                                        score.nVals()-1,
                                        score.getErrCount());
      return p > mErrRate; }

    double mErrRate;
};

template <unsigned K>
KmerDict<K>* ReadFixer<K>::buildDict( unsigned nThreads )
{
    size_t dictAllocSize = mReads.SizeSum()/20;
    KmerDict<K>* pDict = new KmerDict<K>(dictAllocSize);
    pDict->process(mReads,true,false,nThreads,BATCH_SIZE);
    std::cout << "Dictionary:  allocated size=" << dictAllocSize
            << " actual size=" << pDict->size() << std::endl;

    Spectrum<K> spectrum(*pDict);
    if ( !spectrum.isFirstPeakValid() )
        FatalErr("Couldn't parse kmer spectrum.");

    size_t goodIdx = spectrum.getFirstTroughIndex()+1;
    size_t dictSize = spectrum.sum(goodIdx);
    pDict->clean(typename KmerDict<K>::BadKmerCountFunctor(goodIdx));
    AssertEq(dictSize,pDict->size());
    return pDict;
}

template <unsigned K>
void ReadFixer<K>::buildKLocVec( unsigned nThreads )
{
    size_t nReads = mReads.size();
    size_t vecAllocSize = 16*nReads;
    mKLocVec.reserve(vecAllocSize);
    typedef KLocKmerEater Eater;
    typedef vecbvec::const_iterator Itr;
    typedef KmerizationProcessor<K,Eater,Itr> KProc;
    size_t nBatches = (nReads+BATCH_SIZE-1)/BATCH_SIZE;
    Dotter dotter(nBatches);
    Eater eater(&mGraph,&mKLocVec);
    KProc proc(mReads.begin(),mReads.end(),eater,&dotter);
    std::cout << Date() << ": Creating KLocVec." << std::endl;
    parallelFor(0ul,dotter.getNBatches(),proc,nThreads);
    std::cout << "KLocVec:  allocated size=" << vecAllocSize
            << " actual size=" << mKLocVec.size() << std::endl;

    std::cout << Date() << ": Sorting KLocVec." << std::endl;
    ParallelSort(mKLocVec);

    std::cout << Date() << ": KLocVec sorted." << std::endl;
}

template <unsigned K>
void ReadFixer<K>::findReadErrors( size_t readId )
{
    bvec const& read = mReads[readId];
    if ( read.size() <= K )
        return; // Non-structured!

    ReadErrorVec errs;
    KMer<K> kmer(read.begin());
    KmerDictEntry<K> const* pEntry = mpDict->findEntry(kmer);
    bvec::const_iterator rItr(read.begin(K));
    bvec::const_iterator rEnd(read.end());
    while ( rItr != rEnd )
    {
        kmer.toSuccessor(*rItr);
        ++rItr;
        KmerDictEntry<K> const* pEntry2 = mpDict->findEntry(kmer);
        if ( !pEntry || pEntry2 )
            pEntry = pEntry2;
        else if ( !pEntry2 )
        {
            KmerID start(mGraph.getKmerID(pEntry->getKDef()));
            int startOffset = rItr.pos()-K-1;
            while ( rItr != rEnd )
            {
                kmer.toSuccessor(*rItr);
                ++rItr;
                if ( (pEntry2 = mpDict->findEntry(kmer)) )
                {
                    KmerID end(mGraph.getKmerID(pEntry2->getKDef()));
                    int endOffset = rItr.pos()-K;
                    bool rc = KMer<K>(mGraph.getBases(end))!=kmer;
                    fixReadSegment(readId,startOffset,endOffset,rc,
                                       start.val(),end.val(),errs);
                    pEntry = pEntry2;
                    break;
                }
            }
        }
    }
    mErrors[readId] = errs;
    if ( mpCorrectedReads )
    {
        bvec& corrRd = mpCorrectedReads[0][readId];
        bvec tmp;
        typedef ReadErrorVec::const_reverse_iterator Itr;
        for ( Itr itr(errs.crbegin()), end(errs.crend()); itr != end; ++itr )
        {
            ReadError err = *itr;
            unsigned loc = err.getLocation();
            switch ( err.getType() )
            {
            case ReadError::DELETION:
                tmp.assign(corrRd.begin(loc),corrRd.end());
                corrRd.resize(loc).push_back(err.getRefBase()).append(tmp);
                break;
            case ReadError::INSERTION:
                tmp.assign(corrRd.begin(loc+1),corrRd.end());
                corrRd.resize(loc).append(tmp);
                break;
            case ReadError::SUBSTITUTION:
                corrRd.set(loc,err.getRefBase());
                break;
            }
        }
    }
}

template <unsigned K>
void ReadFixer<K>::fixReadSegment( size_t readId, unsigned startOffset,
                                    unsigned endOffset, bool rc,
                                    size_t kStart, size_t kEnd,
                                    ReadErrorVec& errs )
{
    using std::lower_bound;
    using std::min;
    typedef KLocVec::iterator Itr;

    // figure out what we're fixing and how far apart the bracketing kmers are
    unsigned chopLen = min((endOffset-startOffset+K-SEP_TOLERANCE-6)/2,K-3);
    bvec const& readToFix = mReads[readId];
    bvec fixSeg(readToFix.begin(startOffset+chopLen),
                    readToFix.begin(endOffset+K-chopLen));
    unsigned curEditDistance = fixSeg.size();
    int separation = endOffset-startOffset;
    if ( rc )
        separation = -separation;

    // find all the reads that contain kEnd
    Itr end(mKLocVec.end());
    Itr itr(lower_bound(mKLocVec.begin(),end,
                        KLoc(KmerID(kEnd-min(kEnd,RLoc::MAX_READ_LEN-1ul)),
                                RLoc::MAX_READ_LEN)));
    RLocVec rLocVec;
    rLocVec.reserve(200);
    for ( ; itr != end && itr->getKmerID().val()<=kEnd; ++itr )
        if ( itr->inBounds(kEnd) )
            rLocVec.push_back(itr->getRLoc(kEnd));

    std::sort(rLocVec.begin(),rLocVec.end());

    // iterate over all the reads that contain kStart
    vecbvec readSegs;
    readSegs.reserve(rLocVec.size());
    itr = Itr(lower_bound(mKLocVec.begin(),end,
                         KLoc(KmerID(kStart-min(kStart,RLoc::MAX_READ_LEN-1ul)),
                                 RLoc::MAX_READ_LEN)));
    for ( ; itr != end && itr->getKmerID().val()<=kStart; ++itr )
    {
        if ( itr->inBounds(kStart) )
        {
            RLoc rLoc(itr->getRLoc(kStart));
            size_t targetReadId = rLoc.getReadID().val();
            if ( targetReadId == readId )
                continue; // NOT STRUCTURED!

            typedef RLocVec::iterator RItr;
            RLoc startRLoc(rLoc.getReadID(),0);
            RItr rItr(lower_bound(rLocVec.begin(),rLocVec.end(),startRLoc));
            RLoc endRLoc(ReadID(targetReadId+1),0);
            RItr rEnd(lower_bound(rItr,rLocVec.end(),endRLoc));
            for ( ; rItr != rEnd; ++rItr )
            {
                int testSep = rItr->getOffset()-rLoc.getOffset();
                if ( rItr->isRC() )
                    testSep = -testSep;
                if ( sepsNearlyEqual(separation,testSep) )
                {
                    bvec const& targetRead = mReads[targetReadId];
                    bvec seg;
                    if ( rLoc.getOffset() > rItr->getOffset() )
                    {
                        unsigned rdLen = targetRead.size();
                        unsigned start = rdLen-rLoc.getOffset()-K+chopLen;
                        unsigned stop = rdLen-rItr->getOffset()-chopLen;
                        seg = bvec(targetRead.rcbegin(start),
                                    targetRead.rcbegin(stop));
                    }
                    else
                    {
                        unsigned start = rLoc.getOffset()+chopLen;
                        unsigned stop = rItr->getOffset()+K-chopLen;
                        seg = bvec(targetRead.begin(start),
                                    targetRead.begin(stop));
                    }
                    unsigned edDist = editDistance(fixSeg.begin(),fixSeg.end(),
                                                    seg.begin(),seg.end());
                    if ( edDist < curEditDistance )
                    {
                        readSegs.clear();
                        curEditDistance = edDist;
                    }
                    if ( edDist == curEditDistance )
                        readSegs.push_back(seg);
                }
            }
        }
    }
    MultipleAligner ma(mScorer,fixSeg.size());
    ma.addReads(readSegs.begin(),readSegs.end());
    Patcher patcher( mErrRate );
    patcher.patch(ma,fixSeg,startOffset+chopLen,errs);
}

template <unsigned K>
void ReadFixer<K>::KLocKmerEater::operator()( KMer<K> const& kmer, KMerContext,
                                                size_t rId, size_t rOff )
{
    if ( !mKLocVec.capacity() )
        mKLocVec.reserve(1000000);
    KmerDictEntry<K> const* pEntry = mpGraph->getDict().findEntry(kmer);
    if ( pEntry )
    {
        KDef const& kDef = pEntry->getKDef();
        KmerID kmerID = mpGraph->getKmerID(kDef);
        KLoc kLoc(kmerID, rId, rOff, kmer!=KMer<K>(mpGraph->getBases(kmerID)));
        if ( !mLastLoc.extend(kLoc) )
        {
            if ( mLastLoc.getLen() )
                mKLocVec.push_back(mLastLoc);
            mLastLoc = kLoc;
        }
    }
}

void buildRefReads( const vecbvec& ref, RefLocusVec const& refLocs,
                        vecbvec* pVBV )
{
    pVBV->reserve(refLocs.size());
    bvec tmp;
    typedef RefLocusVec::const_iterator Itr;
    for ( Itr itr(refLocs.begin()), end(refLocs.end()); itr != end; ++itr )
    {
        RefLocus const& loc = *itr;
        bvec const& refBV = ref[loc.mRefID];
        bvec::const_iterator beg(refBV.begin(loc.mOffset));
        int wrap = loc.mOffset + loc.mLength - refBV.size();
        if ( wrap <= 0 )
            tmp.assign(beg,beg+loc.mLength);
        else
            tmp.assign(beg,refBV.end()).append(refBV.begin(),refBV.begin(wrap));
        if ( loc.mRC )
            tmp.ReverseComplement();
        pVBV->push_back(tmp);
    }
}


}

int main( int argc, char *argv[] )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_Double_OrDefault_Doc(ERR_RATE, .001, "Likelihood of a sequencing error.");

    CommandArgument_Bool_OrDefault_Doc(USE_SIM_READS, False, "Whether to use simulated reads.");

    // === parameters for simulated reads ===
    CommandArgument_String_OrDefault_Doc(SAMPLE, "rhody", "which sample to use");
     CommandArgument_String_OrDefault_Doc(X, "", 
          "Define regions of genome to use.  If nothing specified use the entire "
          "genome.  A single nonnegative integer gets you a particular "
          "genome contig.  A range like 2:1000-20000 gets you e.g. "
          "bases 1000 to 20000 on reference contig 2.  Or you may provide a "
          "comma-separated list {...} whose elements are genome contigs ids or "
          "ranges as above.  A ParseIntSet list of contig ids may be used.  Finally "
          "one may use \"random(n)\" where n is a positive integer, which selects a "
          "clock-seeded random region of the genome of size n.");
    CommandArgument_Int_OrDefault_Doc(COVERAGE, 50, "Read coverage");
    CommandArgument_Int_OrDefault_Doc(LEN, 250, "Read length mean");
    CommandArgument_Int_OrDefault_Doc(LEN_SIG, 0, "Read length standard deviation");
    CommandArgument_Int_OrDefault_Doc(LEN_MIN, 0, "Read length minimum");
    CommandArgument_Double_OrDefault_Doc(ERR_DEL, 0.0, "Deletion error rate.");
    CommandArgument_Double_OrDefault_Doc(ERR_SUB, 0.001, "Substitution error rate.");
    CommandArgument_Double_OrDefault_Doc(ERR_INS, 0.0, "Insertion error rate.");
    CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, 123456, "Random seeds");
    CommandArgument_String_OrDefault_Doc(RIDS, "", "Which reads to check");

    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0u, "Number of parallel threads.");
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);

    auto_ptr<Generator> reads_generator;
    vecbasevector G;

    if ( ! USE_SIM_READS ) {
        std::vector<CachedBAMFile> bams;
        bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit2"));
        bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit3"));
        bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit4"));
        bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit5"));
        std::cout << Date() << " Reading BAMs." << std::endl;
        auto_ptr<Generator> brg( new BamReadGenerator(bams) );

        G.ReadAll( "/wga/dev/references/Rhodobacter_sphaeroides/genome.fastb" );
        reads_generator =  brg;
    } else {
         // Create genome.
         double GENOME_SUB_PERCENT = 0.0;
         //vecbasevector G, Gplus, G3, G3plus;
         vecbasevector Gplus, G3, G3plus;
         vec<int> Gplus_ext;
         vec<HyperBasevector> GH;
         vec<double> ploidy;
         vec<bool> is_circular;
         VecIntPairVec Glocs, G3locs, G3pluslocs;
         const int LG = 12;
         String X_actual;
         ref_data ref;
         long_logging logc( "" );
         if ( SAMPLE != "unknown" )
         {
             RefTraceControl RTCtrl;
             BuildRefData( "", SAMPLE, X, "", X_actual, ref, GENOME_SUB_PERCENT, logc, RTCtrl );
         } else {
             cout << "SAMPLE=unknown is not allowed." << endl;
             return 1;
         }
         // Simulated the reads
        auto_ptr<ReadSimulatorSimpleGenerator>  rssg( new ReadSimulatorSimpleGenerator(G, is_circular) );
        RandomLength insert_len( LEN, LEN_SIG, LEN_MIN, RANDOM_SEED );
        rssg->computeReadsSizeAndPositionUnpaired (0, COVERAGE, RANDOM_SEED, false, 0.0, insert_len );
        rssg->buildReadsFromPositionsWithErrors(RANDOM_SEED,ERR_DEL, ERR_INS, ERR_SUB );

        reads_generator = rssg;
    }

    // === Fix the reads ===
    std::cout << Date() << " Pathing reads." << std::endl;
    unsigned const K = 25;
    ReadFixer<K> fixer(*reads_generator, NUM_THREADS, ERR_RATE);
    std::cout << Date() << " Correcting reads." << std::endl;
    vecbvec corr;
    fixer.correctReads(NUM_THREADS,&corr);

    // === Evaluation ===
    std::cout << Date() << " Evaluating corrections." << std::endl;
    vecbvec refReads;
    buildRefReads( G, reads_generator->getReadLocs(), &refReads );
    // Analyze overall performance
    ReadErrorAnalyzer rea(refReads, reads_generator->getReads(), corr);
    rea.SetVerbosity(0);
    rea.Analyze();

    // Analyze all the reads and print any reads with missed or added errors
    // rea.AnalyzeAll();
    
    // Check individual reads
    if ( RIDS != "" ) {
        FriendTeller fteller( reads_generator->getReadLocs(), G );
        rea.AddFriendTeller( &fteller );
            rea.SetVerbosity(2);
            vec<int> rids;
            ParseIntSet( RIDS, rids );
            for ( size_t i = 0; i < rids.size(); ++i ) 
                rea.AnalyzeRead( rids[i] );
    }
    std::cout << Date() << ": Done." << std::endl;
}
