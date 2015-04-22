///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FillFragmentsNew.cc
 * \author tsharpe
 * \date Aug 10, 2011
 *
 * \brief
 */
#ifndef NDEBUG
#define LOG(x) std::cout << x << std::endl
#else
#define LOG(x)
#endif

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "String.h"
#include "Vec.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/ReadPather.h"
#include "reporting/PerfStat.h"
#include <algorithm>
#include <fstream>
#include <math.h>

namespace
{

bvec::size_type getMaximumReadLen( String const& fastbFilename )
{
    std::vector<bvec::size_type> lengths = bvec::getSizes(fastbFilename);
    bvec::size_type result = 0;
    if ( lengths.size() )
        result = *std::max_element(lengths.begin(),lengths.end());
    return result;
}

class LibSepStats
{
public:
    LibSepStats() : mNPairs(0), mN(0), mSum(0), mSumSq(0) {}

    // compiler-supplied copying and destructor are OK

    void addSep( int sep ) { mN += 1; mSum += sep; mSumSq += sep*sep; }

    void addPair() { mNPairs += 1; }

    size_t getNPairs() const { return mNPairs; }
    size_t getN() const { return mN; }
    double getMean() const { return mN ? mSum/mN : 0.; }
    double getStdDev() const
    { return mN>1 ? sqrt((mSumSq-mSum*mSum/mN)/(mN-1)) : 0.; }

private:
    size_t mNPairs;
    size_t mN;
    double mSum;
    double mSumSq;
};

template <unsigned K>
class FragmentFiller
{
public:
    FragmentFiller( String const& pairsFile,
                    String const& pathsFile, String const& graphFile,
                    long maxSep, unsigned minOverlap, double minFillFrac )
    : mPM(pairsFile), mPWE(pathsFile,graphFile),
      mMaxSep(maxSep), mMinOverlap(minOverlap), mMinFillFrac(minFillFrac)
    { mSeps.resize(mPM.nLibraries()); mClosures.resize(mPM.nPairs()); }

    // compiler-supplied destructor is OK

    size_t getNFills() const
    { typedef std::vector<LibSepStats>::const_iterator Itr;
      size_t result = 0;
      for ( Itr itr(mSeps.begin()), end(mSeps.end()); itr != end; ++itr )
        result += itr->getN();
      return result; }

    size_t getNPairs() const { return mPM.nPairs(); }

    void processPairs( vec<longlong> const& ids )
    { typedef vec<longlong>::const_iterator Itr;
      for ( Itr itr(ids.begin()), end(ids.end()); itr != end; ++itr )
        processPair(*itr); }

    void processAllPairs()
    { size_t nPairs = mPM.nPairs();
      for ( size_t idx = 0; idx != nPairs; ++idx )
        processPair(idx); }

    void processPair( size_t id );

    PairsManager const& updateLibStats();

    void writeClosures( IncrementalWriter<bvec>& writer ) const
    { typedef std::vector<EdgeList>::const_iterator Itr;
      for ( Itr itr(mClosures.begin()), end(mClosures.end()); itr!=end; ++itr )
        if ( itr->size() )
          writer.add(mPWE.getBases(*itr)); }

private:
    FragmentFiller( FragmentFiller const& ); // unimplemented -- no copying
    FragmentFiller& operator=( FragmentFiller const& ); // unimplemented -- no copying

    bool exploreGraph( EdgeList const& el1, EdgeList const& el2,
                        EdgeDesc const& cur, EdgeList& fill, size_t& nTries,
                        long sep, EdgeList& closure );

    size_t getMaxSegmentsBetween( EdgeDesc const& ed1, EdgeDesc const& ed2 );

    PairsManager mPM;
    PathsWithEvidence<K> mPWE;
    std::vector<LibSepStats> mSeps;
    long mMaxSep; // maximum fragment separation in kmers
    unsigned mMinOverlap; // minimum number of overlapping kmers in filler
    double mMinFillFrac; // minimum fraction of fills before updating lib stats
    std::vector<EdgeList> mClosures;
};

template <unsigned K>
PairsManager const& FragmentFiller<K>::updateLibStats()
{
    size_t nnn = mPM.nLibraries();
    double maxSep = mMaxSep + K - 1;
    for ( size_t idx = 0; idx < nnn; ++idx )
    {
        LibSepStats const& lss = mSeps[idx];
        if ( lss.getNPairs() )
        {
            double fillFrac = 1.*lss.getN()/lss.getNPairs();
            String libName(mPM.getLibraryName(idx));
            std::cout << "Library " << libName << " had "
                      << std::setprecision(3) << 100.*fillFrac << "% of its "
                      << lss.getNPairs() << " pairs filled." << std::endl;
            if ( fillFrac < mMinFillFrac )
                std::cout << "Library statistics unchanged for library "
                          << libName << ": Too few pairs filled." << std::endl;
            else if ( lss.getMean()+3*lss.getStdDev() > maxSep )
                std::cout << "Library statistics unchanged for library "
                          << libName << ": Mean separation is large compared "
                          "to the read lengths." << std::endl;
            else
            {
                int mean = round(lss.getMean());
                int sd = round(lss.getStdDev());
                std::cout << "Library " << libName << " statistics changed to "
                          << mean << "(+-" << sd << ") from "
                          << mPM.getLibrarySep(idx) << "(+-"
                          << mPM.getLibrarySD(idx) << ')' << std::endl;
                mPM.changeLibrarySepSd(idx,mean,sd);
            }
        }
    }
    return mPM;
}

template <unsigned K>
void FragmentFiller<K>::processPair( size_t pairId )
{
    typedef EdgeList::const_iterator ELItr;

    if ( pairId > mPM.nPairs() )
    {
        std::cout << "Ignoring invalid pair id " << pairId << std::endl;
        return; // EARLY RETURN!
    }

    EdgeList el1(mPWE.getEdgeList(mPM.ID1(pairId)));
    EdgeList el2(mPWE.getEdgeList(mPM.ID2(pairId)));
    el2.rc(); // switch 2nd path around so the reads are pointing the same way

    LOG("Working on " << pairId << ": " << el1 << " + " << el2);

    using std::min;
    size_t nnn = min(el1.size(),el2.size());
    if ( !nnn )
    {
        LOG("Giving up:  One or both of the ends are unpathed.");
        return; // EARLY RETURN!
    }
    if ( mPWE.getEdgeListLen(el1) < mMinOverlap ||
         mPWE.getEdgeListLen(el2) < mMinOverlap )
    {
        LOG("Giving up:  Reads are too short to meet MIN_OVERLAP criterion.");
        return; // EARLY RETURN!
    }
    if ( mPWE.getGraph().getEdge(el1.front().getEdgeID()).getComponentID() !=
          mPWE.getGraph().getEdge(el2.front().getEdgeID()).getComponentID() )
    {
        LOG("Giving up:  Reads are pathed on different graph components.");
        return; // EARLY RETURN!
    }

    // first try: look for fills where the pair ends overlap
    ELItr beg1(el1.begin());
    ELItr end1(el1.end());
    ELItr beg2(el2.begin());
    ELItr end2(el2.end());
    ELItr itr(end1);
    EdgeList& closure = mClosures[pairId];
    for ( size_t idx = 1; idx <= nnn; ++idx )
    {
        using std::equal;
        if ( equal(--itr,end1,beg2) )
        {
            if ( closure.size() )
            {
                LOG("Giving up:  Ambiguous closure.");
                closure.clear();
                return; // EARLY RETURN!
            }
            if ( itr == beg1 && el1.getInitialSkip() > el2.getInitialSkip() )
            {
                LOG("Ignoring improper overlap.");
            }
            else
            {
                EdgeList tmp(el1.getInitialSkip(),el2.getFinalSkip());
                tmp.reserve((itr-beg1)+el2.size());
                tmp.insert(tmp.end(),beg1,itr);
                tmp.insert(tmp.end(),beg2,end2);
                closure = tmp;
            }
        }
    }

    // second try: look for graph walks that fill separated fragments
    long sep = el1.getFinalSkip() + el2.getInitialSkip();
    if ( sep <= mMaxSep )
    {
        EdgeList fill;
        size_t nTries = 0;
        exploreGraph(el1,el2,el1.back(),fill,nTries,sep,closure);
    }

    LibSepStats& lss = mSeps[mPM.libraryID(pairId)];
    lss.addPair();
    if ( closure.size() )
    {
        sep = mPWE.getEdgeListLen(closure)
                - mPWE.getEdgeListLen(el1)
                - mPWE.getEdgeListLen(el2);
        if ( sep > mMaxSep )
        {
            LOG("Ignoring oddly lengthy overlap.");
            closure.clear();
        }
        else
        {
            LOG("Recording closure: " << closure);
            lss.addSep(sep - K + 1);
        }
    }
}

template <unsigned K>
bool FragmentFiller<K>::exploreGraph( EdgeList const& el1, EdgeList const& el2,
                                        EdgeDesc const& cur,
                                        EdgeList& fill, size_t& nTries,
                                        long sep, EdgeList& closure )
{
    if ( ++nTries > 100 )
    {
        getMaxSegmentsBetween(el1.back(),el2.front());
        LOG("Giving up.  Tangled graph.");
        return false; // EARLY RETURN!
    }

    bool result = true;
    bool inCycle = false;
    EdgeDesc const& last = el2.front();
    for ( size_t idx = 0; result && idx < 4; ++idx )
    {
        EdgeDesc next = mPWE.getNextEdgeDesc(cur,idx);
        if ( !next.getEdgeID().isNull() )
        {
            if ( next == last )
            {
                if ( closure.size() )
                {
                    LOG("Giving up.  Ambiguous closure.");
                    closure.clear();
                    result = false;
                }
                EdgeList tmp(el1.getInitialSkip(),el2.getFinalSkip());
                tmp.reserve(el1.size()+fill.size()+el2.size());
                tmp.insert(tmp.end(),el1.begin(),el1.end());
                tmp.insert(tmp.end(),fill.begin(),fill.end());
                tmp.insert(tmp.end(),el2.begin(),el2.end());
                closure = tmp;
            }
            else if ( std::find(fill.begin(),fill.end(),next) != fill.end() )
                inCycle = true;
            else
            {
                long edgeLen = mPWE.getEdgeLen(next);
                if ( sep+edgeLen <= mMaxSep )
                {
                    fill.push_back(next);
                    LOG("Consider " << fill);
                    result=exploreGraph(el1,el2,next,fill,nTries,sep+edgeLen,closure);
                    fill.pop_back();
                }
            }
        }
    }
    if ( inCycle && closure.size() )
    {
        LOG("Giving up.  Cyclic closure.");
        closure.clear();
        result = false;
    }
    return result;
}

template <unsigned K>
size_t FragmentFiller<K>::getMaxSegmentsBetween( EdgeDesc const& ed1,
                                                 EdgeDesc const& ed2 )
{
    int result = 0;
    typedef UnipathEvidenceVec::const_iterator Itr;
    UnipathEvidenceVec const& ev1 = mPWE.getEvidence(ed1.getEdgeID());
    Itr itr1(ev1.begin());
    Itr end1(ev1.end());
    UnipathEvidenceVec const& ev2 = mPWE.getEvidence(ed2.getEdgeID());
    Itr itr2(ev2.begin());
    Itr end2(ev2.end());
    int diff;
    while ( itr1 != end1 && itr2 != end2 )
    {
        if ( itr1->getReadID() < itr2->getReadID() )
            ++itr1;
        else if ( itr2->getReadID() < itr1->getReadID() )
            ++itr2;
        else
        {
            Itr beg2(itr2);
            do
            {
                if ( ed1.getStatus()==CanonicalForm::PALINDROME &&
                      ed1.getStatus()==CanonicalForm::PALINDROME )
                {
                    diff = itr2->getSegmentID() - itr1->getSegmentID();
                    if ( diff < 0 ) diff = -diff;
                    if ( diff > result ) result = diff;
                }
                else if ( ((ed1.getStatus()==CanonicalForm::REV) == itr1->isRC()) &&
                          ((ed2.getStatus()==CanonicalForm::REV) == itr2->isRC()) )
                {
                    if ( (diff = itr2->getSegmentID()-itr1->getSegmentID()) > result )
                        result = diff;
                }
                else if ( ((ed1.getStatus()==CanonicalForm::REV) != itr1->isRC()) &&
                          ((ed2.getStatus()==CanonicalForm::REV) != itr2->isRC()) )
                {
                    if ( (diff = itr1->getSegmentID()-itr2->getSegmentID()) > result )
                        result = diff;
                }
            }
            while ( ++itr2 != end2 && itr1->getReadID() == itr2->getReadID() );
            if ( ++itr1 != end1 && itr1->getReadID() == beg2->getReadID() )
                itr2 = beg2;
        }
    }
std::cout << "Found depth max of " << result << " by poring through "
          << ev1.size() << '+' << ev2.size() << " evidence records."
          << std::endl;
    return result;
}

}


int main( int argc, char** argv )
{
    String empty;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_OrDefault(PRE,".");
    CommandArgument_String_OrDefault(DATA,".");
    CommandArgument_String_OrDefault(RUN,".");
    CommandArgument_String_OrDefault(READS_IN,"frag_reads_corr");
    CommandArgument_String_OrDefault(READS_OUT,"filled_reads");
    CommandArgument_String_OrDefault(PAIRS_OUT,empty);
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0,
            "Number of threads to use (defaults to all available processors)");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_CLIQ,64000,
            "Maximum amount of evidence to wade through.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_OVERLAP,0,
            "No. of bases overlap in valid glue.");
    CommandArgument_String_OrDefault_Doc(EXTRA_FILLS,"",
            "extra fastb file of bogus filled fragments to be appended to the "
            "file of real filled fragments - temporary hack");
    CommandArgument_Bool_OrDefault_Doc(UNIQUE_ONLY,False,
            "Emit closures only when unique.");

    CommandArgument_Double_OrDefault_Doc(MIN_FILL_FRACTION,.1,
            "The minimum fraction of successful fills to update lib stats.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_GLUE,100,
            "Maximum number of glue reads to consider.");

    CommandArgument_String_OrDefault(PAIR_IDS,empty);

    EndCommandArguments;

    unsigned const K = 28;

    // Thread control

    NUM_THREADS = configNumThreads(NUM_THREADS);

    String runDir = PRE + '/' + DATA + '/' + RUN + '/';
    String readsFile = runDir + READS_IN + ".fastb";
    String pairsFile = runDir + READS_IN + ".pairs";
    String filledFile = runDir + READS_OUT + ".fastb";

    if ( MIN_OVERLAP < K )
        MIN_OVERLAP = K;
    MIN_OVERLAP -= K - 1;

    if ( PAIRS_OUT == empty )
        PAIRS_OUT = READS_IN + "_cpd";
    PAIRS_OUT = runDir + PAIRS_OUT + ".pairs";

    FragmentFiller<K> ff(pairsFile,
                            PathCollection<K>::getInfoFilename(readsFile),
                            UnipathGraph<K>::getInfoFilename(readsFile),
                            getMaximumReadLen(readsFile)-K+1-2*MIN_OVERLAP,
                            MIN_OVERLAP, MIN_FILL_FRACTION);

    if ( PAIR_IDS != empty )
    {
        vec<longlong> ids;
        ParseLongLongSet(PAIR_IDS, ids);
        ff.processPairs(ids);
    }
    else
    {
        ff.processAllPairs();
        ff.updateLibStats().Write(PAIRS_OUT);
    }

    size_t nFills = ff.getNFills();
    PerfStat::log() << std::fixed << std::setprecision(1)
                    << PerfStat("frac_filled_pairs",
                                "% of fragment pairs that were filled",
                                100.*nFills/ff.getNPairs());

    std::cout << Date() << " Recording results." << std::endl;

    IncrementalWriter<bvec> writer(filledFile.c_str(),nFills);
    ff.writeClosures(writer);

    if ( EXTRA_FILLS != "" )
    {
        VirtualMasterVec<bvec> extras(EXTRA_FILLS.c_str());
        writer.add(extras.begin(),extras.end());
    }

    writer.close();
    PairsManager pm(writer.getNElements());
    pm.Write(filledFile.ReplaceExtension(".fastb",".pairs"));

    std::cout << Date() << " Done." << std::endl;
}
