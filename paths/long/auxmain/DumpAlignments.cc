/*
 * DumpAlignments.cc
 *
 *  Created on: Feb 25, 2013
 *      Author: tsharpe
 */


#include "MainTools.h"
#include "Basevector.h"
#include "Compare.h"
#include "ParallelVecUtilities.h"
#include "feudal/FilesOutputIterator.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/HashSet.h"
#include "feudal/Mempool.h"
#include "feudal/SerfVec.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/KMer.h"
#include "system/WorklistN.h"
#include "system/file/FileReader.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <utility>

namespace
{

class Context
{
    static uint8_t const NO_CONTEXT = 16;
public:
    Context() : mVal(NO_CONTEXT) {}
    Context( uint8_t prefix, uint8_t suffix )
    : mVal(4*prefix+suffix) {}

    Context rc() const { return Context(gRC[mVal]); }
    uint8_t val() const { return mVal; }

    friend bool isDifferent( Context const& c1, Context const& c2 )
    { return c1.mVal != c2.mVal ||
                c1.mVal == NO_CONTEXT || c2.mVal == NO_CONTEXT; }

    friend bool operator<( Context const& c1, Context const& c2 )
    { return c1.mVal < c2.mVal; }

    friend bool operator>( Context const& c1, Context const& c2 )
    { return c1.mVal > c2.mVal; }

    friend bool operator==( Context const& c1, Context const& c2 )
    { return c1.mVal == c2.mVal; }

    friend bool operator!=( Context const& c1, Context const& c2 )
    { return c1.mVal != c2.mVal; }

    friend ostream& operator<<( ostream& os, Context const& context )
    { return os << (int)context.mVal; }

private:
    Context( uint8_t val ) : mVal(val) {}

    uint8_t mVal;
    static uint8_t const gRC[17];
};

uint8_t const Context::gRC[17] = {15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0,16};

struct Loc
{
    static uint8_t const NULL_PREFIX = 255;

    Loc()=default;
    Loc( uint32_t readId, uint8_t offset, uint8_t len, Context context,
            bool isRC )
    : mReadId(readId), mOffset(offset), mLen(len),
      mContext(context), mIsRC(isRC) {}

    uint32_t mReadId; // index of read in input set
    uint8_t mOffset;  // offset within read
    uint8_t mLen;     // length of read
    Context mContext; // preceding and following base code
    bool mIsRC;       // kmer as observed in read is not canonical

    friend bool operator<( Loc const& loc1, Loc const& loc2 )
    { if ( loc1.mReadId < loc2.mReadId ) return true;
      if ( loc1.mReadId > loc2.mReadId ) return false;
      if ( loc1.mOffset < loc2.mOffset ) return true;
      if ( loc1.mOffset > loc2.mOffset ) return false;
      if ( loc1.mLen < loc2.mLen ) return true;
      if ( loc1.mLen > loc2.mLen ) return false;
      if ( loc1.mContext < loc2.mContext ) return true;
      if ( loc1.mContext > loc2.mContext ) return false;
      if ( loc1.mIsRC < loc2.mIsRC ) return true;
      return false; }

    friend ostream& operator<<( ostream& os, Loc const& loc )
    { if ( loc.mIsRC ) os << '~';
      return os << loc.mReadId << '@' << (int)loc.mOffset << ':' << (int)loc.mLen
              << '[' << loc.mContext << ']'; }
};

struct Alignment
{
    Alignment()=default;
    Alignment( uint32_t readId1, uint32_t readId2, int16_t offset,
                bool isRC )
    : mReadId1(readId1), mReadId2(readId2), mOffset(offset), mIsRC(isRC) {}

    uint32_t mReadId1;
    uint32_t mReadId2;
    int16_t mOffset;
    bool mIsRC;

    friend bool operator<( Alignment const& a1, Alignment const& a2 )
    { if ( a1.mReadId1 < a2.mReadId1 ) return true;
      if ( a1.mReadId1 > a2.mReadId1 ) return false;
      if ( a1.mReadId2 < a2.mReadId2 ) return true;
      if ( a1.mReadId2 > a2.mReadId2 ) return false;
      if ( a1.mOffset < a2.mOffset ) return true;
      if ( a1.mOffset > a2.mOffset ) return false;
      if ( a1.mIsRC < a2.mIsRC ) return true;
      return false; }

    friend bool operator==( Alignment const& a1, Alignment const& a2 )
    { return a1.mReadId1 == a2.mReadId1 && a1.mReadId2 == a2.mReadId2
            && a1.mOffset == a2.mOffset && a1.mIsRC == a2.mIsRC; }
};

struct Friend
{
    Friend()=default;
    Friend( uint32_t readId, int16_t offset, bool isRC )
    : mReadId(readId), mOffset(offset), mIsRC(isRC) {}

    uint32_t mReadId;
    int16_t mOffset;
    bool mIsRC;
};

template <unsigned K>
class Dict
{
public:
    struct Entry : public KMer<K>
    {
        Entry( MempoolAllocator<Entry> a ) : mLocs(a) {}

        Entry& operator=( KMer<K> const& kmer )
        { static_cast<KMer<K>&>(*this) = kmer;
          mLocs.clear(); return *this; }

        Entry& operator=( Entry&& entry )
        { static_cast<KMer<K>&>(*this) = entry;
          mLocs = std::move(entry.mLocs);
          return *this; }

        SerfVec<Loc> mLocs;
        friend ostream& operator<<( ostream& os, Entry const& entry )
        { SerfVec<Loc> locs = entry.mLocs;
          std::sort(locs.begin(),locs.end());
          os << static_cast<KMer<K>const&>(entry);
          os << ' ' << locs.size();
          for ( auto itr=locs.begin(), end=locs.end(); itr != end; ++itr )
              os << ' ' << *itr;
          return os; }
    };

    class EntryFactory
    {
    public:
        EntryFactory( MempoolAllocator<Entry> alloc ) : mAlloc(alloc) {}

        template <class X>
        MempoolAllocator<X> alloc( X* ) const
        { return MempoolAllocator<X>(mAlloc); }

        Entry* create( size_t nnn ) const
        { MempoolAllocator<Entry> alloc(mAlloc);
          Entry* result = alloc.allocate(nnn);
          for ( auto itr=result,end=result+nnn; itr != end; ++itr )
            new (itr) Entry(mAlloc);
          return result; }

        void destroy( Entry* pEnt, size_t nnn ) const
        { for ( auto itr=pEnt,end=pEnt+nnn; itr != end; ++itr )
            itr->~Entry();
          MempoolAllocator<Entry> alloc(mAlloc);
          alloc.deallocate(pEnt,nnn); }

    private:
        MempoolAllocator<Entry> mAlloc;
    };

    typedef HashSet<Entry,typename KMer<K>::Hasher,
            std::equal_to<KMer<K>>,EntryFactory> Set;

    Dict( size_t nKmers, size_t pass, size_t nPasses )
    : mPass(pass), mNPasses(nPasses),
      mSet(nKmers/nPasses,
            typename KMer<K>::Hasher(),
            std::equal_to<KMer<K>>(),
            EntryFactory(mAlloc))
    { ForceAssertLt(mNPasses,7919ul); }

    typename Set::const_iterator begin() const { return mSet.begin(); }
    typename Set::const_iterator end() const { return mSet.end(); }
    size_t size() const { return mSet.size(); }

    void processRead( size_t readId, bvec const& read )
    { unsigned len = read.size();
      ForceAssertLt(len,256u);
      if ( len == K )
        processEntry(KMer<K>(read.begin()),readId,0,len);
      else if ( len > K )
      { auto itr(read.begin(K)), end(read.end()-1);
        KMer<K> kkk(read.begin());
        uint8_t readOffset = 0;
        processEntry(kkk,readId,readOffset++,len);
        while ( itr != end )
        { uint8_t prefix = kkk.front();
          kkk.toSuccessor(*itr); ++itr;
          Context context(prefix,*itr);
          processEntry(kkk,readId,readOffset++,len,context); }
        kkk.toSuccessor(*itr);
        processEntry(kkk,readId,readOffset,len); } }

    void dumpLocs()
    { for ( auto oItr=mSet.begin(), oEnd=mSet.end(); oItr != oEnd; ++oItr )
        for ( auto iItr=oItr->begin(),iEnd=oItr->end(); iItr!=iEnd; ++iItr )
          std::cout << *iItr << '\n'; }

    template <class OItr>
    void dumpAlignments( OItr outItr )
    { for ( auto oItr=mSet.begin(), oEnd=mSet.end(); oItr != oEnd; ++oItr )
      { for ( auto iItr=oItr->begin(),iEnd=oItr->end(); iItr!=iEnd; ++iItr )
        { SerfVec<Loc> const& locs = iItr->mLocs;
          size_t nAligns = locs.size();
          if ( nAligns < 2 || nAligns > MAX_ALIGNS )
            continue;
          for ( auto itr1=locs.begin(), end=locs.end(); itr1 != end; ++itr1 )
          { bool productive = false;
            for ( auto itr2=locs.begin(); itr2 != end; ++itr2 )
            { if ( itr1 != itr2 &&
                    isDifferent(itr1->mContext,itr2->mContext) )
              { productive = true;
                uint32_t read1 = itr1->mReadId;
                uint32_t read2 = itr2->mReadId;
                if ( read1 == read2 )
                  continue;
                if ( itr1->mIsRC == itr2->mIsRC )
                { int32_t offset = itr1->mOffset - itr2->mOffset;
                  *outItr = Alignment(read1,read2,offset,false); }
                else
                { int32_t offset = itr1->mOffset-(itr2->mLen-itr2->mOffset-K);
                  *outItr = Alignment(read1,read2,offset,true); }
                ++outItr; } }
            if ( !productive ) break; } } } }

private:
    static size_t const MAX_ALIGNS = 1000;

    struct Applicator
    {
        Applicator( uint32_t readId, uint8_t offset, uint8_t len,
                        Context context, bool isRC )
        : mLoc(readId,offset,len,context,isRC) {}

        void operator()( Entry const& entry )
        { SerfVec<Loc>& locs = const_cast<Entry&>(entry).mLocs;
          if ( !locs.capacity() )
              locs.push_back(mLoc,0.f,1u);
          else
              locs.push_back(mLoc,3.f,60u); }

        Loc mLoc;
    };

    void processEntry( KMer<K> const& kkk, uint32_t readId,
                         uint8_t offset, uint8_t len,
                         Context context = Context() )
    { if ( !kkk.isRev() )
      { size_t hash = mSet.getHasher()(kkk);
        if ( thisPass(hash) )
          mSet.apply(kkk,Applicator(readId,offset,len,context,false)); }
      else
      { KMer<K> kkk2(kkk); kkk2.rc();
        size_t hash = mSet.getHasher()(kkk2);
        if ( thisPass(hash) )
          mSet.apply(kkk2,Applicator(readId,offset,len,context.rc(),true)); } }

    bool thisPass( size_t hash )
    { return (hash/7919) % mNPasses == mPass; }

    size_t mPass;
    size_t mNPasses;
    MempoolOwner<Entry> mAlloc;
    Set mSet;
};

unsigned const K = 21;

struct Batcher
{
    static size_t const BATCH_SIZE = 10000;

    Batcher( Dict<K>& dict, VirtualMasterVec<bvec> const& vbv, Dotter& dotter )
    : mDict(dict), mVBV(vbv), mDotter(dotter) {}

    void operator()( size_t batchNo )
    { size_t idx = batchNo*BATCH_SIZE;
      size_t end = std::min(idx+BATCH_SIZE,mVBV.size());
      auto itr = mVBV.begin()+idx;
      for ( size_t idx = BATCH_SIZE*batchNo; idx != end; ++idx, ++itr )
        mDict.processRead(idx,*itr);
      mDotter.batchDone(); }

    Dict<K>& mDict;
    VirtualMasterVec<bvec> mVBV;
    Dotter& mDotter;
};

size_t writeAlignments( String const& READS,
                        String const& ALIGNS_OUT,
                        String const& ALIGNS_OUT_EXT,
                        unsigned NUM_THREADS,
                        std::vector<String>* pFiles )
{
    VirtualMasterVec<bvec> reads(READS);
    size_t memFree = MemAvailable();
    size_t recsPerFile = memFree/sizeof(Alignment);
    FilesOutput filesOut(ALIGNS_OUT,ALIGNS_OUT_EXT,recsPerFile);
    size_t nKmers = 4*reads.sizeSum()/5;
    size_t memUsePerKmer = 1.5*(sizeof(Dict<K>::Entry)+1) +10.*sizeof(Loc) +.5;
    size_t nPasses = (nKmers*memUsePerKmer+memFree-1)/memFree;
    std::cout << Date() << ": creating alignments for " << reads.size()
            << " reads in " << nPasses << " passes." << std::endl;
    double clock = WallClockTime();
    for ( size_t pass = 0; pass != nPasses; ++pass )
    {
        std::cout << Date() << ": Building dictionary for pass "
                << pass << std::endl;
        Dict<K> dict(nKmers,pass,nPasses);
        size_t nBatches =
                (reads.size()+Batcher::BATCH_SIZE-1)/Batcher::BATCH_SIZE;
        Dotter dotter(nBatches);
        Batcher proc(dict,reads,dotter);
        parallelFor(0ul,nBatches,proc,NUM_THREADS);
        std::cout << Date() << ": dictionary built for pass "
                << pass << ".  Expected " << nKmers/nPasses << " kmers.  Got "
                << dict.size() << '.' << std::endl;
        //dict.dumpLocs();
        dict.dumpAlignments(filesOut.getIterator((Alignment*)nullptr));
        std::cout << Date() << ": writing alignments done for pass "
                << pass << std::endl;
    }
    std::cout << Date() << ": alignments done in "
            << TimeSince(clock) << std::endl;
    *pFiles = filesOut.getFiles();
    return reads.size();
}

void sortFile( String const& filename, vec<Alignment>& aligns, size_t nReads,
                bool UNIQUE )
{
    FileReader fr(filename);
    size_t sz = fr.getSize();
    ForceAssertEq(sz%sizeof(Alignment),0ul);
    aligns.resize( sz/sizeof(Alignment) );
    fr.read( &aligns[0], sz );
    fr.close();

    ParallelSort(aligns);
    if ( UNIQUE )
        aligns.erase(std::unique(aligns.begin(),aligns.end()),aligns.end());

    String friendFilename = filename.ReplaceExtension(".aligns",".friends");
    IncrementalWriter<SerfVec<Friend>> wtr(friendFilename,nReads);
    SerfVec<Friend> friends;
    auto itr = aligns.begin();
    auto end = aligns.end();
    for ( size_t idx = 0; idx != nReads; ++idx )
    {
        while ( itr != end )
        {
            if ( itr->mReadId1 != idx )
                break;
            friends.push_back(Friend(itr->mReadId2,itr->mOffset,itr->mIsRC));
            ++itr;
        }
        wtr.add(friends);
        friends.clear();
    }
    ForceAssert(itr==end);
}

} // end of anonymous namespace

TRIVIALLY_SERIALIZABLE(Loc);
TRIVIALLY_SERIALIZABLE(Alignment);
TRIVIALLY_SERIALIZABLE(Friend);

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(READS, "reads to process");
    CommandArgument_String_OrDefault_Doc(ALIGNS_OUT, "tmp",
      "filename template for output (.1, .2, etc. appended");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_MEMORY_GB, 0,
      "maximum amount of memory to use");
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
      "maximum number of threads to use");
    CommandArgument_Bool_OrDefault(SORT,True);
    CommandArgument_Bool_OrDefault(UNIQUE,True);
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);
    if ( MAX_MEMORY_GB )
        SetMaxMemory( (MAX_MEMORY_GB*1ul) << 30 );

    std::vector<String> files;
    size_t nReads =
            writeAlignments(READS,ALIGNS_OUT,"aligns",NUM_THREADS,&files);
/*
    if ( SORT )
    {
        vec<Alignment> aligns;
        for ( auto itr=files.begin(),end=files.end(); itr != end; ++itr )
            sortFile(*itr,aligns,nReads,UNIQUE);
    }
*/
}

#include "feudal/SmallVecDefs.h"
template class SmallVec<Loc,MempoolAllocator<Loc>>;
template class SmallVec<Friend,MempoolAllocator<Friend>>;
