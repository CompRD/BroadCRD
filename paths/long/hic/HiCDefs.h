///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#ifndef _HICDEFS_H
#define _HICDEFS_H

#include "Vec.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include <utility>
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include <unordered_map>
#include "TokenizeString.h"


class TaggedInt {
public:
    static const unsigned cId = ~0u >> 1;
    static const unsigned cTag = ~cId;

    static const int maxVal = cId;

    TaggedInt() : mIdTag(0u) {};
    TaggedInt(int lineid, bool tag) { setIdTag(lineid,tag); }
    void setId( int lineid ) {
        ForceAssertGe(lineid,0);
        bool tag = isTag();
        mIdTag = (mIdTag & cTag ) | (static_cast<unsigned>(lineid) & cId);
//        setTag(tag);
        ForceAssertEq(tag,isTag());
    }
    void setTag( bool tag ) {
        if ( tag ) mIdTag |= cTag;
        else mIdTag &= cId;
    }
    void flipTag() { setTag( !isTag() ); }
    void setIdTag(int lineid, bool tag) { setId(lineid); setTag(tag); }
    int getId() const { return mIdTag & cId; }
    bool isTag() const { return mIdTag & cTag; }

    friend ostream& operator<<( ostream& os, TaggedInt const& t ) {
        return os << t.getId() << (t.isTag()?"'":"");
    }

    // the following kind of makes you wish that scalar types were
    // first class objects, no?
    friend bool operator<(TaggedInt const& t1, TaggedInt const& t2 ) {
        return (t1.mIdTag < t2.mIdTag);
    }

    friend bool operator==(TaggedInt const& t1, TaggedInt const& t2 ) {
        return (t1.mIdTag == t2.mIdTag);
    }

private:
    unsigned mIdTag;
};
typedef TaggedInt TaggedLineId;
TRIVIALLY_SERIALIZABLE(TaggedInt);

// Location on an edge.
// Please note that there is code in the wild that uses pair<int,int> instead
// of this structure, so you'll have to hunt that code down if you change this.
class EdgeLoc
{
public:
    EdgeLoc() : mEdgeId(-1), mPos(-1) {}
    EdgeLoc( int edgeId, int pos ) : mEdgeId(edgeId), mPos(pos) {}

    int getEdgeId() const { return mEdgeId; }
    int getPos() const { return mPos; }

#if 0
    friend bool operator==( EdgeLoc const& el1, EdgeLoc const& el2 )
    { return el1.mEdgeId == el2.mEdgeId && el1.mPos == el2.mPos; }
    friend bool operator<( EdgeLoc const& el1, EdgeLoc const& el2 )
    { if ( el1.mEdgeId != el2.mEdgeId ) return el1.mEdgeId < el2.mEdgeId;
      return el1.mPos < el2.mPos; }
    friend ostream& operator<<( ostream& os, EdgeLoc const& el )
    { return os << el.mEdgeId << ':' << el.mPos; }
#endif

private:
    int mEdgeId;
    int mPos;
};
TRIVIALLY_SERIALIZABLE(EdgeLoc);

// A pair of EdgeLocs.
class EdgeLocPair : public std::pair<EdgeLoc,EdgeLoc>
{
public:
    EdgeLocPair() = default;
    EdgeLocPair( EdgeLoc const& el1, EdgeLoc const& el2 )
    : std::pair<EdgeLoc,EdgeLoc>(el1,el2) {}

    EdgeLoc const& el1() const { return first; }
    EdgeLoc const& el2() const { return second; }
};
TRIVIALLY_SERIALIZABLE(EdgeLocPair);

class LineLoc
{
public:
    LineLoc() : mLineId(-1), mLeftOffset(-1), mRightOffset(-1) {}
    LineLoc( int edgeId, int left, int right ) : mLineId(edgeId),
           mLeftOffset(left), mRightOffset(right) {}

    int getLineId() const { return mLineId; }
    int getLeftOffset() const { return mLeftOffset; };
    int getRightOffset() const { return mRightOffset; };

    friend ostream& operator<<(ostream& out, LineLoc const& ll) {
        out << ll.mLineId << " [" << ll.mLeftOffset << "|"
                << ll.mRightOffset << "]";
        return out;
    }

private:
    int mLineId;
    int mLeftOffset;
    int mRightOffset;
};
TRIVIALLY_SERIALIZABLE(LineLoc);


// This is the object pickled in a "hic.*.aligns" file.
typedef vec<EdgeLocPair> HICVec;

// Please note that there is code in the wild that uses pair<int,int> instead
// of this structure, so you'll have to hunt that code down if you change this.
class GenomeLoc
{
public:
    GenomeLoc() : mChrId(-1), mOffset(0) {}
    GenomeLoc( int chrId, int offset ) : mChrId(chrId), mOffset(offset) {}

    bool isNull() const { return mChrId == -1; }
    int getChrId() const { return mChrId; }
    int getOffset() const { return mOffset; }

    friend bool operator<( GenomeLoc const& gl1, GenomeLoc const& gl2 )
    { if ( gl1.mChrId != gl2.mChrId ) return gl1.mChrId < gl2.mChrId;
      return gl1.mOffset < gl2.mOffset; }

    friend ostream& operator<<( ostream& os, GenomeLoc const& gl )
    { return os << gl.mChrId << ':' << gl.mOffset; }

private:
    int mChrId;
    int mOffset;
};
TRIVIALLY_SERIALIZABLE(GenomeLoc);

// This is the object in the ".aligns" file written by GapToy.
// Each inner vector shows alternate alignments for an edge.
typedef vec<GenomeLoc> GenomeLocVec;
typedef vec<GenomeLocVec> VecGenomeLocVec;

// a more-or-less random collection of stuff about a line
struct LineInfo
{
    LineInfo()
    : mIsRC(false),mLineIdRC(-1),mLineLen(0),mFirstEdgeId(-1),mLastEdgeId(-1) {}

    GenomeLoc mLoc;
    bool mIsRC; // if line is not oriented with reference
    int mLineIdRC; // id of the RC line
    int mLineLen; // length of the line
    int mFirstEdgeId; // first edge of the line
    int mLastEdgeId; // last edge of the line

    static char const* gFileName; // name for file containing a LineInfoVec
};
TRIVIALLY_SERIALIZABLE(LineInfo);
typedef vec<LineInfo> LineInfoVec;

struct HiCHitRate
{
    HiCHitRate()=default;
    HiCHitRate( int lineId, float hitsPerKb )
    : mLineId(lineId), mHitsPerKb(hitsPerKb) {}

    int mLineId;
    float mHitsPerKb;

    static char const* gFileName; // name for file containing a VecHiCHitRateVec

    friend std::ostream& operator<<( std::ostream& out, HiCHitRate const& h ) {
        out << h.mHitsPerKb << "[" << h.mLineId << "]";
        return out;
    }

};
TRIVIALLY_SERIALIZABLE(HiCHitRate);

typedef SerfVec<HiCHitRate> HiCHitRateVec;
typedef MasterVec<HiCHitRateVec> VecHiCHitRateVec;
extern template class OuterVec<HiCHitRateVec>;

class SepHistogram
{
public:
    // bins are fractional sizes; we calculate the
    // bin centers and add weight to the left and right bins
    SepHistogram() : mLow(0), mHi(0), mBins(0), mBinSize(0.),
            mTotal(0), mAnalyticSum(0.){}
    SepHistogram( float low, float hi, unsigned nbins = 200 ) :
        mLow(low), mHi(hi), mBins(nbins),
        mTotal(nbins), mCounts(nbins,1u), mAnalyticSum(0.) {
        float binsize = ( mHi - mLow ) / static_cast<double>(mBins);
        mHi = mLow + binsize * mBins;   // account for roundoff
        mBinSize = binsize;

        FixAnalyticSum();
        PRINT4(mLow,mHi,mBins,mBinSize);
    }

    void CloneShape( SepHistogram const& that ) {
        this->mLow = that.mLow;
        this->mHi = that.mHi;
        this->mBinSize = that.mBinSize;
        this->mBins = that.mBins;
        this->mCounts.resize(this->mBins,1u);
        this->mTotal=this->mBins;
    }

    float operator()(float value) const {
        if ( value < mLow || value > mHi ) return 0.;
        return mCounts[ValueToBin(value)] / mTotal;
    }

    void FixAnalyticSum() {
        for ( size_t i = 0; i < mBins; ++i )
            mAnalyticSum += 1.0 / this->BinToCenterValue(i);
    }

    float Analytic(float value) const {
        if ( value < mLow || value > mHi ) return 0.;
        double p = 1.0 / value;
        return (p / mAnalyticSum);
    }

    // no clamping
    size_t ValueToBin( float value ) const {
        ForceAssertLe( value, mHi );
        ForceAssertGe( value, mLow );
        return ( value - mLow ) / mBinSize;
    }

    float BinToCenterValue( unsigned bin ) const {
        ForceAssertGe( bin, 0u );
        ForceAssertLt( bin, mBins );
        return mLow + ( static_cast<float>(bin) + 0.5 )*mBinSize;
    }

    float ValueWeight( float value, size_t bin ) const {
        float tlow = bin * mBinSize + mLow;
        float thigh = (bin+1) * mBinSize + mLow;
        if ( value < tlow || value > thigh ) return 0.;
        return (thigh - value ) / mBinSize;
    }

    // returns 1 if one is added and 0 otherwise
    size_t Add( float value ) {
        if ( value < mLow || value >= mHi ) return 0;
        size_t bin = ValueToBin(value);
        if ( bin >= mCounts.size() ) {
            PRINT8("OVER",bin,value,mLow,mHi,mHi-value,mBinSize,mBins);
        }
        ForceAssertLt( bin, mCounts.size() );
        if ( bin == mBins - 1 )
            mCounts[bin] += 1.0;
        else {
            float w0 = ValueWeight(value, bin);
            float w1 = 1.0 - w0;
            mCounts[bin]+=w0;
            mCounts[bin+1] += w1;
        }
        mTotal+=1;
        return 1;
    }

    void Smooth(float bw) {
        ForceAssertGt(bw,0);
        int mult = 4;
        int ksize = mult*2*std::ceil(bw)+1;
        int hsize = (ksize - 1) / 2;
        vec<float> kernel;
        for ( int i = 0; i < ksize; ++i ) {
            float x = i - (ksize-1)/2;
            kernel.push_back( exp( -1 * x*x/(2*bw*bw) ) );
        }
        double sum = Sum(kernel);
        for ( int i = 0; i < ksize; ++i ) kernel[i] /= sum;
//        cout << "smoothing kernel:" << endl;
//        for ( int i = 0; i < ksize; ++i )
//            cout << i << " (" << i-(ksize-1)/2 << "): " << kernel[i] << endl;
        vec<float> ksums(ksize);
        ksums[0] = kernel[0];
        for ( int i = 1; i < ksize; ++i )
            ksums[i] = ksums[i-1] + kernel[i];

        auto tmp(mCounts);
        ForceAssertGt(mCounts.isize(), hsize);
        int i = 0;
        for ( ; i < hsize; ++i ) {
            double val = 0.;
            double sum = 0.;
            for ( int j = -i; j <= hsize; ++j ) {
                val += tmp[i+j]*kernel[j+hsize];
                sum += kernel[j+hsize];
            }
            mCounts[i] = val / sum;
        }
        for ( ; i < mCounts.isize() - hsize; ++i ) {
            double val = 0.;
            for ( int j = -hsize; j <= hsize; ++j )
                val += tmp[i+j]*kernel[j+hsize];
            mCounts[i] = val;
        }
        for ( ; i < mCounts.isize(); ++i ) {
            double val = 0.;
            double sum = 0.;
            for ( int j = -hsize; j <= mCounts.isize()-i-1; ++j ) {
                val += tmp[i+j]*kernel[j+hsize];
                sum += kernel[j+hsize];
            }
            mCounts[i] = val / sum;
        }

    }

    friend double chisq(
            SepHistogram const& query,
            SepHistogram const& model,
            float lower = 0.,
            float upper = std::numeric_limits<float>::max() )
    {
        if ( upper > model.mHi ) upper = model.mHi;
        if ( lower < model.mLow ) lower = model.mLow;
        size_t bupper, blower;
        bupper = query.ValueToBin(upper);
        blower = query.ValueToBin(lower);

        long double qtotal = 0., mtotal = 0.;
        for ( size_t i = blower; i < query.mBins && i <= bupper; ++i ) {
            qtotal += query.mCounts[i];
//            mtotal += model.mCounts[i];
            auto centerVal = query.BinToCenterValue(i);
            auto binHalfWidth = query.mBinSize / 2.0;
            double k2 = binHalfWidth + centerVal;
            double k1 = binHalfWidth - centerVal;
            mtotal += 1/(k1*k1) - 1/(k2*k2);
        }

        long double chis = 0.;
        for ( size_t i = blower; i < query.mBins && i <= bupper; ++i ) {
            long double obse = query.mCounts[i] / qtotal;
            auto centerVal = query.BinToCenterValue(i);
            auto binHalfWidth = query.mBinSize / 2.0;
            double k2 = binHalfWidth + centerVal;
            double k1 = binHalfWidth - centerVal;
            long double expc = ( 1/(k1*k1) - 1/(k2*k2) ) / mtotal;

            if ( expc > 1e-38 ) chis += (obse-expc)*(obse-expc)/expc;
        }

        return chis;
    }


    friend double relative_entropy(
            SepHistogram const& query,
            SepHistogram const& model,
            float lower = 0.,
            float upper = 0.,
            int trace_id = -1)
    {
        bool trace = ( trace_id > -1 );
        // for the moment we require the distributions to be similarly
        // structured.  that will change someday because it's whack
        ForceAssertEq(query.mBinSize, model.mBinSize);   // precise :(
        ForceAssertEq(query.mLow,     model.mLow);
        ForceAssertEq(query.mBins,    model.mBins);
        size_t bupper, blower;
        if ( upper != 0 && upper > model.mHi ) upper = model.mHi;
        if ( lower != 0 && lower < model.mLow ) lower = model.mLow;
        if ( upper == 0. ) bupper = query.mBins - 1;
        else bupper = query.ValueToBin(upper);
        if ( lower == 0. ) blower = 0;
        else blower = query.ValueToBin(lower);

        double qtotal = 0., mtotal = 0.;
        for ( size_t i = blower; i < query.mBins && i <= bupper; ++i ) {
            qtotal += query.mCounts[i];
            mtotal += model.mCounts[i];
        }

        ForceAssertGt(qtotal, 0);
        ForceAssertGt(mtotal, 0);

        long double entropy = 0.;
        double const epsilon = 1e-10;
//        PRINT2(blower,bupper);
        for ( size_t i = blower; i < query.mBins && i <= bupper; ++i ) {
            long double p = model.mCounts[i] / mtotal;
            long double q = query.mCounts[i] / qtotal;
            // you're getting something for free, here,
            // but it shouldn't matter if you've set things up correctly
            if ( q < epsilon ) continue;
#if 0
            if ( (p / q) < epsilon ) {
                cout << "WARNING: skipping query bin " << i << " value "
                        << query.BinToCenterValue(i) <<
                        "p=" << p << ", q=" << q << endl;
                continue;
            }
#endif
            entropy += p * ( log( p ) - log ( q ) );
            if ( trace ) cout << "ENTROPY-" << trace_id << " " << query.mCounts[i] << " "<< entropy << endl;
        }

        return entropy;
    }

    void writeBinary( BinaryWriter& writer ) const {
        writer.write(mLow);
        writer.write(mHi);
        writer.write(mBins);
        writer.write(mBinSize);
        writer.write(mTotal);
        writer.write(mCounts);
    }
    void readBinary( BinaryReader& reader ) {
        reader.read(&mLow);
        reader.read(&mHi);
        reader.read(&mBins);
        reader.read(&mBinSize);
        reader.read(&mTotal);
        reader.read(&mCounts);
    }
    static size_t externalSizeof() { return 0; }

    void Print( std::ostream& out, String head = "" ) const {
        for ( size_t i = 0; i < mCounts.size(); ++i ) {
            if ( head != "" ) out << head;
            out << i*mBinSize + mLow + (mBinSize/2.)<< " " << mCounts[i] << endl;
        }
    }

    void PrintStats( std::ostream& out, String head = "" ) const {
        if ( head != "" ) out << head;
        PRINT5_TO(out,mLow,mHi,mBins,mBinSize,mTotal);
    }

private:
    float mLow;
    float mHi;
    unsigned mBins;
    float mBinSize;
    unsigned mTotal;
    vec<float> mCounts;
    double mAnalyticSum;
};
SELF_SERIALIZABLE(SepHistogram);

class EdgeIdToLineIdMap
{
public:
    static int const NOT_MAPPED = -1;

    void clear() { mMap.clear(); }
    void resize( size_t nEdges ) { mMap.resize(nEdges,NOT_MAPPED); }

    int operator[]( size_t edgeId ) const { return mMap[edgeId]; }
    bool isUnmapped( size_t edgeId ) const
    { return mMap[edgeId] == NOT_MAPPED; }

    void set( size_t edgeId, int lineId )
    { mMap[edgeId] = lineId; }

    void writeBinary( BinaryWriter& writer ) const { writer.write(mMap); }
    void readBinary( BinaryReader& reader ) { reader.read(&mMap); }

    static char const* gFileName; // name for file containing a EdgeIdToLineIdMap
private:
    vec<int> mMap;
};
SELF_SERIALIZABLE(EdgeIdToLineIdMap);

// LineOffsetCalc -- convert a position on an edge to a position in a line.
class LineOffsetCalc {
public:
    typedef size_t EdgeId;
    typedef int LineId;
    typedef size_t SegId;
    typedef size_t PathId;
    typedef std::pair<SegId,PathId> SegPath;

    // READ THIS -- if you need this class to be thread safe,
    // then before you make a query, you must preload the cache
    // for every line you are interested in by calling cachePreload(lineId).
    enum ThreadSafety { NOT_THREAD_SAFE, THREAD_SAFE };

    LineOffsetCalc() = delete;
    LineOffsetCalc( HyperBasevectorX const& hbx, LineVec const& lines,
            vec<int> inv,
            EdgeIdToLineIdMap const& edgeToLineMap,
            HBVXKmerRuler const& ruler,
            ThreadSafety needThreadSafe = NOT_THREAD_SAFE ) : mHBX(hbx), mLines(lines),
            mInv(inv), mEdgeToLineMap(edgeToLineMap), mRuler(ruler),
            mNeedThreadSafe(needThreadSafe)
    {
    };

    LineLoc operator()(EdgeLoc const& el) const {
        return this->operator()(el.getEdgeId(), el.getPos());
    }

    LineLoc operator()(EdgeId e, size_t offset = 0) const {
        int lineId = mEdgeToLineMap[e];
        ForceAssertNe( lineId, EdgeIdToLineIdMap::NOT_MAPPED );
        if (mNeedThreadSafe == NOT_THREAD_SAFE) checkCache(lineId);
        ForceAssertLt(e,mInv.size());
        if ( mEdgeMap.find( lineEdge(lineId,e) ) == mEdgeMap.end() ) e = mInv[e];
        ForceAssertLt(lineId, mLines.isize());
        auto const& line = mLines[lineId];

        SegPath seg_path;
        try {
            seg_path = mEdgeMap.at( lineEdge(lineId, e));
        } catch (std::out_of_range &err) {
#pragma omp critical
            cout << "cached did not include line,edge=" << lineId << "," << e << endl;
            throw;
        }

        for ( size_t segId = 0; segId < line.size(); ++segId ) {
                 if ( segId != seg_path.first )
                     offset += mSegLengthCache.at(lineSeg(lineId,segId));
                 else {
                     auto pathId = seg_path.second;
                     offset += GetPathLength(line[segId][pathId], mRuler, e);
                     break;
                 }
        }
        return LineLoc(lineId,offset,mLineLenCache[lineId] - offset);
    }

    void cachePreload(LineId lineId) const {
        checkCache(lineId);
    }

private:
    std::pair<LineId,SegId> lineSeg(LineId lineid, SegId segid ) const {
        return std::make_pair( lineid, segid );
    }
    std::pair<LineId,EdgeId> lineEdge(LineId lineid, EdgeId edgeid ) const {
        return std::make_pair( lineid, edgeid );
    }
    void checkCache( LineId lineId ) const {
        if ( mLineLenCache.find(lineId) != mLineLenCache.end() ) return;

        int lineLen =  0;
        auto const& line = mLines[lineId];
        for ( SegId segId = 0; segId < line.size(); ++segId ) {
            auto segLen = GetSegmentLength(line[segId], mRuler);
            mSegLengthCache[ lineSeg(lineId, segId) ] = segLen;
            lineLen += segLen;
            for ( size_t pathId = 0; pathId < line[segId].size(); ++pathId )
                for ( auto const& edgeId : line[segId][pathId] )
                    mEdgeMap[ lineEdge(lineId,edgeId) ] =
                            std::make_pair(segId, pathId);
        }
        mLineLenCache[lineId]=lineLen;

    }

    HyperBasevectorX const& mHBX;
    LineVec const& mLines;
    vec<int> mInv;
    EdgeIdToLineIdMap const& mEdgeToLineMap;
    HBVXKmerRuler const& mRuler;
    ThreadSafety mNeedThreadSafe;
    // line,seg -> seg length
    std::map< std::pair<LineId,SegId>, size_t > mutable mSegLengthCache;
    std::map< std::pair<LineId,EdgeId>, SegPath> mutable mEdgeMap;
    std::unordered_map< LineId, int > mutable mLineLenCache;
};


// A pair of LineLocs.
class LineLocPair : public std::pair<LineLoc,LineLoc>
{
public:
    LineLocPair() = default;
    LineLocPair( LineLoc const& el1, LineLoc const& el2 )
    : std::pair<LineLoc,LineLoc>(el1,el2) {}

    LineLoc const& ll1() const { return first; }
    LineLoc const& ll2() const { return second; }
};
TRIVIALLY_SERIALIZABLE(LineLocPair);

typedef vec<LineLocPair> LineLocPairVec;



void getHiCPairs( String const& filename, EdgeIdToLineIdMap const& map,
        HICVec* pPairs, bool skipSelf = true, vec<int> lineid = vec<int>() );

vec<LineLocPair> fromHICVec(HICVec const& hiCPairs, LineOffsetCalc const& offCalc);

int canonicalLine( int lineid, LineInfoVec const& liv );

void createLineInvolution( LineInfoVec const& lineInfo, vec<int>& lineInv );

#if 0
class PerLineEdgeIndex {
public:
    typedef std::pair<size_t,size_t> SegPath;
    typedef size_t EdgeId;

    PerLineEdgeIndex() = delete;
    PerLineEdgeIndex(int lineId, Line const& line, HyperBasevectorX const& hbx) :
        mLineId(lineId), mLine(line), mHBX(hbx), mRuler(mHBX) {
        buildMaps();
    }

    LineLoc operator()(EdgeId e, size_t offset = 0) const {
        SegPath seg_path = mEdgeMap.at(e);
        for ( size_t segId = 0; segId < mLine.size(); ++segId ) {
            if ( segId != seg_path.first )
                offset += mSegLengths[segId];
            else {
                auto pathId = seg_path.second;
                offset += GetPathLength(mLine[segId][pathId], mRuler, e);
                break;
            }
        }
        return LineLoc(mLineId,offset);
    }

    bool has(EdgeId e) const {
        return ( mEdgeMap.find(e) != mEdgeMap.end() );
    }

    int getLineLength() {
        return GetLineLength(mLine, mRuler);
    }

    void Print( ostream& out ) {
        cout << "EdgeMap for LineId " << mLineId << endl;
        for ( auto itr = mEdgeMap.begin(); itr != mEdgeMap.end(); ++itr )
            cout << itr->first << " --> " << "(" << itr->second.first
                << "," << itr->second.second << ")" << endl;
    }

private:
    void buildMaps() {
        mSegLengths.resize( mLine.size() );
        for ( size_t segId = 0; segId < mLine.size(); ++segId ) {
            mSegLengths[segId] = GetSegmentLength(mLine[segId], mRuler);
            for ( size_t pathId = 0; pathId < mLine[segId].size(); ++pathId )
                for ( auto const& edgeId : mLine[segId][pathId] ) {
                    mEdgeMap[edgeId] = std::make_pair( segId, pathId );
                }
        }
    }

    int mLineId;
    Line const& mLine;
    HyperBasevectorX const& mHBX;
    HBVXKmerRuler mRuler;
    // size_t EdgeId -> < SegmentId, PathId >
    std::unordered_map<EdgeId, SegPath> mEdgeMap;
    vec<size_t> mSegLengths;
};
#endif

#endif // _HICDEFS_H
