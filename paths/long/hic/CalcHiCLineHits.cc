///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * CalcHiCLineHits.cc
 *
 *  Created on: Jun 7, 2014
 *      Author: tsharpe
 */

#include "MainTools.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "paths/long/hic/HiCDefs.h"
#include "paths/long/large/Lines.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <unordered_map>
#include <utility>


namespace
{

struct HBVXMboIRuler
{ HBVXMboIRuler(HyperBasevectorX const& hbvx ) : mHBVX(hbvx) {}
  int operator()( int edgeId ) {
      auto const& edge = mHBVX.EdgeObject(edgeId);
      const String MboI = "GATC";
      const int hSize = MboI.size();
      auto const bMboI = BaseVec(MboI);

      int count = 0;
      auto itr = edge.cbegin();

      while ( itr != edge.cend() ) {
          auto hitr = bMboI.cbegin();
          while ( hitr != bMboI.cend() && itr != edge.cend() ) {
              if ( *hitr != *itr ) { break; }
              hitr++; itr++;
          }
          if ( hitr == bMboI.cend() ) { count++; }
          if ( itr != edge.cend() ) itr++;
      }


      return count;
  }
  HyperBasevectorX const& mHBVX;
};

inline void GetLineLengthsMboI( HyperBasevectorX const& hb, LineVec const& lines,
                                vec<int>& llens )
{ GetLineLengths(lines,llens,HBVXMboIRuler(hb)); }


GenomeLoc getLineLocus( Line const& line, VecGenomeLocVec const& aligns )
{
    for ( LineSegment const& seg : line )
        for ( LinePath const& path : seg )
            for ( int edgeId : path )
            {
                vec<GenomeLoc> const& locs = aligns[edgeId];
                if ( !locs.empty() )
                    return locs.front();
            }
    return GenomeLoc();
}

size_t buildMap( String const& dir,
                    EdgeIdToLineIdMap* pMap, LineInfoVec* pInfoVec, vec<int>* pRestr )
{
    String alignsFile = dir + "/a.aligns";
    String linesFile = dir + "/a.s.lines";
    String invFile = dir + "/a.s.inv";
    String hbvxFile = dir + "/a.s.hbx";
    if ( !isReadable(linesFile) ||
            !isReadable(invFile) ||
            !isReadable(hbvxFile) )
    {
        linesFile = dir + "/a.lines";
        invFile = dir + "/a.inv";
        hbvxFile = dir + "/a.hbx";
    }

    LineVec lines;
    BinaryReader::readFile(linesFile,&lines);
    size_t nLines = lines.size();
    if ( !nLines )
    {
        std::cout << "That's funny -- there are no lines." << std::endl;
        exit(1);
    }
    LineInfoVec& infoVec = *pInfoVec;
    infoVec.clear();
    infoVec.resize(nLines);

    size_t nEdges = 0;
    for ( Line const& line : lines )
        for ( LineSegment const& seg : line )
            for ( LinePath const& path : seg )
                for ( unsigned edgeId : path )
                    if ( edgeId >= nEdges )
                        nEdges = edgeId+1;

    vec<int> inv;
    BinaryReader::readFile(invFile,&inv);
    ForceAssertEq(inv.size(),nEdges);

    EdgeIdToLineIdMap& map = *pMap;
    map.clear();
    map.resize(nEdges);

    VecGenomeLocVec aligns;
    BinaryReader::readFile(alignsFile,&aligns);
    ForceAssertLe(aligns.size(),nEdges);
    aligns.resize(nEdges);

    HyperBasevectorX hbvx;
    BinaryReader::readFile(hbvxFile,&hbvx);

    GetLineLengthsMboI(hbvx,lines,*pRestr);
    for ( size_t i = 0; i < pRestr->size(); ++i )
        (*pRestr)[i] += 1;

    vec<int> llens;
    GetLineLengths(hbvx,lines,llens);

    for ( size_t lineId = 0; lineId != nLines; ++lineId )
    {
        Line const& line = lines[lineId];
        LineInfo& info = infoVec[lineId];
        info.mLoc = getLineLocus(line,aligns);
        info.mLineLen = llens[lineId];
        info.mFirstEdgeId = line.front().front().front();
        info.mLastEdgeId = line.back().back().back();
        int firstEdgeIdRC = inv[info.mFirstEdgeId];
        int lastEdgeIdRC = inv[info.mLastEdgeId];
        ForceAssertNe(firstEdgeIdRC,-1);
        ForceAssertNe(lastEdgeIdRC,-1);

        if ( map.isUnmapped(firstEdgeIdRC) )
        {
            if ( !map.isUnmapped(lastEdgeIdRC) )
                FatalErr("In processing " << lineId << " edge " << lastEdgeIdRC
                        << " is already mapped to " << map[lastEdgeIdRC]);
        }
        else
        {
            int lineIdRC = map[firstEdgeIdRC];
            ForceAssertEq(map[lastEdgeIdRC],lineIdRC);
            info.mLineIdRC = lineIdRC;
            LineInfo& infoRC = infoVec[lineIdRC];
            infoRC.mLineIdRC = lineId;
            if ( infoRC.mLoc.isNull() && !info.mLoc.isNull() )
            {
                infoRC.mLoc = info.mLoc;
                infoRC.mIsRC = !info.mIsRC;
            }
            else if ( info.mLoc.isNull() && !infoRC.mLoc.isNull() )
            {
                info.mLoc = infoRC.mLoc;
                info.mIsRC = !infoRC.mIsRC;
            }
            else if ( !info.mLoc.isNull() && !infoRC.mLoc.isNull() &&
                       info.mLoc.getChrId() == infoRC.mLoc.getChrId() )
            {
                info.mIsRC = (info.mLoc.getOffset() < infoRC.mLoc.getOffset());
                infoRC.mIsRC = !info.mIsRC;
            }
        }
        map.set(info.mFirstEdgeId,lineId);
        map.set(info.mLastEdgeId,lineId);
    }

    for ( size_t lineId = 0; lineId != nLines; ++lineId )
    {
        int lineIdToAssign = lineId;
        LineInfo const& info = infoVec[lineId];
        if ( info.mIsRC && info.mLineIdRC != -1 )
            lineIdToAssign = info.mLineIdRC;
        for ( LineSegment const& seg : lines[lineId] )
            for ( LinePath const& path : seg )
                for ( int edgeId : path )
                {
                    map.set(edgeId,lineIdToAssign);
                    edgeId = inv[edgeId];
                    if ( edgeId != -1 )
                        map.set(edgeId,lineIdToAssign);
                }
    }

    return nLines;
}


class HitProc
{
public:
    HitProc( size_t batchSize, EdgeIdToLineIdMap const& edgeToLineMap,
                HICVec const& hiCPairs,
                vec<std::unordered_map<int,unsigned>>& results )
    : mBatchSize(batchSize), mEdgeToLineMap(edgeToLineMap), mHiCPairs(hiCPairs),
      mResults(results) {}

    void operator()( size_t batchNo )
    {
        size_t nLines = mResults.size();
        int minLineId = std::min(nLines,batchNo*mBatchSize);
        int endLineId = std::min(nLines,minLineId+mBatchSize);
        mCounts.clear();
        mCounts.resize(endLineId-minLineId);
        for ( EdgeLocPair const& pair : mHiCPairs )
        {
            int lineId1 = mEdgeToLineMap[pair.el1().getEdgeId()];
            int lineId2 = mEdgeToLineMap[pair.el2().getEdgeId()];
            if ( lineId1 >= minLineId && lineId1 < endLineId )
                mCounts[lineId1-minLineId][lineId2] += 1u;
            if ( lineId2 >= minLineId && lineId2 < endLineId )
                mCounts[lineId2-minLineId][lineId1] += 1u;
        }
        auto rItr = mResults.begin()+minLineId;
        for ( auto& map : mCounts )
        {
            using std::swap;
            swap(map,*rItr);
            ++rItr;
        }
    }

private:
    size_t mBatchSize;
    EdgeIdToLineIdMap const& mEdgeToLineMap;
    HICVec const& mHiCPairs;
    vec<std::unordered_map<int,unsigned>>& mResults;
    vec<std::unordered_map<int,unsigned>> mCounts;
};

inline float perKb( unsigned hits, int len )
{ return 1024.*hits/len; }

}

int main( int argc, char** argv )
{
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Doc(DIR,"assembly dir.");
    CommandArgument_String_Doc(HICPAIRS,"Hi-C edge-pairing file.");
    CommandArgument_String_Doc(OUTDIR,"directory for output.");
    EndCommandArguments;

    EdgeIdToLineIdMap edgeToLineMap;
    LineInfoVec lineInfoVec;
    vec<int> restr;
    size_t nLines = buildMap( DIR, &edgeToLineMap, &lineInfoVec, &restr );

    BinaryWriter::writeFile(OUTDIR+'/'+LineInfo::gFileName,lineInfoVec);
    BinaryWriter::writeFile(OUTDIR+'/'+EdgeIdToLineIdMap::gFileName,
            edgeToLineMap);

    HICVec hiCPairs;
    getHiCPairs(HICPAIRS,edgeToLineMap,&hiCPairs);
    if ( hiCPairs.empty() )
    {
        std::cout << "That's funny -- the HiC pairs file is empty."
                  << std::endl;
        exit(1);
    }

    vec<std::unordered_map<int,unsigned>> hitMaps(nLines);
    size_t const BATCHSIZE = 65536;
    size_t nBatches = (nLines+BATCHSIZE-1)/BATCHSIZE;
    if ( nBatches )
    {
        HitProc proc(BATCHSIZE,edgeToLineMap,hiCPairs,hitMaps);
        parallelFor(0ul,nBatches,proc);
    }

    VecHiCHitRateVec hitRates(nLines);
    auto oItr = hitRates.begin();
    for ( std::unordered_map<int,unsigned> const& map : hitMaps )
    {
        if ( map.size() )
        {
            HiCHitRateVec& oVec = *oItr;
            oVec.reserve(map.size());
            for ( auto const& entry : map )
            {
                LineInfo const& info = lineInfoVec[entry.first];
                float hitRate = perKb(entry.second,info.mLineLen);
//                float hitRate = entry.second / restr[entry.first];
                oVec.push_back(HiCHitRate(entry.first,hitRate));
            }
            std::sort(oVec.begin(),oVec.end(),
                    []( HiCHitRate const& hr1, HiCHitRate const& hr2 )
                    { return hr1.mHitsPerKb > hr2.mHitsPerKb; });
        }
        ++oItr;
    }

    hitRates.WriteAll(OUTDIR+'/'+HiCHitRate::gFileName);
}
