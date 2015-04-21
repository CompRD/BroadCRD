///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/*
 * HiCDefs.cc
 *
 *  Created on: Jun 13, 2014
 *      Author: tsharpe
 */
#include <unordered_set>
#include "paths/long/hic/HiCDefs.h"

char const* HiCHitRate::gFileName = "hic.hitrates";
char const* LineInfo::gFileName = "hic.lineinfo";
char const* EdgeIdToLineIdMap::gFileName = "hic.edgelinemap";

#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class SmallVec<HiCHitRate,MempoolAllocator<HiCHitRate>>;
template class OuterVec<HiCHitRateVec>;

int const EdgeIdToLineIdMap::NOT_MAPPED;

void getHiCPairs( String const& filename, EdgeIdToLineIdMap const& map,
                    HICVec* pPairs, bool skipSelf /* = true */,
                    vec<int> lineIds /* = () */)
{
    std::unordered_set<int> lineIdsMap(lineIds.begin(), lineIds.end());

    cout << Date() << ": loading pairs and " <<
            ( skipSelf  ? "excluding" : "including" )
            << " pairs on the same line" << endl;
    cout << Date() << ": restricting to " <<
            ( lineIds.size() ? ToString(lineIds.size()) : ToString("all") )
            << " lines" << endl;

    BinaryReader::readFile(filename,pPairs);
    pPairs->erase(std::remove_if(pPairs->begin(),pPairs->end(),
            [&map,&lineIdsMap,skipSelf](EdgeLocPair const& pair)
            { int lineId1 = map[pair.el1().getEdgeId()];
              if ( lineId1 == EdgeIdToLineIdMap::NOT_MAPPED ) return true;
              int lineId2 = map[pair.el2().getEdgeId()];
              if ( lineId2 == EdgeIdToLineIdMap::NOT_MAPPED ) return true;
              if ( lineIdsMap.size() &&
                      lineIdsMap.find(lineId1) == lineIdsMap.end() &&
                      lineIdsMap.find(lineId2) == lineIdsMap.end() )
                  return true;
              return (skipSelf && lineId1 == lineId2); }),pPairs->end());
    pPairs->shrink_to_fit();
}

vec<LineLocPair> fromHICVec(HICVec const& hiCPairs, LineOffsetCalc const& offCalc)
{
    vec<LineLocPair> llp(hiCPairs.size());
#pragma omp parallel for
    for ( size_t i = 0; i < hiCPairs.size(); ++i ) {
        auto const& hiCPair = hiCPairs[i];
        llp[i] = LineLocPair( offCalc(hiCPair.el1()), offCalc(hiCPair.el2()) );
    }
    return llp;
}



int canonicalLine( int lineid, LineInfoVec const& liv )
{
    ForceAssertGe(lineid,0);
    ForceAssertLt(lineid,liv.isize());
    LineInfo const& li = liv[lineid];

    if ( li.mIsRC && li.mLineIdRC != -1 ) return li.mLineIdRC;
    else return lineid;
}


void createLineInvolution( LineInfoVec const& lineInfo, vec<int>& lineInv )
{
    lineInv.clear();
    lineInv.resize( lineInfo.size() );
    for ( int line = 0; line < lineInfo.isize(); ++line ) {
        if ( lineInfo[line].mLineIdRC >= 0 ) {
            lineInv[line] = lineInfo[line].mLineIdRC;
        }
    }
}
