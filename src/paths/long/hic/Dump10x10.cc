///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * Dump10x10.cc
 *
 *  Created on: Jun 18, 2014
 *      Author: tsharpe
 */
#include "MainTools.h"
#include "paths/long/hic/HiCDefs.h"
#include <sstream>
#include <unordered_map>

namespace
{

inline bool isTarget( LineInfo const& info, int minLen )
{
    return !info.mIsRC && info.mLineLen >= minLen && !info.mLoc.isNull();
}

void setTargets( LineInfoVec const& lineInfoVec, int minLen,
                    vec<int>* pTargets )
{
    pTargets->clear();
    pTargets->reserve(lineInfoVec.size());

    int lineId = 0;
    auto end = lineInfoVec.end();
    for ( auto itr=lineInfoVec.begin(); itr != end; ++itr,++lineId )
        if ( isTarget(*itr,minLen) )
            pTargets->push_back(lineId);

    std::sort(pTargets->begin(),pTargets->end(),
                [&lineInfoVec]( int id1, int id2 )
                { return lineInfoVec[id1].mLoc < lineInfoVec[id2].mLoc; });
}

}
int main( int argc, char** argv )
{
    RunTime();

    String empty;
    BeginCommandArguments;
    CommandArgument_String_Doc(DIR,"HiC-to-line data dir.");
    CommandArgument_String_OrDefault_Doc(LOCUS,empty,"Genomic locus.");
    CommandArgument_Int_OrDefault_Doc(LINE,-1,"Starting line.");
    CommandArgument_Int_OrDefault(MIN_LINE_LEN,2000);
    EndCommandArguments;

    if ( LINE == -1 && LOCUS == empty )
        FatalErr("You must supply either a LINE arg or a LOCUS arg.");

    LineInfoVec lineInfoVec;
    BinaryReader::readFile(DIR+'/'+LineInfo::gFileName,&lineInfoVec);

    if ( LINE != -1 )
    {
        if ( unsigned(LINE) >= lineInfoVec.size() )
            FatalErr("Invalid line -- there are only "
                        << lineInfoVec.size() << " lines.");
        if ( !isTarget(lineInfoVec[LINE],MIN_LINE_LEN) )
            FatalErr("Line " << LINE << " is not a valid target.");
    }

    VecHiCHitRateVec hitRates(DIR+'/'+HiCHitRate::gFileName);
    ForceAssertEq(lineInfoVec.size(),hitRates.size());

    vec<int> targets;
    setTargets(lineInfoVec,MIN_LINE_LEN,&targets);

    // TODO: use sorted lineinfo and locus inputs to determine this
    vec<int> lines;
    size_t const MAX_ROWS = 30ul;
    lines.reserve(MAX_ROWS);

    auto beg = targets.begin();
    auto end = targets.end();
    auto itr = end;
    if ( LINE != -1 )
        itr = std::find(beg,end,LINE);
    else if ( LOCUS != empty )
    {
        std::istringstream iss(LOCUS);
        int chrId, offset;
        iss >> chrId; iss.get(); iss >> offset;
        if ( !iss )
            FatalErr("Couldn't turn " << LOCUS << " into a genomic locus.");
        GenomeLoc loc(chrId,offset);
        itr = std::lower_bound(beg,end,loc,
                                [&lineInfoVec]( int id, GenomeLoc const& loc )
                                { return lineInfoVec[id].mLoc < loc; });
    }
    ForceAssert(itr != end);

    int chrId = lineInfoVec[*itr].mLoc.getChrId();
    size_t nnn = MAX_ROWS/2;
    beg = itr;
    while ( nnn-- && beg != targets.begin() )
        if ( lineInfoVec[*--beg].mLoc.getChrId() != chrId )
        {
            ++beg;
            break;
        }
    nnn = MAX_ROWS - MAX_ROWS/2;
    end = itr;
    while ( nnn-- && ++end != targets.end() )
        if ( lineInfoVec[*end].mLoc.getChrId() != chrId )
            break;
    while ( beg != end )
    {
        lines.push_back(*beg);
        ++beg;
    }

    std::cout << std::fixed << std::setprecision(1);
    size_t const MAX_COLS = 10ul;
    size_t excess = lines.size()-std::min(lines.size(),MAX_COLS);
    auto colsBeg = lines.begin()+(excess/2);
    auto colsEnd = lines.end()-((excess+1)/2);
    std::cout << "       ";
    for ( auto colsItr = colsBeg; colsItr != colsEnd; ++colsItr )
        std::cout << std::setw(8) << *colsItr;
    std::cout << '\n';
    std::cout << "       ";
    for ( auto colsItr = colsBeg; colsItr != colsEnd; ++colsItr )
        std::cout << "========";
    std::cout << '\n';
    for ( int lineId : lines )
    {
        if ( lineId == *colsBeg )
        {
            std::cout << "       ";
            for ( auto colsItr = colsBeg; colsItr != colsEnd; ++colsItr )
                std::cout << "--------";
            std::cout << '\n';
        }
        std::cout << std::setw(7) << lineId;
        HiCHitRateVec const& rateVec = hitRates[lineId];
        std::unordered_map<int,float> map(2*rateVec.size());
        for ( HiCHitRate const& hitRate : rateVec )
            map[hitRate.mLineId] = hitRate.mHitsPerKb;
        auto end = map.end();
        for ( auto colsItr = colsBeg; colsItr != colsEnd; ++colsItr )
        {
            int lineId2 = *colsItr;
            auto itr = map.find(lineId2);
            std::cout << std::setw(8);
            if ( itr == end )
                std::cout << '0';
            else
                std::cout << itr->second;
        }
        LineInfo const& info = lineInfoVec[lineId];
        std::cout << std::setw(10) << info.mFirstEdgeId << ".."
                  << std::setw(10) << std::left << info.mLastEdgeId
                  << std::right << std::setw(9) << info.mLineLen
                  << ' ' << info.mLoc << '\n';
        if ( lineId == colsEnd[-1] )
        {
            std::cout << "       ";
            for ( auto colsItr = colsBeg; colsItr != colsEnd; ++colsItr )
                std::cout << "--------";
            std::cout << '\n';
        }
    }

    std::cout << '\n';
    for ( auto colsItr = colsBeg; colsItr != colsEnd; ++colsItr )
    {
        int lineId = *colsItr;
        std::cout << std::setw(8) << lineId;
        HiCHitRateVec const& rateVec = hitRates[lineId];
        int nCols = 0;
        for ( auto itr=rateVec.begin(),end=rateVec.end(); itr != end; ++itr )
        {
            int lineId2 = itr->mLineId;
            if ( std::find(lines.begin(),lines.end(),lineId2) == lines.end() &&
                    lineInfoVec[lineId2].mLineLen >= MIN_LINE_LEN )
            {
                std::cout << std::setw(10) << itr->mLineId
                          << std::setw(7) << itr->mHitsPerKb;
                if ( ++nCols == 5 ) break;
            }
        }
        std::cout << '\n';
    }
}
