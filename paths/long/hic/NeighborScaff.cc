///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Oct 8, 2014 - <crdhelp@broadinstitute.org>
//


using namespace std;
#include "paths/long/hic/HiCDefs.h"
#include "paths/long/hic/HiCExperiment.h"
#include "paths/long/hic/ScaffoldPile.h"


void neighborScaff( int scaff_seed,
        ScaffoldPile const& scaffPile,
        HiCExperiment const& hic,
        int DEPTH,
        int MINLINESIZE,
        vec<int>& testStarts,
        size_t head,
        Bool verb)
{
    testStarts.clear();
    std::unordered_map<int,float> testLineMap;

    // this is canonicalized w.r.t. whatever canonicalization is used for lines
    int cscaff_seed = canonicalLine(scaff_seed, hic.mLineInfoVec);
    // this is canonicalized w.r.t. whatever is actually present in the scaffold pile
    // after whatever flipping has happened.

    int scaff_start = scaffPile.start(cscaff_seed);
    vec<int> scaffLines = scaffPile(cscaff_seed);

    if ( head == 0 ) head = scaffLines.size();
    vec<int> cScaffLines;
    for ( size_t i = 0; i < head && i < scaffLines.size(); ++i ) {
        cScaffLines.push_back( canonicalLine(scaffLines[i], hic.mLineInfoVec) );
    }
    // walk from the end for 'head' values, but ensure that you don't hit the
    // 'head' values that we just added above
    for ( size_t i = scaffLines.size(); scaffLines.size() - i < head
            && i > head && i > 0; --i ) {
        cScaffLines.push_back( canonicalLine(scaffLines[i-1], hic.mLineInfoVec) );
    }

    for ( auto const l : scaffLines )
        cScaffLines.push_back( canonicalLine(l,hic.mLineInfoVec));

    for ( auto const cseed : cScaffLines ) {
        auto const& rates = hic.mHitRates[cseed];
        if ( rates.size() == 0 ) continue;
        size_t si = 0;
        for ( ; si < rates.size(); ++si ) {
            auto len = hic.mLineInfoVec[rates[si].mLineId].mLineLen;
            if ( len >= MINLINESIZE ) break;
        }
        if ( si == rates.size() ) continue;
        auto rateThresh = rates[si].mHitsPerKb * 0.10;
        if ( verb ) {
            cout << "seed=" << cseed << ", hits=";
            for ( size_t i = si; i < rates.size() &&
                rates[i].mHitsPerKb >= rateThresh; i++ )
                cout << rates[i].mLineId << " (" << rates[i].mHitsPerKb << ") ";
            cout << endl;
            PRINT3(cseed, rates.size(), rateThresh );
        }
        for ( auto const& hit : rates ) {
            auto lineId = hit.mLineId;
            if ( hit.mHitsPerKb < rateThresh ) continue;
            if ( Member( cScaffLines, lineId )) continue;
            auto lineLen = hic.mLineInfoVec[lineId].mLineLen;
            if ( lineLen < MINLINESIZE ) continue;
            try {
                if ( testLineMap.at(lineId) > hit.mHitsPerKb )
                    testLineMap[lineId] = hit.mHitsPerKb;
            } catch ( std::out_of_range ) {
                testLineMap[lineId] = hit.mHitsPerKb;
            }
        }
    }

    if ( verb ) {
        for ( auto itr = testLineMap.cbegin();
                itr != testLineMap.cend(); ++itr ) {
            cout << "line " << itr->first << " max hits " << itr->second <<
                    endl;
        }

    }

    vec<float> allTestRates;
    vec<int> allTestLines;
    for ( auto itr = testLineMap.begin(); itr != testLineMap.end(); ++itr ) {
        allTestLines.push_back( itr->first );
        allTestRates.push_back( itr->second );
    }

    ReverseSortSync(allTestRates, allTestLines);

    set<int> seen;
    seen.insert(scaff_start);
    testStarts.push_back(scaff_start);
    for ( auto line : allTestLines ) {
        if (!scaffPile.has_key_or_inv(line)) continue;
        int start_line = scaffPile.start(line);
        if ( seen.find(start_line) == seen.end() ) {
            seen.insert(start_line);
            testStarts.push_back(start_line);
        }
    }

    if ( testStarts.isize() > DEPTH+1 ) testStarts.resize( DEPTH+1 );

}
