///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>
#include "MainTools.h"
#include "paths/long/hic/HiCDefs.h"
#include "VecUtilities.h"
#include <iosfwd>
#include "system/System.h"
#include <queue>
#include <unordered_set>
#include "paths/long/hic/TestTree.h"
#include "paths/long/hic/HiCExperiment.h"
#include "paths/long/hic/Scores.h"
#include "graph/Digraph.h"
#include "paths/long/hic/LineGraph.h"
#include "paths/long/hic/CanonicalTestOrder.h"
#include "paths/long/hic/LineMapReader.h"
#include "feudal/BinaryStream.h"
#include "paths/long/hic/ScaffoldPile.h"
#include "paths/long/hic/NeighborScaff.h"
#include <array>


namespace {
template < typename T >
void printCols( size_t ncols, T const& s )
{
    std::ostringstream out;
    out << s;
    while ( out.str().size() < ncols ) {
        out << " ";
    }
    cout << out.str();
}

};



int main( int argc, char* argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(ADIR,"assembly dir.");
    CommandArgument_String_Doc(HICDIR,"hic.* files from CalcHiCLineHits");
    CommandArgument_String_Doc(LINEPAIRS,"Hi-C LINE-pairing file.");
    CommandArgument_String_Doc(SEPDIST,"hic.*.sepdist file");
    CommandArgument_Int(DEPTH);
    CommandArgument_Int_OrDefault_Doc(START_ID, -1, "starting line ID");
    CommandArgument_Int_OrDefault_Doc(NUM_LINES, -1, "starting line ID");
    CommandArgument_Bool_OrDefault(VERBOSE,True);
    CommandArgument_Bool_OrDefault(FILTERSELF,False);
    CommandArgument_Int_OrDefault_Doc(MINLINESIZE,6000,"ignore lines smaller than this");
    EndCommandArguments;

    String mapFile = ADIR + "/" + "a.lines.map";
    cout << Date() << ": reading line map from " << mapFile << endl;
    auto lineMap = LineMapReader(mapFile).LineLoci();
    cout << Date() << ": read " << lineMap.size() << " line alignments." << endl;

    HiCExperiment hic(
            ADIR + "/a.hbx",          // these are currently being loaded for no reason...
            ADIR + "/a.inv",
            HICDIR + "/" + EdgeIdToLineIdMap::gFileName,
            HICDIR + "/" + LineInfo::gFileName,
            ADIR + "/" + "a.lines",
            HICDIR + "/hic.hitrates",
            SEPDIST
            );

    // stash away line lengths
    // note that this could just be a vec... a relic right now that needs
    // to be done away with.
    unordered_map<int,int> lineLengths;
    for ( int lineId = 0; lineId < hic.mLineInfoVec.isize(); ++lineId )
        lineLengths[lineId] = hic.mLineInfoVec[lineId].mLineLen;


    cout << Date() << ": creating line involution" << endl;
    vec<int> lineInv;
    createLineInvolution(hic.mLineInfoVec, lineInv);

    cout << Date() << ": loading linePairs from " << LINEPAIRS << endl;
    LineLocPairVec linePairs;
    BinaryReader::readFile(LINEPAIRS, &linePairs);
    cout << Date() << ": " << linePairs.size() << " pairs loaded." << endl;

    if ( FILTERSELF ) {
        linePairs.erase( std::remove_if( linePairs.begin(), linePairs.end(),
                [](LineLocPair const& llp) -> bool {
                return (llp.ll1().getLineId() == llp.ll2().getLineId());
        } ), linePairs.end() );
        linePairs.shrink_to_fit();
        cout << Date() << ": " << linePairs.size()  << " after filtering." << endl;
    }

    ScaffoldPile pile(lineInv);
    for ( auto itr = lineMap.begin(); itr != lineMap.end(); ++itr ) {
        if ( itr->ref_gap ) continue;
        if ( itr->len < MINLINESIZE ) continue;
        pile.add({itr->lineId});
    }
//    for ( int lineId = 0; lineId < hic.mLineInfoVec.isize(); ++lineId ) {
//        if ( !hic.mLineInfoVec[lineId].mIsRC && hic.mLineInfoVec[lineInv[lineId]].mIsRC ) {
//            pile.add({lineId});
//        }
//    }

    vec<unsigned> worst_hist(DEPTH+1,0u);
    vec<unsigned> best_hist(DEPTH+1,0u);
    int last1_line=-1;
    String last1_chr;
    int last0_line=-1;
    String last0_chr;
    unsigned not_found_left = 0u;
    unsigned not_found_right = 0u;
    unsigned not_found = 0u;
    auto itr = lineMap.begin();
    for ( ; itr != lineMap.end() ; ++itr ) {
        auto mapLine = *itr;

        // skip if ref gap
        if ( mapLine.ref_gap ) continue;
        // skip if < MINLINESIZE
        if ( mapLine.len < MINLINESIZE ) continue;

        cout << itr->chr << ":" << itr->start << "-" << itr->stop << endl;

        if ( last1_line != -1 && last0_line != -1 && last1_chr == last0_chr && last1_chr == mapLine.chr ) {
            vec<int> testStarts;
            neighborScaff( last0_line, pile, hic, DEPTH, MINLINESIZE, testStarts, false );
            ForceAssertLe(testStarts.isize(), DEPTH+1);

            // find index of left in testStarts
            int found_left = -1;
            for ( size_t i = 0; i < testStarts.size(); ++i )
                if ( testStarts[i] == last1_line ) { found_left = i; break; }

            // find index or right in testStarts
            int found_right = -1;
            for ( size_t i = 0; i < testStarts.size(); ++i )
                if ( testStarts[i] == mapLine.lineId ) { found_right = i; break; }

            if ( found_left != -1 || found_right != -1 ) {
                int best, worst;
                if ( found_left != -1 ) {
                    best = worst = found_left;  // assume only left found
                    if ( found_right != -1 ) {  // okay, but we found right
                        if ( found_right < found_left ) // see which is better
                            best = found_right;
                        else
                            worst = found_right;
                    }
                } else {        // only have right found
                    best = worst = found_right;
                }

                if ( found_left == -1 ) not_found_left++;
                if ( found_right == -1 ) not_found_right++;
                if ( found_left == -1 && found_right == -1 ) not_found++;

                best_hist[best]++;
                worst_hist[worst]++;
            }

        }

        last1_chr = last0_chr;
        last1_line = last0_line;
        last0_chr = mapLine.chr;
        last0_line = mapLine.lineId;
    }

    cout << "best hist: " << printSeq(best_hist) << endl;
    cout << "worst hist: " << printSeq(worst_hist) << endl;
    PRINT3(not_found_left, not_found_right, not_found);


    return 0;
}



