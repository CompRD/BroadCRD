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
    CommandArgument_Bool_OrDefault(VERBOSE,False);
    CommandArgument_Bool_OrDefault(FILTERSELF,False);
    CommandArgument_Int_OrDefault_Doc(MINLINESIZE,6000,"ignore lines smaller than this");
    CommandArgument_IntSet_Doc(SEED, "seed scaffold, as an ordered set");
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
    for ( int lineId = 0; lineId < hic.mLineInfoVec.isize(); ++lineId ) {
        if ( !hic.mLineInfoVec[lineId].mIsRC && hic.mLineInfoVec[hic.mLineInfoVec[lineId].mLineIdRC].mIsRC ) {
//            PRINT3(lineId,lineInv[lineId],hic.mLineInfoVec[lineId].mLineIdRC);
//            cout << "adding " << lineId << endl;
            pile.add({lineId});
        }
    }

    for ( auto const lineId : SEED )
        if ( hic.mLineInfoVec[lineId].mIsRC )
            FatalErr("we're not set up for RC lines yet");

    ForceAssertGt(SEED.size(), 0u);
    pile.merge(SEED);
#if 0
    for ( size_t i = 1; i < SEED.size(); ++i ) {
        cout << "about to join seed " << SEED[0] << " to " << SEED[i] << endl;
        pile.join(SEED[0],SEED[i]);
    }
#endif

    vec<int> testStarts;
    neighborScaff( SEED[0], pile, hic, DEPTH, MINLINESIZE, testStarts, 0, VERBOSE );

    cout << "neighbors for " << pile.scaffold(SEED[0]) << ":" << endl;
    cout << printSeq(testStarts) << endl;


    return 0;
}



