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

    String mapFile = ADIR + "/" + "a.s.lines.map";
    cout << Date() << ": reading line map from " << mapFile << endl;
    auto lineMap = LineMapReader(mapFile).LineLoci();
    cout << Date() << ": read " << lineMap.size() << " line alignments." << endl;

    HiCExperiment hic(
            ADIR + "/a.s.hbx",          // these are currently being loaded for no reason...
            ADIR + "/a.s.inv",
            HICDIR + "/" + EdgeIdToLineIdMap::gFileName,
            HICDIR + "/" + LineInfo::gFileName,
            ADIR + "/" + "a.s.lines",
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




    // initialize last with chr, line
    String last_chr = "";
    int last_line = -1;
    int last_len = 0;
    int last_stop = 0;

    // for each line
    auto itr = lineMap.begin();
    if ( START_ID != -1 ) {
        while ( itr != lineMap.end() && itr->lineId != START_ID ) itr++;
        if ( itr == lineMap.end() ) {
            FatalErr("specified START_ID was not found");
        }
    }

    vec<unsigned> hist(4,0u);
    for ( int num_lines = NUM_LINES; itr != lineMap.end() && num_lines--; ++itr ) {
        auto mapLine = *itr;

        // skip if ref gap
        if ( mapLine.ref_gap ) continue;
        // skip if < MINLINESIZE
        if ( mapLine.len < MINLINESIZE ) continue;
        // if last_chr != -1, test against last
        auto const metric = Scores::LOG_LIK_A;
        if ( last_chr == mapLine.chr ) {
            auto scores = testLinePair( {last_line,false}, {mapLine.lineId,false},
                    lineLengths, linePairs, hic.mSepDist );
            vec<int> idx(scores.size(), vec<int>::IDENTITY );
            SortSync( scores, idx, [metric]( Scores const& s1, Scores const& s2 ) {
                return s1(metric) < s2(metric);
            });
            double ratio = 0.;
            winningPairOrderOrient( scores, ratio, metric );
            double ratio2 = winningPairRatio(scores, 0, metric);

            // zero is always the correct configuration; find out where it
            // ended up in the list
            size_t win = 0;
            if ( idx[1] == 0 ) win = 1;
            else if ( idx[2] == 0 ) win = 2;
            else if ( idx[3] == 0 ) win = 3;
            else ForceAssertEq(idx[0],0);
            hist[win]++;
            ostringstream out;
            out << mapLine.chr << ":" << mapLine.start << "-" << mapLine.stop;
            printCols(30,out.str());
            out.str("   ");
            out << last_line << "," << mapLine.lineId;
            printCols(20,out.str());
            cout << "   " << printSeq(idx) << "   \t";
            printCols(10,ToString(last_len));
            cout << " ";
            printCols(10,ToString(mapLine.len));
            cout << " ";
            printCols(10,ToString(mapLine.start - last_stop));
            cout << "   ";
            if ( win != 0 ) cout << "(";
            printCols(10, ToString(ratio));
            if ( win != 0 ) cout << ")";
            cout << "   ";
            if ( win != 0 ) cout << "(";
            printCols(10, ToString(ratio2));
            if ( win != 0 ) cout << ")";

            cout << endl;
        }
        // set last to cur
        last_chr = mapLine.chr;
        last_line = mapLine.lineId;
        last_len = mapLine.len;
        last_stop = mapLine.stop;
    }

    cout << "histogram of where the correct config ended up: " << printSeq(hist) << endl;

    return 0;
}



