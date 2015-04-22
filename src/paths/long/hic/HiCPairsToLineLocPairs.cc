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
#include "paths/long/hic/HiCExperiment.h"




int main( int argc, char* argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(ADIR,"assembly dir.");
    CommandArgument_String_Doc(HICDIR,"hic.* files from CalcHiCLineHits");
    CommandArgument_String_Doc(HICPAIRS,"Hi-C edge-pairing file.");
    CommandArgument_String_OrDefault_Doc(OUTHEAD, "",
            "output head -- default is input with .aligns replaced by .line_aligns" );
    EndCommandArguments;

    HiCExperiment hic(
            ADIR + "/a.hbx",
            ADIR + "/a.inv",
            HICDIR + "/" + EdgeIdToLineIdMap::gFileName,
            HICDIR + "/" + LineInfo::gFileName,
            ADIR + "/" + "a.lines"
            );

    cout << Date() << ": reading HiC pairs from "  << HICPAIRS << endl;
    HICVec hiCPairs;
    getHiCPairs( HICPAIRS, hic.mEdgeToLineMap, &hiCPairs,
             false /* don't skip intra-line pairs */ );

    cout << Date() << ": done reading " << hiCPairs.size() << " pairs"<< endl;

    cout << Date() << ": calculating edge pairs -> line pairs" << endl;
    LineOffsetCalc offCalc( hic.mHBX, hic.mLines, hic.mInv,
            hic.mEdgeToLineMap, hic.mRuler, LineOffsetCalc::THREAD_SAFE );

    cout << Date() << ": preloading cache of line lengths and (line,edge)->segment map" << endl;
    for ( size_t testLine = 0; testLine < hic.mLines.size(); ++testLine )
        offCalc.cachePreload(testLine);

    cout << Date() << ": computing line pairs" << endl;
    LineLocPairVec linePairs = fromHICVec(hiCPairs, offCalc);

    if ( OUTHEAD == "" ) OUTHEAD = HICPAIRS.ReplaceExtension(".aligns", "");
    String outFileName = OUTHEAD + ".line_aligns";
    ForceAssertNe(outFileName, HICPAIRS);

    cout << Date() << ": writing results to " << outFileName << endl;
    BinaryWriter::writeFile(outFileName, linePairs);


    return 0;
}

