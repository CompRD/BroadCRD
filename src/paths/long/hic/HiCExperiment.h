///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Jul 28, 2014
//

#ifndef __HICEXPERIMENT_H
#define __HICEXPERIMENT_H

// utility class to package a bunch of data that we frequently use
struct HiCExperiment {
    HyperBasevectorX mHBX;
    HBVXKmerRuler mRuler;
    vec<int> mInv;
    EdgeIdToLineIdMap mEdgeToLineMap;
    LineInfoVec mLineInfoVec;
    LineVec mLines;
    VecHiCHitRateVec mHitRates;
    SepHistogram mSepDist;

    HiCExperiment() = delete;
    HiCExperiment( String const& hbxFile,
            String const& invFile = "",
            String const& mapFile = "",
            String const& lineInfoFile = "",
            String const& linesFile = "",
            String const& hitsFile = "",
            String const& sepDistFile = ""
            ) : mHBX(loadHBX(hbxFile)), mRuler(mHBX) {
        if ( invFile != "" ) loadInvolution(invFile);
        if ( mapFile != "" ) loadEdgeToLineMap(mapFile);
        if ( lineInfoFile != "" ) loadLineInfo(lineInfoFile);
        if ( linesFile != "" ) loadLines(linesFile);
        if ( hitsFile != "" ) loadHitRates(hitsFile);
        if ( sepDistFile != "" ) loadSepHist(sepDistFile);
    }

    HyperBasevectorX loadHBX(String const& hbxFile) {
        cout << Date() << ": reading HyperBasevectorX(treme)" << endl;
        HyperBasevectorX hbx;
        BinaryReader::readFile( hbxFile, &hbx );
        return hbx;
    }

    void loadInvolution(String const& invFile) {
        cout << Date() << ": reading involution" << endl;
        BinaryReader::readFile(invFile, &mInv);
    }

    void loadEdgeToLineMap(String const& mapFile ) {
        cout << Date() << ": reading edge to line id map" << endl;
        BinaryReader::readFile(mapFile, &mEdgeToLineMap );
    }

    void loadLineInfo(String const& lineInfoFile ) {
        cout << Date() << ": reading lineInfo" << endl;
        BinaryReader::readFile( lineInfoFile, &mLineInfoVec );
    }

    void loadLines(String const& linesFile) {
        cout << Date() << ": reading lines"  << endl;
        LineVec lines;
        BinaryReader::readFile(linesFile, &mLines);
    }

    void loadHitRates(String const& hitsFile) {
        cout << Date() << ": reading hitrates" << endl;
        mHitRates = VecHiCHitRateVec(hitsFile);
    }

    void loadSepHist(String const& sepDistFile, float smooth = 2.) {
        cout << Date() << ": reading distribution" << endl;
        BinaryReader::readFile(sepDistFile, &mSepDist);
        mSepDist.Smooth(smooth);
        mSepDist.FixAnalyticSum();
    }

};

#endif
