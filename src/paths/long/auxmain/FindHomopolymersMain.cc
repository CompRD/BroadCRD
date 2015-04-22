///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FindHomopolymers.cc
 *
 *  Created on: Mar 26, 2013
 *      Author: tsharpe
 */

#include "MainTools.h"
#include "IteratorRange.h"
#include "paths/long/FindHomopolymers.h"

int main( int argc, char** argv )
{
    RunTime( );
    BeginCommandArguments;
    CommandArgument_String(SHBV);
    CommandArgument_UnsignedInt_OrDefault(MIN_RUN_LEN,20u);
    EndCommandArguments;

    SupportedHyperBasevector shbv;
    BinaryReader::readFile(SHBV,&shbv);
    HomopolymerFinder hf(MIN_RUN_LEN);
    VecHomopolymerRun runs;
    hf.find(&runs,shbv);
    std::sort(runs.begin(),runs.end());
    std::cout << "All Runs:" << std::endl;
    std::cout << rangePrinter(runs.begin(),runs.end(),"\n") << std::endl;
    VecHomopolymerRunVec clusters;
    HomopolymerFinder::cluster(&clusters,runs);
    std::cout << "\n\nClusters:" << std::endl;
    for ( auto itr=clusters.begin(),end=clusters.end(); itr != end; ++itr )
    {
        std::cout << rangePrinter(itr->begin(),itr->end(),";") << std::endl;
    }
}
