///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * TestPackedQVec.cc
 *
 *  Created on: Aug 18, 2014
 *      Author: tsharpe
 */

#include "MainTools.h"
#include "feudal/PQVec.h"
#include <iostream>

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String(QUALB);
    EndCommandArguments;

    vecqvec vqv(QUALB);
    VecPQVec vPQv;
    convertAssignParallel(vqv.begin(),vqv.end(),vPQv);
    auto itr = vqv.begin();
    for ( qvec qv : vPQv )
    {
        ForceAssert(qv == *itr);
        ++itr;
    }
    std::cout << "Done." << std::endl;
}
