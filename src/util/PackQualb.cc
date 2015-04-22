///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * PackQualb.cc
 *
 *  Created on: Aug 12, 2014
 *      Author: tsharpe
 */

#include "MainTools.h"
#include "feudal/PQVec.h"
#include "feudal/VirtualMasterVec.h"
#include <iostream>

int main( int argc, char** argv )
{
    RunTime();
    String empty;
    BeginCommandArguments;
    CommandArgument_String(QUALB);
    CommandArgument_String_OrDefault(QUALP,empty);
    EndCommandArguments;

    if ( QUALP == empty )
        QUALP = QUALB.ReplaceExtension(".qualb",".qualp");

    vecqvec vqv(QUALB);
    VecPQVec vpqv;
    convertAssignParallel(vqv.begin(),vqv.end(),vpqv);
    vpqv.WriteAll(QUALP);
    std::cout << "Done." << std::endl;
}
