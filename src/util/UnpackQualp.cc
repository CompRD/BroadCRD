///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * UnpackQualp.cc
 *
 *  Created on: Aug 29, 2014
 *      Author: tsharpe
 */

#include "MainTools.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/PQVec.h"
#include "feudal/VirtualMasterVec.h"
#include <iostream>

int main( int argc, char** argv )
{
    RunTime();
    String empty;
    BeginCommandArguments;
    CommandArgument_String(QUALP);
    CommandArgument_String_OrDefault(QUALB,empty);
    EndCommandArguments;

    if ( QUALB == empty )
        QUALB = QUALP.ReplaceExtension(".qualp",".qualb");

    VirtualMasterVec<PQVec> vpqv(QUALP);
    IncrementalWriter<qvec> vqv(QUALB,vpqv.size());
    qvec qv;
    for ( PQVec const& pqvec : vpqv )
    {
        pqvec.unpack(&qv);
        vqv.add(qv);
    }
    vqv.close();
    std::cout << "Done." << std::endl;
}
