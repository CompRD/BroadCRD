///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Dec 13, 2013 - <crdhelp@broadinstitute.org>
//

#include "paths/long/ReadPath.h"
#include "MainTools.h"
#include "feudal/VirtualMasterVec.h"



int main( int argc, char* argv[] )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc( IN_PATHS, "ReadPathsVec file");
    CommandArgument_Int_OrDefault_Doc( ID, -1, "optional read id" );
    EndCommandArguments;

//    VirtualMasterVec<ReadPath> vpaths( IN_PATHS );
    ReadPathVec vpaths;
    vpaths.ReadAll(IN_PATHS);
    cout << "VirtualMasterVec of " << vpaths.size() << " elements" << endl;

    if ( ID >= 0 ) {
	cout << ID << ": " << vpaths[ID] << endl;
    } else {
	for ( size_t i = 0; i < vpaths.size(); ++i )
	    cout << i << ": " << vpaths[i] << endl;
    }

    return 0;
}
