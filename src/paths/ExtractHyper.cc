/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/HyperKmerPath.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HYPER_IN);
     CommandArgument_Int(K)
     CommandArgument_String_OrDefault(HYPER_OUT, "");
     CommandArgument_Int_Doc(COMPONENT,
       "Extract this component.");
     EndCommandArguments;

     // Load HyperKmerPath

     cout << Date() << "Loading HyperKmerPath" << endl;

     HyperKmerPath hyperIn;
     BinaryReader::readFile(HYPER_IN, &hyperIn);

    
     HyperKmerPath hyperOut;
     cout << Date( ) << ": Extracting HyperKmerPath component " << COMPONENT << endl;
     hyperOut = HyperKmerPath(hyperIn, COMPONENT);
 
     // Write Extracted Hyper

     if (HYPER_OUT == "")
       HYPER_OUT = HYPER_IN + ".out";

     cout << Date( ) << ": Writing HyperKmerPath component " << endl;
    
     BinaryWriter::writeFile( HYPER_OUT, hyperOut );

 }
