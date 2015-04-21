///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "math/Functions.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/BubbleFreeN50.h"

int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(INPUT, "input.{hbv,shbv}");
     CommandArgument_Int_OrDefault_Doc(MIN_LEN, 1000,
          "minimum edge length in bases to be included in N50 calculation");
     EndCommandArguments;


     // load (Supported-)HyperBasevector
     HyperBasevector hbv;
     if ( INPUT.EndsWith(".hbv") )
          BinaryReader::readFile( INPUT, &hbv );
     else if ( INPUT.EndsWith(".shbv") )  {
          SupportedHyperBasevector shbv;
          BinaryReader::readFile( INPUT, &shbv );
          hbv = shbv;
     } else
          FatalErr("file must end in .hbv or .shbv");

     BubbleFreeN50 calc(hbv, MIN_LEN );

     cout << "before popping bubbles we have " << calc.PreNEdges()
             << " edges" << endl;
     cout << "bubble-full edge N50 is " << calc.PreN50() << endl;
     cout << "after popping bubbles we have " << calc.PostNedges()
             << " edges" << endl;
     cout << "bubble-free edge N50 is " << calc.PostN50() << endl;
}
