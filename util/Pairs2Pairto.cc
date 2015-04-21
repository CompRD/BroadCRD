///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Converts a .pairs file to a .pairto file.";

#include "MainTools.h"
#include "PairsManager.h"
#include "ReadPairing.h"

int main( int argc, char *argv[] )
{  
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc( IN, 
    "Input filename (.pairs file format)");
  CommandArgument_String_OrDefault_Doc( OUT, "",
    "Output filename (.pairto file format)" );
  EndCommandArguments;

  String in_head = IN.SafeBefore(".pairs");

  String pairs_file = in_head + ".pairs";
  String pairto_head = (OUT.empty() ? in_head : OUT.SafeBefore(".pairto") );

  cout << Date() << " Reading pairs file: " << pairs_file << endl;
  PairsManager pairs( pairs_file );

  cout << Date() << " Converting pairing format" << endl;
  vec<read_pairing> pairto = pairs.convert_to_read_pairings();;

  cout << Date() << " Writing pairto files: " << pairto_head << ".pairto*" << endl;
  WritePairs( pairto, pairs.nReads(), pairto_head );

  cout << Date() << " Done!" << endl;
}
