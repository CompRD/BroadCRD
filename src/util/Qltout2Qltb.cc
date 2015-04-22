/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "FastIfstream.h"

/*
 * Qltout2Qltb 
 *
 * Pack alignments into binary format
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String_Doc( IN, "Name of the input qltout file" );
  CommandArgument_String_OrDefault_Doc( OUT, "" , "Name of the output qltb file. Defaults to IN\\.{qltout}.qltb" );
  CommandArgument_Bool_OrDefault_Doc( SORT, True, "Sort alignments by query id then errors" );
  EndCommandArguments;

  if ( OUT.empty() ) {
    OUT = IN.SafeBefore(".qltout")+".qltb";
  }

  if ( ! IsRegularFile(IN) ) {
    if ( IsRegularFile(IN+".qltout") ) IN+=".qltout";
    else {
      cout << "Can not locate input data: neither " << IN << " nor " << (IN+".qltout") << " file exists" << endl;
    }
  }

  String line;
  fast_ifstream in(IN);

  vec<look_align> aligns;

  look_align la;
  ulonglong count = 0;

  const String tag("QUERY");
  bool good = false;
  while(1) {
    good = getline_if_match( in, line , tag);
    if ( in.fail( ) ) break;
    if ( good ) {
      la.ReadParseable(line);

      // insert remapped alignment; collector should recognize the attempt to insert the 
      // same align again and should not consider this as inserting a "non-unique" alignment
      aligns.push_back(la);
      count++;
    }
  }
  if ( SORT ) sort(aligns.begin(), aligns.end(), order_lookalign_QueryErrors());
  BinaryWriter::writeFile(OUT,aligns);
  cout << count << " aligns written to " << OUT << endl;
}
