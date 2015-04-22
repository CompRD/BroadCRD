/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FirstLookupTable: call FirstLookup.

#include "Basevector.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "lookup/FirstLookup.h"

int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(SEQS);
  CommandArgument_String_Abbr(LOOKUP_TABLE, L);
  CommandArgument_String(OUT);
  CommandArgument_String_OrDefault( NUM_THREADS, 1 );
  EndCommandArguments;
  
  vecbasevector seqs(SEQS);
  vec<look_align> aligns;
  FirstLookup( seqs, LOOKUP_TABLE, aligns );
  Ofstream( out, OUT );
  for ( int i = 0; i < aligns.isize( ); i++ )
    aligns[i].PrintParseable(out);
}
