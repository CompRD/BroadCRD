// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// TruncateFasta: for each fasta record in a file, delete all bases
// after base N.  Incidentally and pointlessly, this code changes
// ambiguous bases to random bases.  If RENUMBER is set (default), the
// names from the fasta file are not read and are replaced with
// sequence numbers in the output.  Otherwise names are read from the
// fasta and used in the output.

#include "Basevector.h"
#include "MainTools.h"
#include "FetchReads.h"

int main( int argc, char *argv[] )
{
  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_UnsignedInt(N);
  CommandArgument_Bool_OrDefault(RENUMBER, True);
  EndCommandArguments;

  vecbasevector b;
  vecString names; // used if !RENUMBER
  if (RENUMBER) {
    FetchReads( b, 0, IN );
  } else {
    vecqualvector q; // unused
    FetchReads(b, q, &names, 0, IN, 0, 0, cout, True);
  }

  for ( size_t i = 0; i < b.size( ); i++ ) {
    if ( b[i].size( ) > N )
      b[i].resize(N);
    if (RENUMBER)
      b[i].Print( cout, i );
    else
      b[i].Print( cout, names[i] );
  }
}
