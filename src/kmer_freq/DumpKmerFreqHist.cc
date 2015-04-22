/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Dump to stdout a histogram of the kmer frequencies in TABLE.

#include "MainTools.h"
#include "kmers/KmerShape.h"
#include "kmer_freq/KmerFrequencyTable.h"

int main( int argc, char** argv ) {

  RunTime();

  BeginCommandArguments;
  CommandArgument_String(TABLE);
  CommandArgument_String(K);
  CommandArgument_Bool_OrDefault(PRINT_OVERFLOW,True);
  EndCommandArguments;

  KmerFrequencyTable table( KmerShapeId( K ), TABLE );

  vec<longlong> hist;
  table.GetHistogram(hist);

  copy( hist.begin() + 1, hist.end() - ( PRINT_OVERFLOW ? 0 : 1 ),
        ostream_iterator<longlong>( cout, "\n" ) );
}
