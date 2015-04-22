/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Coverage_Solexa.  Plot coverage distribution for a Solexa data set.  This only
// uses perfect reads.  Pileups are not filtered out.

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "bias/CoveragePlot.h"
#include "lookup/LookAlign.h"
#include "solexa/FourBase.h"
#include "solexa/SolexaTools.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HEAD);
     CommandArgument_Int(WINDOW);
     CommandArgument_Int_OrDefault(GRANULARITY, 1);
     CommandArgument_Int_OrDefault(TIG, 0);
     CommandArgument_Int_OrDefault(START, 0);
     CommandArgument_Int_OrDefault(STOP, -1);
     CommandArgument_Double_OrDefault(COVERAGE_DIVIDER, 1.0);
     CommandArgument_String(OUT);
     CommandArgument_String_OrDefault(QLTOUT_SUFFIX, "");
     EndCommandArguments;

     // Set up data structures.

     vecbasevector reads, readsrc, ref, rref;
     VecFourBaseVec I;
     vec<look_align> aligns;
     vec< vec<int> > aligns_index; 
     LoadSolexaData( HEAD, reads, readsrc, I );
     String aligns_file = HEAD + ".qltout";
     if ( QLTOUT_SUFFIX != "" ) aligns_file += "." + QLTOUT_SUFFIX;
     String ref_file = HEAD + ".reference.fasta";
     FetchReads( ref, 0, ref_file );
     rref = ref;
     ReverseComplement(rref);
     LoadLookAligns( aligns_file, aligns, aligns_index, reads.size( ) );
     LoadSolexaData( HEAD, reads, readsrc, I );
     int nreads = reads.size( );
     ForceAssertLt( static_cast<size_t>(TIG), ref.size( ) );

     // Find best alignment for each read.

     vec<int> best;
     GetBestAligns( reads, aligns, aligns_index, best );

     // Select perfect reads.

     vec<look_align> xaligns;
     for ( int id = 0; id < nreads; id++ )
     {    if ( best[id] < 0 || aligns_index[id].empty( ) ) continue;
          const look_align& la = aligns[ aligns_index[id][ best[id] ] ];
          if ( la.Errors( ) > 0 ) continue;
          xaligns.push_back(la);    }

     // Make the plot.

     CoveragePlot( reads, xaligns, ref, TIG, START, STOP, WINDOW, GRANULARITY, 
          COVERAGE_DIVIDER, OUT );    }
