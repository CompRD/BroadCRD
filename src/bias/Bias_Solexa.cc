/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Bias_Solexa.
//
// HEAD is to be in ParseStringSet format.

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "bias/StartBias.h"
#include "bias/UniformBias.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "solexa/FourBase.h"
#include "solexa/SolexaTools.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HEAD);
     CommandArgument_Bool_OrDefault(BRIEF, True);
     CommandArgument_Bool_OrDefault(SUMMARY_FORM, False);
     CommandArgument_Int_OrDefault(MAX_ERRORS, 0);
     CommandArgument_String_OrDefault(QLTOUT_SUFFIX, "");
     CommandArgument_String_OrDefault(REF, "");
     CommandArgument_String_OrDefault(MARKED_REFERENCE_FILE, "");
     EndCommandArguments;

     // Set up data structures.

     vecbasevector bases, ref, rref;
     VecFourBaseVec I;
     vec<look_align> aligns;
     vec< vec<int> > aligns_index;
     String QLTOUT = "qltout";
     if ( QLTOUT_SUFFIX != "" ) QLTOUT += "." + QLTOUT_SUFFIX;
     LoadSolexaData( HEAD, QLTOUT, REF, bases, I, ref, rref, aligns, aligns_index );
     int nreads = bases.size( );

     // Find best alignment for each read.

     vec<int> best;
     GetBestAligns( bases, aligns, aligns_index, best );

     // Generate bias stats.

     if ( !BRIEF ) cout << "\n";
     vec<look_align> xaligns;
     for ( int id = 0; id < nreads; id++ )
     {    if ( best[id] < 0 || aligns_index[id].empty( ) ) continue;
          const look_align& la = aligns[ aligns_index[id][ best[id] ] ];
          if ( la.Errors( ) > MAX_ERRORS ) continue;
          xaligns.push_back(la);    }
     StartBias( bases, ref, rref, xaligns, BRIEF, MARKED_REFERENCE_FILE,
          SUMMARY_FORM );
     if ( BRIEF && !SUMMARY_FORM ) cout << "\n";    }
