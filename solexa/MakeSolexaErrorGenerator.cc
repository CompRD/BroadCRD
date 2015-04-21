/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/** MakeSolexaErrorGenerator creates error generator tables based on a given
    set of real read alignments.
    The tables produced by this code can be imported into an ErrorGenerator
    object and then used to produce simulated reads with an error profile
    matching the orginal data.
*/

#include "MainTools.h"
#include "math/Functions.h"
#include "Basevector.h"
#include "FetchReads.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "solexa/FourBase.h"
#include "solexa/SolexaTools.h"


int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HEAD);

     CommandArgument_Double_OrDefault(MIN_RATIO, 2.0);
     CommandArgument_Bool_OrDefault(FILTER_PILEUPS, True);
     // Use only best alignment for each Solexa read, or use all?
     CommandArgument_Bool_OrDefault(BEST_ONLY, True);
     // Exclude error profiles from low quality reads?
     CommandArgument_Bool_OrDefault(EXCLUDE_LOW_QUALITY, False);
     // If NUM_BASES=0 then don't trim profiles
     CommandArgument_Int_OrDefault(NUM_BASES, 30);
     // File to store error generator
     CommandArgument_String(TABLE_OUT);
     // Show Error Mode Statistics
     CommandArgument_Bool_OrDefault(SHOW_ERROR_MODES, False);
     // Show By Base Error Profile
     CommandArgument_Bool_OrDefault(SHOW_ERROR_PROFILE, False);
     // Optional By Base Error Profile Filename
     CommandArgument_String_OrDefault(ERROR_PROFILE_OUT, "");

     EndCommandArguments;

    if (!BEST_ONLY && (EXCLUDE_LOW_QUALITY || FILTER_PILEUPS) )
      InputErr("If BEST_ONLY=False, then EXCLUDE_LOW_QUALITY and FILTER_PILEUPS"
	       " must be False too");
    
     cout << "Importing data\n";

     // Set up data structures.
     vecbasevector bases, basesrc, ref, rref;
     vecqualvector quals;
     VecFourBaseVec I;
     vec<look_align> aligns;
     vec< vec<int> > aligns_index;

     // Load the solexa data
     LoadSolexaData( HEAD, bases, basesrc, I, ref, rref, aligns, aligns_index);
     int nreads = bases.size( );

     // For each read, find the minimum quality of its first 10 bases.
     vec<float> minq10;
     GetMinQ10( I, nreads, minq10 );

     // Find best alignment for each read.
     cout << "Aligning Reads\n";
     vec<int> best;
     GetBestAligns( bases, aligns, aligns_index, best );

     // Filter out read pileups.
     if (FILTER_PILEUPS)
     {    int nfiltered = FilterPileups( bases, ref, aligns, aligns_index,
               minq10, best, MIN_RATIO );
          cout << nfiltered << " reads filtered out by pileups\n";    }

     // Exclude Low Quality Reads
     if (EXCLUDE_LOW_QUALITY)
       for( int id = 0; id < nreads; id++ )
	 if ( minq10[id] < MIN_RATIO)
	   best[id] = -1;

     // Write out list of best aligns
     BinaryWriter::writeFile( HEAD + ".bestaligns.temp", best );

     // Build Error Generator Tables
// MakeDepend: dependency MakeErrorGenerator
     SystemSucceed( "MakeErrorGenerator READS_IN=" + HEAD + ".fastb"
 		    + " QUALS_IN=" + HEAD + ".new.qualb"
		    + " ALIGNS_IN=" + HEAD + ".qltout"
		    + " REF_IN=" + HEAD + ".reference.fasta"
		    + " BEST_ALIGNS_IN=" + HEAD + ".bestaligns.temp"
		    + " BEST_ONLY=" + (BEST_ONLY ? "True" : "False" )
		    + " NUM_BASES=" + ToString(NUM_BASES)
		    + " SHOW_ERROR_MODES=" + (SHOW_ERROR_MODES ? "True" : "False" )
		    + " SHOW_ERROR_PROFILE=" + (SHOW_ERROR_PROFILE ? "True" : "False" )
		    + " ERROR_PROFILE_OUT=" + ERROR_PROFILE_OUT
		    + " TABLE_OUT=" + TABLE_OUT);

}

