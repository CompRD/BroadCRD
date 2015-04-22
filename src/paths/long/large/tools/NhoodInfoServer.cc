///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// NhoodInfo.
//
// Notes on server-mode operation:
// - Invoked by setting SERVER_DIR.
// - Runs in parallel on SERVER_DIR/n, where n = 1,...,N, and N = max threads.
// - Expects request generator to choose n at random.
// - Each thread first deletes old files.
// - Then it looks for something.req, a request.  This file is then deleted.
// - Generates something.txt and if successful something.png.
// - Runs until killed.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"

#include "paths/long/large/tools/NhoodInfoCore.h"


int main( int argc, char *argv[] )
{    
     BeginCommandArgumentsNoHeader;
     CommandDoc( "Generate a DOT file corresponding to a region in an assembly." );
     CommandArgument_String_OrDefault_Doc(DIR_IN, ".", 
          "location of the DISCOVAR final assembly directory containing a.s.hbx or a.hbx files");
     CommandArgument_Bool_Abbr_OrDefault_Doc(LABEL_CONTIGS, LC, False, 
          "label contigs");
     CommandArgument_String_OrDefault_Doc(SERVER_DIR, "", 
          "if nonempty, run in server mode on this directory, which must exist");
     CommandArgument_Bool_OrDefault_Doc(SEQ_LOOKUP, True, 
          "build assembly lookup table to allow fast lookup by sequence - "
          "slow initialization. Default is True for server.");
     EndCommandArguments;

     // Check dot executable, etc.

     TestDot( );

     // Load.

     NhoodInfoEngine engine;

     engine.Initialize(DIR_IN, SEQ_LOOKUP, True, False);

     // Execute server mode.
     cout << "Entering server mode." << endl;

     engine.RunAsServer(SERVER_DIR, LABEL_CONTIGS);

     // Done.

     cout << Date( ) << ": done" << endl;    
}
