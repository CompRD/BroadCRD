///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* InterleavedFastb2Pairs
 * 
 * A module to generate .pairs file from fastb containing interleaved paired sequences (we assume
 * that read pairs have the same average separation and standard deviation = belong to one library).
 *
 ******************************************************************************/

#include "MainTools.h"
#include "PairsManager.h"


int main(int argc, char **argv)
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String_Doc(IN, "input fastb file");
  CommandArgument_Int_Doc(SEP, "expected separation of reads");
  CommandArgument_Int_Doc(STDEV, "exptected standard deviation of separations");
  CommandArgument_String_Doc(LIB, "library name");

  EndCommandArguments;
  
  ForceAssert( IsRegularFile(IN) );
  


  String pairsFileOut = IN;
  pairsFileOut.ReplaceBy( ".fastb", ".pairs" );
  PRINT( pairsFileOut );
  
  if ( IsRegularFile( pairsFileOut ) )
    FatalErr("pairsFileOut=" + pairsFileOut + " already exists, please remove it if you are certain that you want it overwritten.");
    

  const size_t nreads = MastervecFileObjectCount( IN );
  PRINT(nreads);

  // ---- Create a PairsManager with all the read pairings.
  //      Assumes all paired reads have the same length within a library
  cout << Date() << ": Creating PairsManager." << endl;
  PairsManager pairs(nreads);
  
  pairs.addLibrary( SEP, STDEV, LIB );

  ForceAssert( nreads % 2 == 0 );

  size_t npairs = nreads / 2;

  for ( size_t pi = 0; pi < npairs; pi++ )
    pairs.addPairToLib( 2 * pi, 2 * pi + 1, 0 );


  //---- Save PairsManager.
  pairs.Write( pairsFileOut);


  // ---- Output some statistics
  cout << Date() << ": Library statistics:" << endl << endl;
  pairs.printLibraryStats(cout);
  cout << endl;

  command.PrintTheCommandPretty( cout );

  cout << Date() << ": Done!" << endl;
}
