///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


/* GORDIAN KNOT: From the name of a legendary knot tied to a pole near the
 * temple of Zeus in Gordium.  It was prophesied that whoever loosed the knot
 * would become ruler of all Asia.  Alexander the Great solved the puzzle by
 * slicing through the knot and took it as a sign of Zeus's favor.  He then
 * proceeded to conquer much of the known world.
 *         -- Wiktionary
 *
 *
 *
 *
 * In CutTheGordianKnot, we attempt to disentangle the knottiness in an assembly
 * produced via the RunAllPathsLG pipeline.  If we can do this, we can do
 * anything!
 *
 * This module is a wrapper for the GordianKnot class.
 *
 *
 * Josh Burton
 * November 2009
 *
 ******************************************************************************/



#include "MainTools.h"
#include "PairsManager.h"
#include "paths/GordianKnot.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "feudal/BinaryStream.h"




int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  
  // File names and directories
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String_OrDefault( HYPER_IN, "hyper" );
  CommandArgument_String_OrDefault( HYPER_OUT, "hyper.unknotted" );
  CommandArgument_String_OrDefault( FRAG_READS, "frag_reads_corr" );
  
  // Heuristic parameters for the GordianKnot
  CommandArgument_Bool_OrDefault_Doc( CHEAT, False, "Call RemoveNongenomicAdjacencies instead of RemoveUnsupportedAdjacencies" );
  CommandArgument_Int_OrDefault( N_READS_TO_SUPPORT, 1 );
  CommandArgument_Int_OrDefault_Doc( N_PAIRS_TO_SUPPORT, 1, "In GordianKnot::RemoveUnsupportedAdjacencies, require N_READS_TO_SUPPORT reads *or* N_PAIRS_TO_SUPPORT read pairs to verify an adjacency.  If these are set to -1, then no support will ever be given." );
  CommandArgument_Double_OrDefault_Doc( MIN_HAIR_COV, 0.5, "In GordianKnot::Shave, if a hair edge has coverage more than MIN_HAIR_COV times the average HBV coverage, keep it" );
  CommandArgument_Double_OrDefault_Doc( MIN_BUBBLE_COV, 0.2, "In GordianKnot::PopBubbles, if an edge has coverage more than MIN_BUBBLE_COV times the average HBV coverage, keep it" );
  //CommandArgument_Double_OrDefault_Doc( LOOP_COV_TOLERANCE, 0.25, "In GordianKnot::UnravelLoops, only unravel a loop if its coverage is within LOOP_COV_TOLERANCE of a multiple of the coverage of nearby edges" );
  
  // Flags to control different kinds of output
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  CommandArgument_Bool_OrDefault( WRITE, True );
  
  EndCommandArguments;
  
  cout << Date() << ": Beginning CutTheGordianKnot" << endl;
  
  // Directory names.
  const String data_dir = PRE + "/" + DATA;
  const String  run_dir = PRE + "/" + DATA + "/" + RUN;
  const String  sub_dir = PRE + "/" + DATA + "/" + RUN + "/ASSEMBLIES/" + SUBDIR;
  
  // Input files.
  cout << Date() << ": Loading input files" << endl;
  
  vecbasevector frag_reads( run_dir + "/" + FRAG_READS + ".fastb" );
  PairsManager frag_pairs( run_dir + "/" + FRAG_READS + ".pairs" );
  
  HyperKmerPath hkp( sub_dir + "/" + HYPER_IN );
  KmerBaseBroker kbb( sub_dir, hkp.K() );
  
  // Main algorithm: create a GordianKnot object and call GordianKnot functions.
  cout << Date() << ": Creating the GordianKnot and attempting to cut it" << endl;
  Mkdir777( sub_dir + "/temp" );
  GordianKnot knot( hkp, &frag_reads, &frag_pairs, &kbb, sub_dir + "/temp", VERBOSE );
  knot.Shave( MIN_HAIR_COV );
  if ( CHEAT )
    knot.RemoveNongenomicAdjacencies( data_dir + "/genome" );
  else
    knot.RemoveUnsupportedAdjacencies( N_READS_TO_SUPPORT, N_PAIRS_TO_SUPPORT );
  knot.PopBubbles( MIN_BUBBLE_COV );
  //knot.UnravelLoops( LOOP_COV_TOLERANCE ); // doesn't do anything useful
  hkp = knot.GetHKP( );
  
  
  // Write output files.
  if (WRITE) {
    const String hyper_out = sub_dir + "/" + HYPER_OUT;
    BinaryWriter::writeFile( hyper_out, hkp );
    Ofstream( dot, hyper_out + ".dot" );
    hkp.PrintSummaryDOT0w(dot);
    hkp.DumpFasta( hyper_out + ".fasta", kbb );
    hkp.DumpFastb( hyper_out + ".fastb", kbb );
  }
  
  
  cout << Date() << ": Done with CutTheGordianKnot!" << endl;
  return 0;
}
