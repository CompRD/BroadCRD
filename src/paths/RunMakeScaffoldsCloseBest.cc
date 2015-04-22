///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "math/Functions.h"
#include "pairwise_aligners/AlignConsecutiveContigs.h"
#include "paths/Alignlet.h"
#include "paths/FixScaffoldsCore.h"
#include "paths/MakeScaffoldsCloseBest.h"
#include "paths/MakeScaffoldsScoredBest.h"
#include "paths/RegapSupers.h"
#include "paths/SaveScaffoldGraph.h"
#include "paths/ScaffoldsUtils.h"

/**
 * RunMakeScaffoldsCloseBest
 *
 * A single run of MakeScaffoldsCloseBest.  The code assumes that no
 * previous filtering (as in RemoveHighCNAligns) has run, and it does
 * its own internal filtering.
 * 
 * The following example (from a real assembly, that of Mouse 80Mb)
 * shows why we want to enforce a MAX_LINKS allowed value:
 *
 *               c1    R    R    R    R    c2 (~R)
 *         ... ----> [===][===][===][===] ---> ...
 * last contig in s1                      first contig in s2
 *
 * Contig c2 contains almost all of the repeat R, and there are
 * hundreds of links between the contigs c1 and c2 (all of which were
 * filtered out by RenoveHighCNAligns).  If we used them now we would
 * join supers s1 and s2 with a gap of about 2 Kb, while the real gap
 * is about 32 Kb.
 *
 * Arguments are documented in paths/MakeScaffoldsCloseBest.h.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads2" );
  CommandArgument_String_OrDefault( SCAFFOLDS_IN, "linear_scaffolds0.clean.patched" );
  CommandArgument_String_OrDefault( SCAFFOLDS_OUT, "linear_scaffolds0.clean.patched.ext" );
  CommandArgument_Int_OrDefault( MIN_LINKS, 2 );
  CommandArgument_Int_OrDefault( MAX_LINKS, 14 );
  CommandArgument_Int_OrDefault( MAX_OVERLAP, 15000 );
  CommandArgument_Int_OrDefault( MIN_SCAFFOLD_LEN, 100000 );
  CommandArgument_String_OrDefault( SCAFFOLD_GRAPH_OUT, "" );
  CommandArgument_String_OrDefault( SCAFFOLD_GRAPH_DOT, "" );
  CommandArgument_Int_OrDefault( VERBOSITY, 0 );
  CommandArgument_Int_OrDefault( MIN_LINKS_TO_PRINT, std::numeric_limits<int>::max( ) );
  EndCommandArguments;
  
  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String pairs_file = run_dir + "/" + READS + ".pairs";

  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file =  sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  String contigs_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
  String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
  String final_head = sub_dir + "/" + SCAFFOLDS_OUT;
  
  // Load.
  cout << Date( ) << ": loading contigs fasta" << endl;
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  size_t n_pairs = pairs.nPairs( );
  
  cout << Date( ) << ": done loading\n" << endl;
  
  // Run MakeScaffoldsCloseBest.
  ReportScaffoldsBrief( supers, MIN_LINKS, 0, cout );
  MakeScaffoldsCloseBest( supers, contigs, aligns, index, pairs, "", "",
			  ToString( MIN_LINKS ),
			  ToString( MAX_LINKS ),
			  MAX_OVERLAP,
			  MIN_SCAFFOLD_LEN,
			  False,  // SUCK_SCAFFOLDS
			  False,  // SHAVE_GRAPH
			  SCAFFOLD_GRAPH_OUT,
			  SCAFFOLD_GRAPH_DOT,
			  VERBOSITY,
			  MIN_LINKS_TO_PRINT );
  ReportScaffoldsBrief( supers, MIN_LINKS, 1, cout );
  cout << endl;

  // Save.
  cout << Date( ) << ": saving" << endl;
  SaveScaffoldAssembly( final_head, supers, contigs );
  WriteSummary( final_head + ".summary", supers);
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

