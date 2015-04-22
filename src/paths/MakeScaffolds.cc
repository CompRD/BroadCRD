///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeScaffolds.  Build scaffolds from assembly graph (run this after
// AlignPairsToHyperLG). See paths/MakeScaffodsCore.h for details.

// WARNING: This code is now replaced by MakeScaffoldsLG (the core
// algorithm in MakeScaffolds now resides in MakeScaffoldsCloseBest.h)

#include "MainTools.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "feudal/FeudalTools.h"
#include "paths/MakeScaffoldsCloseBest.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  
  // Input (from <PRE>/<DATA>/<RUN>):
  CommandArgument_String_OrDefault( READS, "all_reads" );
  
  // Input (from <PRE>/<DATA>/<RUN>/<SUBDIR>):
  CommandArgument_String_OrDefault( ALIGNS, "" );   // defaulted to READS
  CommandArgument_String_OrDefault( HYPER, "hyper_plus_edited" );
  CommandArgument_String_OrDefault( SCAFFOLDS_IN, "initial_scaffolds" );
  
  // Output saved in <PRE>/<DATA>/<RUN>/<SUBDIR>/<ASSEMBLIES>:
  //   <SCAFFOLDS_OUT>.{contigs.fasta,assembly.fasta,superb}
  CommandArgument_String_OrDefault( SCAFFOLDS_OUT, "scaffolds" );
  
  // If False do not save output assembly.
  CommandArgument_Bool_OrDefault( WRITE, True );

  // The following args are documented in paths/MakeScaffoldsCloseBest.h
  CommandArgument_String_OrDefault(MIN_PAIR_SEPS, "");
  CommandArgument_String_OrDefault(MAX_PAIR_SEPS, ""); 
  CommandArgument_String_OrDefault( MIN_LINKS, "{20,6,4,2}" );
  CommandArgument_Bool_OrDefault( SUCK_SCAFFOLDS, False );
  CommandArgument_Bool_OrDefault( SHAVE_GRAPH, False );
  CommandArgument_Bool_OrDefault( SCAFFOLD_GRAPH_ALT, False );
  CommandArgument_String_OrDefault( SCAFFOLD_GRAPH_OUT, "" );
  CommandArgument_String_OrDefault( SCAFFOLD_GRAPH_DOT, "" );
  CommandArgument_Int_OrDefault( VERBOSITY, 0 );
  CommandArgument_Int_OrDefault( MIN_LINKS_TO_PRINT, 1000000000 );
  CommandArgument_Double_OrDefault( CONNECTIVITY_LIMIT, 1000000 );
  CommandArgument_Int_OrDefault( MAX_OVERLAP2, std::numeric_limits<int>::max());
  EndCommandArguments;
  
  // Default for ALIGNS.
  if ( ALIGNS == "" ) ALIGNS = READS;
  
  // Define heuristic constants. These could be made into arguments.
  const int dev_mult = 3;
  const int pair_seps_bin_size = 100000; // large value (effectively no binning)
  
  // Dir and file names.

  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String reads_file = run_dir + "/" + READS + ".fastb";
  String pairs_file = run_dir + "/" + READS + ".pairs";
  String contigs_in_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
  String scaffolds_in_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file = sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  
  String out_suffix = SCAFFOLDS_OUT;
  String contigs_out_file = sub_dir + "/" + out_suffix + ".contigs.fasta";
  String assembly_out_file = sub_dir + "/" + out_suffix + ".assembly.fasta";
  String scaffolds_out_file = sub_dir + "/" + out_suffix + ".superb";
  
  // Load scaffold files (vec<fastavector> and superb.)
  cout << Date( ) << ": loading scaffolds" << endl;
  vec<superb> scaffolds;
  ReadSuperbs( scaffolds_in_file, scaffolds );

  vec<fastavector> scaffold_contigs;
  LoadFromFastaFile( contigs_in_file, scaffold_contigs );
  size_t ntigs = scaffold_contigs.size();
  
  // Check that the files are consistent.
  size_t n_tigs_in_scaffolds = 0;
  for ( size_t i = 0; i < scaffolds.size(); i++ ) {
    n_tigs_in_scaffolds += scaffolds[i].Ntigs();
    for ( int j = 0; j < scaffolds[i].Ntigs(); j++ ) {
      int tig = scaffolds[i].Tig(j);
      ForceAssertEq( (int)scaffold_contigs[tig].size(), scaffolds[i].Len(j) );
    }
  }
  ForceAssertGe( ntigs, n_tigs_in_scaffolds );
  
  // Load alignments and convert to lightweight form.
  // aligns0_index works as this: say aligns0_index[id] = ii, then:
  //  ii = -3 :   look_aligns not informative for this read 
  //              (e.g. both reads align to the same scaffold)
  //  ii = -2 :   no look_aligns for this read
  //  ii = -1 :   more than one look_align for this read (id)
  //  ii >= 0 :   aligns0[ii] is the only look_align for this read
  
  cout << Date( ) << ": loading aligns" << flush;
  vec<alignlet> aligns0;
  BinaryReader::readFile( aligns_file, &aligns0 );
  vec<int> aligns0_index;
  BinaryReader::readFile( index_file, &aligns0_index );
  longlong nreads = MastervecFileObjectCount(reads_file);
  cout << " (" << aligns0.size( ) << " aligns found)" << endl;

  // Load read pairs.
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( run_dir + "/" + READS + ".pairs" );
  
  // Run MakeScaffodsCloseBest
  
  MakeScaffoldsCloseBest( scaffolds,
			  scaffold_contigs,
			  aligns0,
			  aligns0_index,
			  pairs,
			  MIN_PAIR_SEPS,
			  MAX_PAIR_SEPS,
			  MIN_LINKS,
			  "",   // MAX_LINKS
			  MAX_OVERLAP2,
			  0,   // MIN_SCAFFOLD_LEN
			  SUCK_SCAFFOLDS,
			  SHAVE_GRAPH,
			  SCAFFOLD_GRAPH_OUT,
			  SCAFFOLD_GRAPH_DOT,
			  VERBOSITY,
			  MIN_LINKS_TO_PRINT );
  
  // Save results.
  
  if ( !WRITE ) exit( 0 );
  
  cout << "\n" << Date( ) << ": writing output files" << endl;
     
  WriteSuperbs( scaffolds_out_file , scaffolds );
  Ofstream( contig_fasta, contigs_out_file );
  for ( size_t i = 0; i < ntigs; i++ )
    scaffold_contigs[i].Print( contig_fasta, "contig_" + ToString(i) );
  contig_fasta.close();
  
  WriteScaffoldedFasta( assembly_out_file, scaffold_contigs, scaffolds );
  
  cout << Date( ) << ": done" << endl;
  
}
