///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FixScaffolds. Break scaffolds at uncovered spots (see FixScaffoldsCore.h
// for details).

#include "MainTools.h"
#include "Fastavector.h"
#include "Charvector.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Superb.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "math/NStatsTools.h"
#include "paths/Alignlet.h"
#include "paths/FixScaffoldsCore.h"
// MakeDepend: dependency Fastb

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_String_OrDefault( READS, "scaffold_reads" );
     CommandArgument_String_OrDefault( ALIGNS, "reads" );
     // input: <SCAFFOLDS_IN>.{contigs.fasta,superb}
     CommandArgument_String_OrDefault( SCAFFOLDS_IN, "linear_scaffolds" );
     // output: <SCAFFOLDS_OUT>.{contigs.fasta,assembly.fasta,superb}
     CommandArgument_String_OrDefault( SCAFFOLDS_OUT, "fixed_scaffolds" );
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Int_OrDefault(MIN_REACH_AWAY, 4);
     CommandArgument_String_OrDefault(TRACE_COVERAGE_OF, "");
     EndCommandArguments;

     // Define heuristic constants.

     const int min_dist_from_end = 1000;
     const double dev_mult = 4.0;
     const int blink_dev_mult = 1;
     const int min_scaffold = 1000;
     const int trim_back = 80;
     const int min_contig = 1000;
     const int min_spread = 10;

     // Dir and file names.

     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     String reads_file = run_dir + "/" + READS + ".fastb";
     String pairs_file = run_dir + "/" + READS + ".pairs";
     String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
     String index_file = sub_dir + "/" + ALIGNS + ".qltoutlet.index";
     
     String contigs_in = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
     String scaffolds_in = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     
     String contigs_out = sub_dir + "/" + SCAFFOLDS_OUT + ".contigs.fasta";
     String assembly_out = sub_dir + "/" + SCAFFOLDS_OUT + ".assembly.fasta";
     String scaffolds_out = sub_dir + "/" + SCAFFOLDS_OUT + ".superb";
     String log_file = sub_dir + "/" + SCAFFOLDS_OUT + ".FixScaffolds.log";
     
     // Parse input.
     
     vec<int> trace_ids;
     ParseIntSet( TRACE_COVERAGE_OF, trace_ids );

     // The log stream.

     ofstream log( log_file.c_str( ) );
     PrintCommandPretty( log );
     cout << "Sending verbose log to " << log_file << "\n" << endl;
     
     // Load.
     
     cout << Date( ) << ": loading fastavector" << endl;
     vec<fastavector> contigs;
     LoadFromFastaFile( contigs_in, contigs );
     
     cout << Date( ) << ": loading scaffolds" << endl;
     vec<superb> scaffolds;
     ReadSuperbs( scaffolds_in, scaffolds );

     cout << Date( ) << ": loading aligns" << endl;
     vec<alignlet> aligns0;
     BinaryReader::readFile( aligns_file, &aligns0 );
     vec<int> aligns0_index;
     BinaryReader::readFile( index_file, &aligns0_index );
     size_t nreads = MastervecFileObjectCount( reads_file );
     ForceAssertEq( nreads, aligns0_index.size( ) );
     
     cout << Date( ) << ": " << ToStringAddCommas( aligns0.size( ) ) 
          << " aligns found" << endl;
     
     cout << Date( ) << ": loading pairs" << endl;
     PairsManager pairs( pairs_file );
     pairs.makeCache( );

     // Run FixScaffoldsCore (contigs, scaffolds, and aligns will be changed).
    
     FixScaffoldsCore( trace_ids, pairs, MIN_REACH_AWAY, contigs,
		       scaffolds, aligns0, aligns0_index, aligns0_index, log );
     
     // Write remaining output files.

     if ( WRITE )
     {
         cout << Date( ) << ": writing output files" << endl;
         Ofstream( cg_out, contigs_out );
         for (int ii=0; ii<contigs.isize( ); ii++)
           contigs[ii].Print( cg_out, "contig_" + ToString( ii ) );
         cg_out.close();

         WriteSuperbs( scaffolds_out, scaffolds );

         WriteScaffoldedFasta( assembly_out, contigs, scaffolds );

         SystemSucceed( "Fastb NH=True PRE=/ IN=" + assembly_out );
     }

     cout << Date() << ": done" << endl;
     Scram();
     
}
