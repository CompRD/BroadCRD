///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "PairsManager.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/LongReadsAligner.h"
#include "util/RunCommand.h"

// MakeDepend: dependency AlignPacBio

/**
 * BasicInputStats
 *
 * Generate a report on the input data of an AllPaths assembly. The
 * argument names match (some of) those of AllPathsLG.
 *
 * WARNING: if the file long_reads_orig.fastb is found in RUN, then
 * the code will align the long reads against the contigs of the final
 * assembly. This will take a while (log is reported in the file
 * AlignPacBio.out in SUBDIR).
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;  
  CommandArgument_String( PRE );
  CommandArgument_String( REFERENCE_NAME );
  CommandArgument_String( DATA_SUBDIR );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  EndCommandArguments;  
  
  // Dir and file names.
  String ref_dir = PRE + "/" + REFERENCE_NAME;
  String data_dir = ref_dir + "/" + DATA_SUBDIR;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String gsize_file = ref_dir + "/genome.size";
  String validate_file = data_dir + "/ValidateAllPathsInputs_core_stats.out";
  String libstats_file = run_dir + "/makeinfo/jump_reads_ec.distribs.out.SamplePairedReadDistributions";
  String longreads_file = run_dir + "/long_reads_orig.fastb";
  String longreads_out_file = sub_dir + "/long_reads_orig_on_final.contigs.out";
  String longreads_qlt_file = sub_dir + "/long_reads_orig_on_final.contigs.qlt";

  // These are needed.
  vec<String> needed;
  needed.push_back( gsize_file );
  needed.push_back( validate_file );
  needed.push_back( libstats_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Load genome size.
  size_t gsize = FirstLineOfFile( gsize_file ).Int( );
  if ( gsize < 1 ) {
    cout << "Fatal error: the file genome.size reports a genome size of "
	 << gsize << " bases. Leaving now.\n" << endl;
    return 1;
  }

  // These are filled upon request (one per type).
  vec< vec<PM_LibraryStats> > stats( 3 );

  // Basic reads input table.
  vec< vec<String> > table;
  
  vec<String> line = MkVec( ToString( " type" ),
			    ToString( "name" ),
			    ToString( "cov" ),
			    ToString( "length" ),
			    ToString( "mean" ),
			    ToString( "i_sep" ),
			    ToString( "i_dev" ) );
  table.push_back( line );

  {
    ifstream in( validate_file.c_str( ) );
    while( in ) {
      String oneline;
      getline( in, oneline );
      if ( ! in ) break;
      
      vec<String> tokens;
      Tokenize( oneline, ',', tokens );
      if ( tokens.size( ) < 1 ) continue;

      bool isF = tokens[0] == "frag";
      bool isJ = tokens[0] == "jump";
      bool isLJ = tokens[0] == "long_jump";
      if ( ! ( isF || isJ || isLJ ) ) continue;
      String str_type = isF ? "[f]" : ( isJ ? "[j]" : "[J]" );
      int id_type = isF ? 0 : ( isJ ? 1 : 2 );

      if ( stats[id_type].size( ) < 1 ) {
	String pairs_file = run_dir + "/";
	if ( id_type == 0 )      pairs_file += "frag_reads_corr_cpd.pairs";
	else if ( id_type == 1 ) pairs_file += "jump_reads_ec.pairs";
	else                     pairs_file += "long_jump_reads_ec.pairs";

	long n_reads_tot = 0;
	ReadPairsManagerLibInfo( pairs_file, n_reads_tot, stats[id_type] );
      }
      
      double cov = SafeQuotient( tokens[2].Int( ), long( gsize ) );
      
      String range;
      if ( tokens[3] == tokens[4] ) range = tokens[3];
      else range = "[" + tokens[3] + ", " + tokens[4] + "]";

      int sep = -666;
      int dev = -666;
      if ( tokens[1] != "Unpaired" ) {
	for (size_t ii=0; ii<stats[id_type].size( ); ii++) {
	  const PM_LibraryStats &stat = stats[id_type][ii];
	  if ( stat.name != tokens[1] ) continue;
	  sep = stat.sep;
	  dev = stat.sd;
	  break;
	}
      }
      String str_sep = sep == -666 ? "na" : ToString( sep );
      String str_dev = dev == -666 ? "na" : ToString( dev );

      line.clear( );
      line = MkVec( str_type,
		    tokens[1],
		    ToString( cov, 1 ),
		    range,
		    tokens[5],
		    str_sep,
		    str_dev );
      table.push_back( line );
    }
    in.close( );
  }

  cout << "BASIC READ STATISTICS\n"
       << "\n"
       << " type: read type ([f]: fragment, [j]: jump, [J]: long jump)\n"
       << " name: library name\n"
       << " cov: coverage (based on a genome size of " << gsize << " bases)\n"
       << " length: length (or range of lengths) of reads\n"
       << " mean: mean read length\n"
       << " i_sep: estimated mean of insert separations\n"
       << " i_dev: estimated deviation of insert separations\n"
       << "\n";
  
  PrintTabular( cout, table, 2, "rrrrrrr" );
  cout << endl;

  // Physical coverage table.
  {
    ifstream in( libstats_file.c_str( ) );

    String title = "[Libraries physical coverage (cummulative) table]";
    String in_line = "PERFSTAT: BLOCK_START " + title;
    String out_line = "PERFSTAT: BLOCK_STOP";
    
    cout << "LIBRARIES PHYSICAL COVERAGE (CUMULATIVE) TABLE\n\n";

    String line;
    while ( in ) {
      getline( in, line );
      if ( !in ) break;
      if ( line != in_line ) continue;
      while ( in ) {
	getline( in, line );
	if ( ! in || line == out_line ) break;
	cout << line << "\n";
      }
    }
    in.close( );
  }
  cout << endl;

  // PacBio coverage table.
  if ( IsRegularFile( longreads_file ) ) {
    const int SAMPLE_SIZE = 10000;
    const String contigs_file = sub_dir + "/final.contigs.fastb";
    const String aligns_file = sub_dir + "/longreads_on_finalcontigs.qltout";

    String log_file = sub_dir + "/AlignPacBio.out";
    String theCommand
      = "AlignPacBio READS=" + longreads_file
      + " REFERENCE=" + contigs_file
      + " QLTOUT=" + aligns_file
      + " SAMPLE_SIZE=" + ToString( SAMPLE_SIZE );
    RunCommandWithLog( theCommand, log_file );

    vec<look_align> aligns;
    LoadLookAligns( aligns_file, aligns );
    vecbvec queries( longreads_file );
    vecbvec targets( contigs_file );
    
    LongReadsCoverages( targets, queries, aligns, SAMPLE_SIZE, cout );
  }

}
