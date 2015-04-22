///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "SupersHandler.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/ReadLoc.h"
#include "paths/assisted/CInsertsDB.h"
#include "util/RunCommand.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * InsertsCloudStats
 *
 * Print basic stats on the local clouds of inserts.
 *
 * ASSEMBLY: head name of assembly
 * MIN_SEP: min sep for inserts
 * MAX_SEP: min sep for inserts
 * NUM_THREADS: use all if 0
 * VERBOSE: print a detailed list of all inserts (large output)
 * ARCHIVE: if true, send output to file (rather than cout)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String( ASSEMBLY );
  CommandArgument_Int_OrDefault( MIN_SEP, 50 );
  CommandArgument_Int_OrDefault( MAX_SEP, 15000 );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  CommandArgument_Bool_OrDefault( ARCHIVE, False );
  EndCommandArguments;
  
  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );
  
  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String assembly_head = sub_dir + "/" + ASSEMBLY;
  String locs_head = assembly_head + ".contigs";
  String locs_file = locs_head + ".readlocs";
  String supers_file = assembly_head + ".superb";
  String log_file = assembly_head + ".InsertsCloudStats.out";

  vec<String> needed;
  needed.push_back( supers_file );
  needed.push_back( locs_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  ofstream archive_out;
  if ( ARCHIVE ) {
    archive_out.open( log_file.c_str( ) );
    PrintCommandPretty( archive_out );
    cout << " Sending output to " << log_file << "\n" << endl;
  }
  ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;
  
  // Load.
  out << Date( ) << ": loading supers" << endl;
  shandler supers( -1, supers_file );
  
  read_locs_on_disk locs_parser( locs_head, run_dir );
  
  // Loop over all supers.
  vec< vec<String> > table;

  if ( ! VERBOSE ) {
    vec<String> legend
      = MkVec( String( "  s_id" ),
	       String( "s_len" ),
	       String( "#tot" ),
	       String( "#valid" ),
	       String( "#separated" ),
	       String( "#loners" ),
	       String( "cloud_mean_size" ) );
    table.push_back( legend );
  }

  out << Date( ) << ": loop over " << supers.Size( ) << " supers" << endl;
  if ( VERBOSE ) out << endl;
  for (int sid=0; sid<supers.Size( ); sid++) {
    if ( ! VERBOSE ) Dot( out, sid );
    CInsertsDB insdb( MIN_SEP, MAX_SEP, &sid, &supers, &locs_parser );

    vec< vec<String> > vtable;
    if ( VERBOSE ) vtable.reserve( insdb.Size( ) );
    
    vec<int> cloud_sizes;
    cloud_sizes.reserve( insdb.Size( ) );
    for (int jj=0; jj<(int)insdb.Size( ); jj++) {
      vec<int> icloud;
      insdb.BuildCloud( jj, icloud );
      if ( VERBOSE ) vtable.push_back( insdb.LineInfo( jj ) );
      else cloud_sizes.push_back( icloud.isize( ) );
    }
    String str_cloud_mean = "na";
    if ( cloud_sizes.size( ) > 0 ) {
      double cloud_mean = Mean( cloud_sizes );
      String strm = ToString( cloud_mean, 1 );
      String strd = ToString( StdDev( cloud_sizes, cloud_mean ), 1 );
      str_cloud_mean = strm + " +/- " + strd;
    }

    if ( VERBOSE ){
      const superb& sup = supers[sid];
      out << "SUPER " << sid
	   << "  " << ToStringAddCommas( sup.TrueLength( ) ) << " bases,"
	   << "  " << ToStringAddCommas( sup.Ntigs( ) ) << " contigs\n";
      PrintTabular( out, vtable, 3, "rrl" );
      out << endl;
    }
    else {
      vec<String> line
	= MkVec( ToString( sid ),
		 ToStringAddCommas( supers[sid].TrueLength( ) ),
		 ToStringAddCommas( insdb.NLocsTotal( ) ),
		 ToStringAddCommas( insdb.NLocsOfType( 0 ) ),
		 ToStringAddCommas( insdb.NLocsOfType( 1 ) ),
		 ToStringAddCommas( insdb.NLocsOfType( 2 ) ),
		 str_cloud_mean );
      table.push_back( line );
    }
    
  }
  if ( ! VERBOSE ) out << "\n" << endl;

  // Print table.
  if ( ! VERBOSE ) {
    PrintTabular( out, table, 3, "rrrrrrr" );
    out << endl;
  }
  
  // Done.
  out << Date( ) << ": done" << endl;
  
}

