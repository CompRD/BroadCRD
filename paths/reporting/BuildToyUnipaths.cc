///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "VecUtilities.h"
#include "kmers/SupportedKmerShapes.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "util/RunCommand.h"
// MakeDepend: dependency BuildUnipathAdjGraph
// MakeDepend: dependency CommonPather
// MakeDepend: dependency Fastb
// MakeDepend: dependency MakeRcDb
// MakeDepend: dependency Unipather

/**
 * BuildToyUnipaths
 *
 * Build tiny uinpaths and accessory files. Useful for presentations,
 * or to teach what unipaths are.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;

  // Specify input genome: either a random chunk from a given fastb...
  CommandArgument_String_OrDefault( FASTB, "" );
  CommandArgument_Int_OrDefault( ID, "" );
  CommandArgument_Int_OrDefault( START, "" );
  CommandArgument_Int_OrDefault( LEN, "" );

  // ... or a full fasta file.
  CommandArgument_String_OrDefault( FASTA, "" );

  // Kmer size.
  CommandArgument_Int( K );

  // Base output dir.
  CommandArgument_String( OUT_BASE );

  // Do not use cached data, re-generate all from scratch.
  CommandArgument_Bool_OrDefault( FORCE, False );

  EndCommandArguments;

  // Dir and file names.
  String strK = ToString( K );
  String str_id = ToString( ID );
  String str_start = ToString( START );
  String str_len = ToString( LEN );

  String fasta_head = str_id + "." + str_start + "." + str_len;

  String out_dir = OUT_BASE + ( FASTA == "" ? "/fasta_head" : "" );

  String target_head = out_dir + "/target";
  String target_source_file = out_dir + "/target.source";
  String target_fastb_file = out_dir + "/target.fastb";
  String target_paths_head = out_dir + "/target.paths";
  String target_paths_file = out_dir + "/target.paths.k" + strK;
  String target_unipaths_file = out_dir + "/target.unipaths.k" + strK;
  String target_pathsdb_file = out_dir + "/target.pathsdb.k" + strK;
  String target_adjgraph_file = out_dir + "/target.unipath_adjgraph.k" + strK;
  
  String all_kmers_file = out_dir + "/all_kmers.k" + strK + ".txt";
  String unipaths_file = out_dir + "/unipaths.k" + strK + ".txt";

  Mkpath( out_dir );

  // Generate vecbvec of target sequence.
  if ( FORCE || ! IsRegularFile( target_fastb_file ) ) {
    if ( FASTA != "" ) {
      String theCommand
	= String( "Fastb" )
	+ String( " FILE=" + FASTA )
	+ String( " PRE= " );
      RunCommand( theCommand );
    }
    else {
      vecbvec all_target( FASTB );
      vecbvec select;
      select.push_back( bvec( all_target[ID], START, LEN ) );
      select.WriteAll( target_fastb_file );
      
      ofstream out( target_source_file.c_str( ) );
      out << "source: " << FASTB << "\n"
	  << "ID:     " << ID << "\n"
	  << "START:  " << START << "\n"
	  << "LEN:    " << LEN << "\n"
	  << endl;
      out.close( );
    }
  }
  
  // CommonPather (target).
  if ( FORCE || ! IsRegularFile( target_paths_file ) ) {
    String theCommand
      = String( "CommonPather" )
      + String( " K=" + strK )
      + String( " READS_IN=" + target_fastb_file )
      + String( " PATHS_OUT=" + target_paths_head );
    RunCommand( theCommand );
  }

  // MakeRcDb (target).
  if ( FORCE || ! IsRegularFile( target_pathsdb_file ) ) {
    String theCommand
      = String( "MakeRcDb" )
      + String( " K=" + strK )
      + String( " READS=" + target_head )
      + String( " PRE= DATA= RUN=" );
    RunCommand( theCommand );
  }
  
  // Unipather (target).
  if ( FORCE || ! IsRegularFile( target_unipaths_file ) ) { 
    String theCommand
      = String( "Unipather" )
      + String( " K=" + strK )
      + String( " READS=" + target_head )
      + String( " PRE= DATA= RUN=" );
    RunCommand( theCommand );
 }

  // BuildUnipathAdjGraph (target).
  if ( FORCE || ! IsRegularFile( target_adjgraph_file ) ) { 
    String theCommand
      = String( "BuildUnipathAdjGraph" )
      + String( " K=" + strK )
      + String( " READS=" + target_head )
      + String( " PRE= DATA= RUN=" );
    RunCommand( theCommand );
  }
  
  // Load.
  vecbvec target_bases( target_fastb_file );
  vecKmerPath paths( target_paths_file );
  vecKmerPath unipaths( target_unipaths_file );
  KmerBaseBroker kbb( String( "" ), K, target_head );

  // All target's kmers (sorted).
  {
    vec<String> kmers;
    for (int ii=0; ii<(int)paths.size( ); ii++) {
      const KmerPath &kpath = paths[ii];
      const int klen = kpath.KmerCount( );
      longlong pos = 0;
      for (int jj=0; jj<kpath.NSegments( ); jj++) {
	const KmerPathInterval &seg = kpath.Segment( jj );
	for (longlong kk=seg.Start( ); kk<=seg.Stop( ); kk++) {
	  const bvec &as_bases = kbb.Bases( kk );
	  kmers.push_back( as_bases.ToString( ) );
	}
      }
    }
    sort( kmers.begin( ), kmers.end( ) );
    kmers.erase( unique( kmers.begin( ), kmers.end( ) ), kmers.end( ) );
    
    ofstream out( all_kmers_file.c_str( ) );
    out << "Found " << kmers.size( ) << " distinct kmers\n\n";
    for (int ii=0; ii<kmers.isize( ); ii++)
      out << kmers[ii] << "\n";
    out << endl;
    out.close( );
  }
  
  // Print unipaths of target.
  {
    vec< vec<String> > table;
    for (int ii=0; ii<(int)unipaths.size( ); ii++) {
      bvec bases( kbb.Seq( unipaths[ii] ) );
      vec<String> line = MkVec( ToString( ii ), bases.ToString( ) );
      table.push_back( line );
    }

    ofstream out (unipaths_file.c_str( ) );
    out << "Original sequence(s):\n\n";
    for (int ii=0; ii<(int)target_bases.size( ); ii++)
      out << target_bases[ii].ToString( ) << "\n";
    out << "\nFound " << unipaths.size( ) << " unipaths\n\n";
    PrintTabular( out, table, 3, "rl" );
    out << endl;
    out.close( );
  }
  
}

