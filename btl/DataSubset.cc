///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

#include "Basevector.h"
#include "ParseSet.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "system/System.h"
#include "util/RunCommand.h"

/**
 * DataSubset
 *
 * Build a subset from a given data set (fastb/qualb/pairs), using
 * only the pairs containing the specified read ids.
 */
int main( int argc, char* argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( HEAD_IN );
  CommandArgument_String( HEAD_OUT );
  CommandArgument_String( READ_IDS );   // parsed with ParseLongLongSet
  EndCommandArguments;

  // Dir and file names.
  String in_bases_file = HEAD_IN + ".fastb";
  String in_quals_file = HEAD_IN + ".qualb";
  String in_pairs_file = HEAD_IN + ".pairs";

  String out_dir = Dirname( HEAD_OUT );
  String out_bases_file = HEAD_OUT + ".fastb";
  String out_quals_file = HEAD_OUT + ".qualb";
  String out_pairs_file = HEAD_OUT + ".pairs";
  String out_select_pids_file = HEAD_OUT + ".select_pids";
  String log_file = HEAD_OUT + ".DataSubset.log";

  vec<String> needed;
  needed.push_back( in_bases_file );
  needed.push_back( in_quals_file );
  needed.push_back( in_pairs_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( out_dir );
  
  // Log stream.
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << "Sending log to " << log_file << "\n" << endl;
  
  // Load reads and pairs, and select pair ids / read ids.
  log << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( in_pairs_file );

  vec<longlong> given_rids;
  ParseLongLongSet( READ_IDS, given_rids );
  
  vec<longlong> sel_pids;
  sel_pids.reserve( given_rids.size( ) );
  for (size_t ii=0; ii<given_rids.size( ); ii++)
    sel_pids.push_back( pairs.getPairID( given_rids[ii] ) );
  sort( sel_pids.begin( ), sel_pids.end( ) );
  sel_pids.erase( unique( sel_pids.begin(), sel_pids.end() ), sel_pids.end() );

  vec<int> sel_rids;   // SparseRead( ) wants ints
  sel_rids.reserve( 2 * sel_pids.size( ) );
  for (size_t ii=0; ii<sel_pids.size( ); ii++) {
    sel_rids.push_back( (int)pairs.ID1( sel_pids[ii] ) );
    sel_rids.push_back( (int)pairs.ID2( sel_pids[ii] ) );
  }
  sort( sel_rids.begin( ), sel_rids.end( ) );
  const size_t nsel_rids = sel_rids.size( );

  // Save ids of selected pairs.
  {
    ofstream sel_out( out_select_pids_file.c_str( ) );
    for (size_t ii=0; ii<sel_pids.size( ); ii++)
      sel_out << sel_pids[ii] << "\n";
    sel_out.close( );
  }
  
  // Load bases and quals.
  log << Date( ) << ": loading " << nsel_rids << " reads" << endl;
  vecbvec bases;
  vecqvec quals;
  bases.SparseRead( in_bases_file, sel_rids, 0 );
  quals.SparseRead( in_quals_file, sel_rids, 0 );
  
  // Build output data.
  log << Date( ) << ": building subset dataset" << endl;
  vecbvec out_bases;
  vecqvec out_quals;
  PairsManager out_pairs( nsel_rids );
  out_bases.reserve( nsel_rids );
  out_quals.reserve( nsel_rids );
		     
  vec<PM_LibraryStats> libs = pairs.getLibraryStats( );
  for (int lib_id=0; lib_id<libs.isize( ); lib_id++) {
    int sep = libs[lib_id].sep;
    int dev = libs[lib_id].sd;
    String name = libs[lib_id].name;
    out_pairs.addLibrary( sep, dev, name );
  }
  
  longlong current_id = 0;
  for (size_t sid=0; sid<sel_pids.size( ); sid++) {
    longlong pid = sel_pids[sid];
    longlong id1 = pairs.ID1( pid );
    longlong id2 = pairs.ID2( pid );
    size_t lib_id = pairs.libraryID( pid );

    out_bases.push_back( bases[id1] );
    out_bases.push_back( bases[id2] );
    out_quals.push_back( quals[id1] );
    out_quals.push_back( quals[id2] );
    out_pairs.addPairToLib( current_id, current_id+1, lib_id );
  
    current_id += 2;
  }

  // Save.
  log << Date( ) << ": saving subset dataset" << endl;
  out_bases.WriteAll( out_bases_file );
  out_quals.WriteAll( out_quals_file );
  out_pairs.Write( out_pairs_file );
  
  // Done.
  String date = Date( );
  cout << date << ": DataSubset done" << endl;
  log << date << ": DataSubset done" << endl;
  log.close( );

}
