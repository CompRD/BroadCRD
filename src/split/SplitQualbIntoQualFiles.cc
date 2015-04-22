/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Qualvector.h"
#include "ParseSet.h"
#include "STLExtensions.h"
#include "String.h"
#include "Vec.h"
#include "math/Functions.h"
#include "system/ParsedArgs.h"
#include "system/System.h"

/**
 * SplitQualbIntoQualFiles
 *
 * Split a qualb file into fixed-size qual files. See also its twin
 * module, SplitFastbIntoFastaFiles.
 *
 * RUN_DIR: where the fastb file is
 * OUT_DIR: as a subdir of RUN_DIR
 * HEAD: head of qualb and ids file
 * SELECT: use these ids only (full path name, parsed with ParseIntSet)
 * BATCH_SIZE: how many reads in each qual file
 * NAMES: load ids file and use read names in the qual output
 * OVERWRITE: safety pin (must set to True to overwrite and existing OUT_DIR)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( RUN_DIR );
  CommandArgument_String_OrDefault( OUT_DIR, "Reads.qual" );
  CommandArgument_String_OrDefault( HEAD, "reads" );
  CommandArgument_String_OrDefault( SELECT, "" );
  CommandArgument_UnsignedInt_OrDefault( BATCH_SIZE, 200000 );
  CommandArgument_Bool_OrDefault( NAMES, False );
  CommandArgument_Bool_OrDefault( OVERWRITE, False );
  EndCommandArguments;

  // Dir and file names.
  String quals_file = RUN_DIR + "/" + HEAD + ".qualb";
  String ids_file = RUN_DIR + "/" + HEAD + ".ids";
  String out_dir = RUN_DIR + "/" + OUT_DIR;
  String out_base = out_dir + "/" +  HEAD + "_";

  // Clean up from old runs.
  if ( IsDirectory( out_dir ) ) {
    if ( OVERWRITE ) {
      cout << Date( ) << ": removing quals from old run" << endl;
      String comm = "rm -rf " + out_dir + "/*.qual";
      System( comm );
    }
    else {
      cout << " " << OUT_DIR << " already exists: cannot proceed.\n" << endl;
      return 0;
    }
  }

  Mkdir777( out_dir );
  
  // Count reads and batches.
  int tot_count = MastervecFileObjectCount( quals_file );

  vec<int> sel;
  if ( SELECT == "" ) {
    sel.reserve( tot_count );
    for (int ii=0; ii<tot_count; ii++)
      sel.push_back( ii );
  }
  else {
    ParseIntSet( SELECT, sel );
    ForceAssert( is_sorted( sel.begin( ), sel.end( ) ) );
  }
  int n_reads = sel.size( );
  int n_batches = 1;
  while ( n_batches * (int)BATCH_SIZE < n_reads )
    n_batches++;

  // Loop over n_batches batches.
  for (int batch_id=0; batch_id<n_batches; batch_id++) {
    vec<int> read_ids;
    int begin_id = (int)BATCH_SIZE * batch_id;
    int end_id = Min( n_reads, (int)BATCH_SIZE * ( batch_id + 1 ) );
    int loc_batch_size = end_id - begin_id;
    read_ids.reserve( loc_batch_size );
    for (int ii=begin_id; ii<end_id; ii++)
      read_ids.push_back( sel[ii] );
    
    cout << Date( ) << ": batch "
	 << batch_id << ", of "
	 << loc_batch_size << " entries: ["
	 << begin_id << ", "
	 << end_id << ")" << endl;

    vecString ids;
    if ( NAMES )
      ids.SparseRead( ids_file, read_ids, 0 );
    
    vecqualvector quals;
    quals.SparseRead( quals_file, read_ids, 0 );
    
    String out_file = out_base + ToString( batch_id ) + ".qual";
    ofstream out( out_file.c_str( ) );
    for (int ii=begin_id; ii<end_id; ii++) {
      int rid = sel[ii];
      String read_name = NAMES ? ids[rid] : "read_" + ToString( rid );
      Print( out, quals[rid], read_name );
    }
    out.close( );
  }

}
