/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/CRefMerger.h"
#include "paths/SaveScaffoldGraph.h"
#include "paths/UnibaseUtils.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * MergeContigsOnReference
 *
 * Align contigs to a reference, and merge consecutive contigs, if
 * they overlap by >= MIN_OVERLAP bases (perfect matches only). If a
 * repetitive contig is involved in multiple overlaps, all are taken
 * into account. For example, if unique sequences A and B both align
 * repetitive sequence R, then we merge both A + R -> A', and B + R ->
 * B'.
 *
 * K: needed by the aligner (GetAlignsFast)
 * CONTIGS: full path name of input fastb
 * REF_HEAD: it loads <REF_HEAD>.{fastb,lookup}
 * ASSEMBLY_OUT: full path name of output assembly head
 * UNIBASES_K: if not empty, input contigs are unibases (remove rc copies)
 * MIN_CLEN: only save contigs >= MIN_CLEN
 * MIN_OVERLAP: minimum (perfect) overlap required for merging
 * SWBAND_RATIO: sw band, defined as ( overlap / SWBAND_RATIO )
 * MAX_GAP: max gap size allowed
 * MIN_GAP: min gap size allowed
 * MIN_GAP_DEV: min value for gaps' dev
 * NUM_THREADS: use all available if 0
 * FW_ONLY: discard contigs that align rc on reference
 * FORCE: do not load cached aligns, regenerate them
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( CONTIGS );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String( ASSEMBLY_OUT );
  CommandArgument_String_OrDefault( UNIBASES_K, "" );
  CommandArgument_Int_OrDefault( MIN_CLEN, 1000 );
  CommandArgument_Int_OrDefault( MIN_OVERLAP, 12 );
  CommandArgument_Int_OrDefault( SWBAND_RATIO, 6 );
  CommandArgument_Int_OrDefault( MAX_GAP, 20000 );
  CommandArgument_Int_OrDefault( MIN_GAP, -10000 );
  CommandArgument_Int_OrDefault( MIN_GAP_DEV, 20 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0);
  CommandArgument_Bool_OrDefault( FW_ONLY, True );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // WARNING! The code is temporary broken if UNIBASES_K is not defined.
  if ( UNIBASES_K == "" ) {
    cout << "FATAL ERROR - At this time the argument UNIBASES_K must be\n"
	 << "given: CRerfManager does not currently support rc aligns.\n"
	 << "Notice that contigs could own both fw and rc aligns, which means\n"
	 << "that we cannot just flip contigs if they align rc.\n"
	 << "\n"
	 << "LEAVING NOW.\n"
	 << endl;
    return 1;
  }

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Dir and file names.
  String tmp_dir = ASSEMBLY_OUT + ".tmp";

  String alignsFile = tmp_dir + "/aligns.qlt";
  String lookupFile = REF_HEAD + ".lookup";
  String targetFile = REF_HEAD + ".fastb";
  
  // Needed.
  vec<String> needed;
  needed.push_back( CONTIGS );
  needed.push_back( lookupFile );
  needed.push_back( targetFile );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( tmp_dir );
  
  // Load.
  vec<look_align> aligns;
  GetAlignsFast( K, CONTIGS, lookupFile, alignsFile, aligns, !FORCE, tmp_dir );

  cout << Date( ) << ": loading contigs" << endl;
  vecbvec contigs( CONTIGS );

  cout << Date( ) << ": loading reference" << endl;
  vecbvec targets( targetFile );
  
  // Deal with rc copies (only if input consists of unibases).
  if ( UNIBASES_K != "" ) {
    cout << Date( ) << ": unibases in input, removing rc copies" << endl;

    // If an unibase has at least a fw align, then it is a keeper.
    vec<bool> keepers( contigs.size( ), false );
    vec<bool> aligned( contigs.size( ), false );
    for (size_t ii=0; ii<aligns.size( ); ii++) {
      if ( aligns[ii].Fw1( ) ) keepers[ aligns[ii].query_id] = True;
      aligned[ aligns[ii].query_id] = True;
    }

    // Delete an unibase if it is not a keeper, but its rc is.
    vec<bool> deleters( contigs.size( ), false );
    vec<int> to_rc;
    UnibaseInvolution( contigs, to_rc );
    for (size_t ii=0; ii<aligns.size( ); ii++) {
      int cid = aligns[ii].query_id;
      if ( ( ! keepers[cid] ) && ( keepers[ to_rc[cid] ] ) )
	deleters[cid] = true;
    }
    
    // If neither the unibase nor its rc are aligned, remove one of them.
    for (int ii=0; ii<to_rc.isize( ); ii++) {
      if ( ii >= to_rc[ii] || aligned[ii] || aligned[ to_rc[ii] ] ) continue;
      deleters[ii] = true;
    }

    // Clear deleters contigs.
    for (size_t ii=0; ii<contigs.size( ); ii++)
      if ( deleters[ii] )
	contigs[ii].resize( 0 );
    
    // Remove rc aligns, and aligns of deleted contigs.
    vec<look_align> select;
    select.reserve( aligns.size( ) );
    for (size_t ii=0; ii<aligns.size( ); ii++) {
      if ( deleters[ aligns[ii].query_id ] ) continue;
      if ( ! aligns[ii].IsProper( ) ) continue;
      if ( aligns[ii].Rc1( ) ) continue;
      select.push_back( aligns[ii] );
    }
    swap( select, aligns );
  }
  
  // In any case have to remove rc and improper aligns.
  else {
    cout << Date( ) << ": filtering aligns" << endl;
    vec<look_align> select;
    select.reserve( aligns.size( ) );
    for (size_t ii=0; ii<aligns.size( ); ii++) {
      if ( ! aligns[ii].IsProper( ) ) continue;
      if ( aligns[ii].Rc1( ) ) continue;
      select.push_back( aligns[ii] );
    }
    swap( select, aligns );
  }

  // Merge contigs.
  CRefMerger merger( MIN_OVERLAP, SWBAND_RATIO, MAX_GAP, MIN_GAP, MIN_GAP_DEV,
		     targets, contigs, aligns );
  merger.Merge( &cout );
  merger.Save( ASSEMBLY_OUT, MIN_CLEN, &cout );
  
  // Done.
  cout << Date( ) << ": MergeContigsOnReference done" << endl;
  
}
