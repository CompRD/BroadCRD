/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "lookup/LookAlign.h"
#include "math/NStatsTools.h"
#include "paths/PdfEntry.h"
#include "paths/UnibaseUtils.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * PrepareUnibases
 *
 * Remove rc copies of unibases, and align what is left to a given
 * reference (if given). It will also save the estimated copy numbers
 * of the fw unibases.
 *
 * HEAD_UNI: head of unibases
 * OUT_DIR: where output will be saved
 * OUT_BASE: base name for output files, in OUT_DIR
 * K: kmer size (for unibases)
 * HEAD_REF: if given align fw unibases to reference (head of reference)
 * NUM_THREADS: needed by GetAlignsFast (if HEAD_REF is given)
 * KEEP_NXX: discard short unibases shorter than N<xx>
 * FORCE: arg to GetAlignsFast (if HEAD_REF is given)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( HEAD_UNI );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String( OUT_BASE );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_String_OrDefault( HEAD_REF, "" );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Int_OrDefault( KEEP_NXX, 98 );
  CommandArgument_Int_OrDefault( KEEP_SIZE, 1000 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  EndCommandArguments;

  // Dir and file names.
  String strK = ToString( K );
  String cnF = HEAD_UNI + ".unipaths.predicted_count.k" + strK;
  String unibasesF = HEAD_UNI + ".unibases.k" + strK;
  String lookupF = HEAD_REF + ".lookup";
 
  String out_head = OUT_DIR + "/" + OUT_BASE;
  String logF = out_head + ".PrepareUnibases.log";
  String uniF = out_head + ".fastb";
  String uni_cnsF = out_head + ".cns";
  String uni_pdfsF = out_head + ".pdfs";
  String orig_idsF = out_head + ".orig.ids";
  String uni_qltF = out_head + ".on_reference.qlt";

  vec<String> needed;
  needed.push_back( cnF );
  needed.push_back( unibasesF );
  if ( HEAD_REF != "" ) needed.push_back( lookupF );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( Dirname( out_head ) );

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Log stream.
  ofstream log( logF.c_str( ) );
  PrintCommandPretty( log );

  // Load fw unibases (generate them if needed).
  vec<int> uni_cns;
  VecPdfEntryVec uni_pdfs;
  vecbvec uni;
  if ( FORCE || ! ( IsRegularFile( uniF ) && IsRegularFile( orig_idsF ) ) ) {
    vec<int> orig_ids;

    cout << Date( ) << ": loading unibases" << endl;
    vecbvec unibases( unibasesF );
    
    cout << Date( ) << ": loading copy number estimates" << endl;
    VecPdfEntryVec cn_pdfs( cnF );
    
    cout << Date( ) << ": removing rc unibases" << endl;
    int n_keepers = 0;
    vec<bool> keepers( unibases.size( ), false );
    vec<int> to_rc;
    UnibaseInvolution( unibases, to_rc );
    for (int id1=0; id1<(int)unibases.size( ); id1++) {
      int id2 = to_rc[id1];
      if ( id1 <= id2 ) {
	keepers[id1] = true;
	n_keepers++;
      }
    }
    
    if ( KEEP_NXX < 100 || KEEP_SIZE > 0) {
      cout << Date( ) << ": discarding short unibases";
      int n_tot = n_keepers;
      
      // Compute unibases' Nxx.
      vec<int> lens;
      lens.reserve( n_keepers );
      for (int id=0; id<(int)unibases.size( ); id++)
	if ( keepers[id] )
	  lens.push_back( unibases[id].size( ) );
      
      vec<int> sel( 1, KEEP_NXX );
      vec<int> idx( 1, 0 );
      BasicNStats( lens, idx, &sel );
      ForceAssert( idx.size( ) == 1 );
      int Nxx = lens[idx[0]-1];
      cout << ", N" << KEEP_NXX << "=" << Nxx << endl;

      // Discard all unibases shorter than the Nxx size.
      for (int id=0; id<(int)unibases.size( ); id++) {
	if ( keepers[id] && (int)unibases[id].size( ) < Nxx
	     && (KEEP_SIZE == 0 || (int)unibases[id].size( ) < KEEP_SIZE)) {
	  keepers[id] = false;
	  n_keepers--;
	}
      }
      cout << Date( ) << ": kept " << n_keepers << " out of " << n_tot << endl;
    }
    
    int estim = -1;
    orig_ids.reserve( n_keepers );
    uni_cns.reserve( n_keepers );
    uni_pdfs.reserve( n_keepers );
    uni.reserve( n_keepers );
    for (int ii=0; ii<keepers.isize( ); ii++) {
      if ( ! keepers[ii] ) continue;
      GetMostLikelyValue( estim, cn_pdfs[ii] );
      orig_ids.push_back( ii );
      uni_cns.push_back( estim );
      uni_pdfs.push_back( cn_pdfs[ii] );
      uni.push_back( unibases[ii] );
    }
    
    cout << Date( ) << ": saving fw unibases" << endl;
    uni.WriteAll( uniF );
    WRITE( uni_cnsF, uni_cns );
    WRITE( orig_idsF, orig_ids );
    uni_pdfs.WriteAll( uni_pdfsF );
  }
  else {
    cout << Date( ) << ": loading fw unibases" << endl;
    uni.ReadAll( uniF );
    uni_pdfs.ReadAll( uni_pdfsF );
    READX( uni_cnsF, uni_cns );
  }

  // No reference, we are done.
  if ( HEAD_REF == "" ) {
    cout << Date( ) << ": no reference given - done" << endl;
    return 0;
  }

  // Align or load contigs to reference.
  vec<look_align> hits;
  GetAlignsFast( 96, uniF, lookupF, uni_qltF, hits, !FORCE, OUT_DIR ); 

  cout << Date( ) << ": sorting " << hits.size( ) << " aligns" << endl;
  order_lookalign_TargetBeginEnd sorter;
  sort( hits.begin( ), hits.end( ), sorter );
  
  // Log mapping.
  cout << Date( ) << ": logging alignments of unibases" << endl;
  log << "ALIGNMENTS OF UNIBASES ON REF, WITH ESTIMATED CN\n" << endl;
  for (size_t ii=0; ii<hits.size( ); ii++) {
    int estim_cn = uni_cns[ hits[ii].query_id ];
    double prob_cn = uni_pdfs[ hits[ii].query_id ][estim_cn].Prob( );

    log << hits[ii].query_id << "\t"
	<< hits[ii].query_length << "\tcn: "
	<< estim_cn << " ("
	<< ToString( 100. * prob_cn, 1 ) << "%)\t"
	<< ( hits[ii].rc1 ? "rc" : "fw" ) << " on "
	<< hits[ii].target_id << "\t["
	<< hits[ii].pos2( ) << ", "
	<< hits[ii].Pos2( ) << ")_"
	<< hits[ii].target_length << "\n";
					   }
  log << endl;
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

