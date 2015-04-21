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
#include "paths/PdfEntry.h"
#include "paths/UnibaseUtils.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * UniChains
 *
 * Generate chains of unibases. These are linear sequences of
 * unibases, where two unibases adjacent in the sequence are either
 * linked by a certain number of links, or close to each other on a
 * given reference (or set of references).
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( HEAD_UNI );
  CommandArgument_String( HEAD_REF );
  CommandArgument_String( OUT_DIR );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  EndCommandArguments;

  // Dir and file names.
  String strK = ToString( K );
  String cnF = HEAD_UNI + ".unipaths.predicted_count.k" + strK;
  String unibasesF = HEAD_UNI + ".unibases.k" + strK;
  String lookupF = HEAD_REF + ".lookup";
 
  String uniF = OUT_DIR + "/uni.k" + strK + ".fastb";
  String uni_cnsF = OUT_DIR + "/uni.cns";
  String orig_idsF = OUT_DIR + "/orig_unibases.ids";
  String uni_qltF = OUT_DIR + "/uni_on_reference.qlt";

  vec<String> needed;
  needed.push_back( cnF );
  needed.push_back( unibasesF );
  needed.push_back( lookupF );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( OUT_DIR );

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Load fw unibases (generate them if needed).
  vec<int> uni_cns;
  vecbvec uni;
  if ( ! ( IsRegularFile( uniF ) && IsRegularFile( orig_idsF ) ) ) {
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
    
    int estim = -1;
    orig_ids.reserve( n_keepers );
    uni_cns.reserve( n_keepers );
    uni.reserve( n_keepers );
    for (int ii=0; ii<keepers.isize( ); ii++) {
      if ( ! keepers[ii] ) continue;
      GetMostLikelyValue( estim, cn_pdfs[ii] );
      orig_ids.push_back( ii );
      uni_cns.push_back( estim );
      uni.push_back( unibases[ii] );
    }
    
    cout << Date( ) << ": saving fw unibases" << endl;
    uni.WriteAll( uniF );
    WRITE( uni_cnsF, uni_cns );
    WRITE( orig_idsF, orig_ids );
  }
  else {
    cout << Date( ) << ": loading fw unibases" << endl;
    uni.ReadAll( uniF );
    READX( uni_cnsF, uni_cns );
  }

  // Align or load contigs to reference (logging within).
  vec<look_align> hits;
  GetAlignsFast( 96, uniF, lookupF, uni_qltF, hits, !FORCE, OUT_DIR ); 

  cout << Date( ) << ": sorting " << hits.size( ) << " aligns" << endl;
  order_lookalign_TargetBeginEnd sorter;
  sort( hits.begin( ), hits.end( ), sorter );
  
  // Print info.
  cout << "\n";
  for (size_t ii=0; ii<hits.size( ); ii++)
    if ( uni_cns[ii] < 2 )
      cout << hits[ii].query_id << "\t"
	   << hits[ii].query_length << "\t"
	   << ( hits[ii].rc1 ? "rc" : "fw" ) << " on "
	   << hits[ii].target_id << "\t["
	   << hits[ii].pos2( ) << ", "
	   << hits[ii].Pos2( ) << ")_"
	   << hits[ii].target_length << "\n";
  cout << endl;
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

