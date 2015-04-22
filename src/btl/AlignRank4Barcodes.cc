///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

#include "Basevector.h"
#include "Intvector.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "btl/CBarcodes.h"
#include "btl/Minial.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"
#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/**
 * AlignRank4Barcodes
 *
 * Load full rank barcodes, and align the genomic portions to the
 * reference.
 */
int main( int argc, char* argv[] )
{
  RunTime( );

  BeginCommandArguments;

  // Core input.
  CommandArgument_String( REFERENCE_LOOKUP );
  CommandArgument_String( CAP_FASTB_FILE );
  CommandArgument_String( READS_HEAD );
  CommandArgument_String( BARCODES );
  CommandArgument_String( CAP_ALIGNS );
  CommandArgument_String( OUT_DIR );

  CommandArgument_Int_OrDefault( K, 26 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( FORCE, True );

  EndCommandArguments;

  // Check args.
  if ( K != 26 ) {
    cout << "Fatal error: for now only allowed value for K is 26.\n" << endl;
    return 1;
  }

  // File names.
  String reads_fastb_file = READS_HEAD + ".fastb";
  String pairs_file = READS_HEAD + ".pairs";

  vec<String> needed;
  needed.push_back( REFERENCE_LOOKUP );
  needed.push_back( CAP_FASTB_FILE );
  needed.push_back( reads_fastb_file );
  needed.push_back( pairs_file );
  needed.push_back( BARCODES );
  needed.push_back( CAP_ALIGNS );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  String log_file = OUT_DIR + "/main.log";
  String select_fastb_file = OUT_DIR + "/select.fastb";
  String aligns_file = OUT_DIR + "/on_reference.qlt";
  String barcodes_aligns_file = OUT_DIR + "/barcodes.qlt";

  Mkpath( OUT_DIR );

  // Thread control (needed by FastAlignShortReads).
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Log stream.
  cout << "Sending log to " << log_file << endl;
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );

  // Load full rank barcodes.
  log << Date( ) << ": loading barcodes... " << flush;
  vec<CBarcodes> barcodes;
  LoadVecBarcodes( BARCODES, barcodes );
  log << ToStringAddCommas( barcodes.size( ) ) << " barcodes found" << endl;

  log << Date( ) << ": selecting full rank barcodes... " << flush;
  {
    vec<CBarcodes> select;
    select.reserve( barcodes.size( ) );
    for (size_t ii=0; ii<barcodes.size( ); ii++)
      if ( barcodes[ii].Rank( ) >= 4 )
	select.push_back( barcodes[ii] );
    swap( barcodes, select );
  }
  log << ToStringAddCommas( barcodes.size( ) ) << " selected" << endl;
  
  // Load pairs and cap aligns.
  size_t n_reads = MastervecFileObjectCount( reads_fastb_file );
  
  vecbvec cap_bases( CAP_FASTB_FILE );
  ForceAssertEq( (int)cap_bases.size( ), 1 );
  const int cap_len = cap_bases[0].size( );

  log << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  
  log << Date( ) << ": loading cap aligns" << endl;
  READ( CAP_ALIGNS, vec<minial>, caps );

  // Generate a map read_id to cap_align_id;
  log << Date( ) << ": building reads to cap aligns map" << endl;
  vec<int64_t> caps_idx( n_reads, -1 );
  for (int64_t ii=0; ii<(int64_t)caps.size( ); ii++) {
    int64_t rid = caps[ii].ReadId( );
    if ( caps_idx[rid] == -1 ) caps_idx[rid] = ii;
    else caps_idx[rid] = -2;
  }
  
  // Collect ids of reads to be aligned (SparseRead wants int ids).
  log << Date( ) << ": seleting ids of reads for alignment" << endl;
  vec<int> rids;
  vec<int> seglens;   // note that -1 means "all reads" (no cap found)
  rids.reserve( barcodes.size( ) );
  for (size_t ii=0; ii<barcodes.size( ); ii++) {
    const CBarcodes &bar = barcodes[ii];
    const size_t bc_pairs = bar.NPairs( );
    for (size_t jj=0; jj<bc_pairs; jj++) {
      size_t pid = bar.PairId( jj );
      size_t id1 = pairs.ID1( pid );
      size_t id2 = pairs.ID2( pid );
      
      const minial *al1 = ( caps_idx[id1] > -1 ) ? &caps[ caps_idx[id1] ] : 0;
      const minial *al2 = ( caps_idx[id2] > -1 ) ? &caps[ caps_idx[id2] ] : 0;
      const bool fwcap_1 = ( al1 && ! al1->Rc( ) );
      const bool fwcap_2 = ( al2 && ! al2->Rc( ) );
      if ( ! ( fwcap_1 || fwcap_2 ) ) continue;   // this should not happen
      
      const minial *sel_al = ( fwcap_1 ) ? al2 : al1;
      const int sel_id = ( fwcap_1 ) ? id2 : id1;
      const int seglen = sel_al ? sel_al->offset_ : -1;
      if ( seglen != -1 && seglen < K ) continue;

      rids.push_back( sel_id );
      seglens.push_back( seglen );
    }
  }
  
  // Load selected reads, chop when needed, and save to file.
  log << Date( ) << ": loading selected reads" << endl;
  vecbvec reads;
  reads.SparseRead( reads_fastb_file, rids, 0 );

  log << Date( ) << ": removing cap from reads" << endl;
  for (size_t ii=0; ii<rids.size( ); ii++) {
    const int rid = rids[ii];
    const int start = seglens[ii];
    if ( start < 0 ) continue;   // no cap found, no need to curtail

    bvec chunk( reads[rid], 0, start );
    reads[rid] = chunk;
  }

  log << Date( ) << ": saving selected reads" << endl;
  reads.WriteAll( select_fastb_file );
  
  // Align reads.
  log << Date( ) << ": aligning reads" << endl;
  vec<look_align> aligns;
  GetAlignsFast( K, select_fastb_file, REFERENCE_LOOKUP,
		 aligns_file, aligns, ! FORCE, OUT_DIR );
  
  // Map reads to aligns.
  log << Date( ) << ": mapping reads to aligns" << endl;
  vec<int> to_align( n_reads, -1 );
  for (int ii=0; ii<aligns.isize( ); ii++) {
    int rid = aligns[ii].query_id;
    if ( to_align[rid] == -1 ) to_align[rid] = ii;
    else to_align[rid] = -2;
  }

  // Output stream for aligns, split by matching barcodes.
  ofstream aout( barcodes_aligns_file.c_str( ) );

  // Loop over all barcodes.
  log << Date( ) << ": dumping aligns by barcodes" << endl;
  for (size_t ii=0; ii<barcodes.size( ); ii++) {
    const CBarcodes &bar = barcodes[ii];
    const size_t bc_pairs = bar.NPairs( );

    vec<look_align> bar_aligns;
    bar_aligns.reserve( bc_pairs );
    for (size_t jj=0; jj<bc_pairs; jj++) {
      size_t pid = bar.PairId( jj );
      int idx1 = to_align[ pairs.ID1( pid ) ];
      int idx2 = to_align[ pairs.ID2( pid ) ];
      if ( idx1 > -1 ) bar_aligns.push_back( aligns[idx1] );
      if ( idx2 > -1 ) bar_aligns.push_back( aligns[idx2] );
    }
    order_lookalign_TargetBegin sorter;
    sort( bar_aligns.begin( ), bar_aligns.end( ) );
    
    bar.PrintOneLine( aout );
    for (int jj=0; jj<bar_aligns.isize( ); jj++)
      bar_aligns[jj].PrintParseable( aout );
    aout << "\n";
  }  

  aout.close( );

  // Done.
  String date = "\n" + Date( );
  cout << date << ": AlignRank4Barcodes done" << endl;
  log << date << ": AlignRank4Barcodes done" << endl;
  log.close( );
  
}
