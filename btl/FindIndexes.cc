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
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/KmerAligner.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/FindVector.h"
#include "util/RunCommand.h"

/**
 * IsLowQuality
 *
 * Decide if a read is low quality.
 */
bool IsLowQuality( const double &LOW_QUAL_THRESHOLD,
		   const int &LOW_QUAL,
		   const qvec &qual )
{
  int n_bad = 0;
  int n_len = qual.size( );
  for (int ii=0; ii<n_len; ii++)
    if ( (int)qual[ii] < LOW_QUAL ) n_bad++;
  
  return ( SafeQuotient( n_bad, n_len ) > LOW_QUAL_THRESHOLD );
}

/**
 * HitsCount
 */
int HitsCount( const vec<minial> &minis,
	       const vec<int64_t> &firstmini,
	       const int64_t &rid,
	       const int *vid = 0 )
{
  int count = 0;
  if ( firstmini[rid] < 0 ) return count;
  
  for (size_t ii=firstmini[rid]; ii<minis.size( ); ii++) {
    if ( minis[ii].ReadId( ) != rid ) break;
    if ( vid && minis[ii].vector_id_ != *vid ) continue;
    count++;
  }
  
  return count;
}

/**
 * FindIndexesCore
 *
 * Note that fw_aligners and rc_aligners may not be passed as const
 * args, because FindPossibleAligments is not a const method.
 *
 * al: align of cap to read (if not null).
 * raw_barcode: ids if indexes on read (-1 if not found)
 * expected_rank: analogue to IntVecRank in CBarcodes.h
 */
void FindIndexesCore( const int BARCODES_LEN,
		      const int PRIMER_LEN,
		      const int PRIMER_SLACK,
		      const int MAX_MISMATCHES,
		      const int MAX_INDELS,
		      const int MIN_ALIGN_LEN,
		      const vec< vecbvec > &fw_barcodes,
		      const vec< vecbvec > &rc_barcodes,
		      vec< KmerAligner<8> > &fw_aligners,
		      vec< KmerAligner<8> > &rc_aligners,
		      const bvec &read,
		      const minial *al,
		      IntVec &raw_barcode,
		      int &expected_rank )
{
  // Default.
  raw_barcode.clear( );
  raw_barcode.resize( 4, -1 );
  
  // Cap not found, nothing to do.
  if ( ! al ) return;

  // Consts.
  const bool cap_rc = al->Rc( );
  const int seg_len = BARCODES_LEN + PRIMER_LEN;
  
  // Expected offset for the four indexes (without slack).
  vec<int> expected( 4, 0 );
  vec<char> indext( 4, ' ');
  for (int ii=0; ii<4; ii++) {
    indext[ii] = ( ii==0 ? 'A' : ( ii==1 ? 'B' : ( ii==2 ? 'C' : 'D' ) ) );
    if ( cap_rc ) {
      int cap_end = (int)read.size( ) - al->offset_;
      expected[ii] = cap_end + PRIMER_LEN + ( ii * seg_len );
    }
    else {
      int cap_start = al->offset_;
      expected[ii] = cap_start - ( ( ii + 1 ) * seg_len );
    }
  }
  
  // Expected rank.
  int rlen = read.size( );
  expected_rank = 0;
  for (int ii=0; ii<4; ii++) {
    if ( expected[ii] >=0 && expected[ii] + BARCODES_LEN <= rlen )
      expected_rank = ii+1;
    else
      break;
  }
  
  // Loop over all four indexes (A, B, C, and D).
  for (int type_id=0; type_id<4; type_id++) {
    vec<int> ids;
    vec<int> pos;

    // This is not a bug. The lab routinely calls "fw" the barcodes
    //  that are actually rc, and "rc" the ones that are actually fw.
    if ( cap_rc ) fw_aligners[type_id].FindPossibleAlignments( read, pos, ids );
    else rc_aligners[type_id].FindPossibleAlignments( read, pos, ids );

    // Band for Smith-Waterman.
    const int sw_band = ( type_id + 1 ) * PRIMER_SLACK;
    
    // Clean and compactify placements.
    vec<int> select;
    select.reserve( ids.size( ) );
    for (int jj=0; jj<ids.isize( ); jj++) {
      int off = -pos[jj];
      if ( Abs( off - expected[type_id] ) > sw_band ) continue;
      
      select.push_back( ids[jj] );
    }

    sort( select.begin( ), select.end( ) );
    select.erase( unique( select.begin( ), select.end( ) ), select.end( ) );
    
    // helper[jj]: errors, barcode_id.
    vec< pair<int,int> > helper;

    // Align candidates.
    for (int jj=0; jj<select.isize( ); jj++) {
      const int bid = select[jj];

      // Not a bug, see comment on fw and rc barcodes above.
      const bvec &bar
	= cap_rc
	? fw_barcodes[type_id][bid]
	: rc_barcodes[type_id][bid];
      
      // Align barcode to read.
      align al;
      int err = 0;
      SmithWatBandedA( read, bar, expected[type_id], sw_band, al, err );
      vector<int> mutgaps = al.MutationsGap1Gap2( read, bar );
      
      if ( mutgaps[0] > MAX_MISMATCHES ) continue;
      if ( mutgaps[1] + mutgaps[2] > MAX_INDELS ) continue;
      if ( al.Pos1( ) - al.pos1( ) < MIN_ALIGN_LEN ) continue;
      
      int n_errs = mutgaps[0] + mutgaps[1] + mutgaps[2];
      helper.push_back( make_pair( n_errs, bid ) );
    }

    // No align found.
    if ( helper.size( ) < 1 ) continue;

    // A tie (bail out).
    sort( helper.begin( ), helper.end( ) );
    if ( helper.size( ) > 1 && helper[0].first == helper[1].first ) continue;

    // A winner.
    raw_barcode[type_id] = helper[0].second;
    
  } // loop over all four indexes
  
}

/**
 * FindIndexes
 *
 * Analyze an Indexed Jump Library.  The pipeline consists of these
 * basic steps: (1) find cap on reads; (2) place indexes on reads; and (3)
 * build pairs from full-rank fragment pairs.
 *
 * REMARKS
 *
 * 1) orientation of indexes vs cap. If the cap is fw on read, then
 *    the indexes are rc on read (and vice versa), as this:
 *
 *          D    C    B    A   CAP
 *        <--- <--- <--- <--- ---->
 *        -------------------------------------->
 *                                        read
 *
 *    Note in the picture the "gap" between cap and index A, and
 *    between adjacent indexes, correspondig the four bases primer
 *    before the index.
 */
int main( int argc, char* argv[] )
{
  RunTime( );

  BeginCommandArguments;

  // Core input.
  CommandArgument_String( VECTORS_DIR );
  CommandArgument_String( READS_HEAD );

  // Output dir.
  CommandArgument_String_OrDefault( OUT_DIR, READS_HEAD + ".FindIndexes" );

  // Size of primers and barcodes, with some slack.
  CommandArgument_UnsignedInt_OrDefault( BARCODES_LEN, 20 );
  CommandArgument_UnsignedInt_OrDefault( PRIMER_LEN, 4 );
  CommandArgument_UnsignedInt_OrDefault( PRIMER_SLACK, 1 );

  // Discard pairs if one (or both) of the read has too many low qual bases.
  CommandArgument_Double_OrDefault( LOW_QUAL_THRESHOLD, 0.75 );
  CommandArgument_Int_OrDefault( LOW_QUAL, 20 );

  // Set FULL_RANK to 2 (or 3), if only using barcodes A, B (and C).
  CommandArgument_Int_OrDefault( FULL_RANK, 4 );
  
  // prioritize using the genomic portion of barcode read
  CommandArgument_Bool_OrDefault( PRIORITIZE_BARCODE_GDNA, False );

  // Discard pairs if genomic portion is too short.
  CommandArgument_Int_OrDefault( MIN_GENLEN, 24 );
  
  // Args for various aligners.
  CommandArgument_UnsignedInt_OrDefault( MAX_MISMATCHES, 4 );
  CommandArgument_UnsignedInt_OrDefault( MAX_INDELS, 1 );
  CommandArgument_Double_OrDefault( MAX_ER, 0.01 );
  CommandArgument_UnsignedInt_OrDefault( BANDWIDTH, 2 );
  
  // Args for KmerAligner, and for the align of barcodes on reads.
  CommandArgument_UnsignedInt_OrDefault( KA_STEP, 1 );
  CommandArgument_UnsignedInt_OrDefault( MIN_ALIGN_LEN, 15 );
  
  // Do not use cached data, if FORCE=True.
  CommandArgument_Bool_OrDefault( FORCE, False );
  

  // Show progress in log.
  CommandArgument_Int_OrDefault( DOTTER, 1000000 );

  // do PCR removal or not
  CommandArgument_Bool_OrDefault( DISCARD_DUPLICATES, True );

  EndCommandArguments;
  
  // Dir and file names.
  String cap_fastb_file = VECTORS_DIR + "/cap.fastb";

  String barcodesA_fastb_file = VECTORS_DIR + "/barcodesA.fastb";
  String barcodesB_fastb_file = VECTORS_DIR + "/barcodesB.fastb";
  String barcodesC_fastb_file = VECTORS_DIR + "/barcodesC.fastb";
  String barcodesD_fastb_file = VECTORS_DIR + "/barcodesD.fastb";

  String barcodesA_ids_file = VECTORS_DIR + "/barcodesA.ids";
  String barcodesB_ids_file = VECTORS_DIR + "/barcodesB.ids";
  String barcodesC_ids_file = VECTORS_DIR + "/barcodesC.ids";
  String barcodesD_ids_file = VECTORS_DIR + "/barcodesD.ids";
  
  String reads_fastb_file = READS_HEAD + ".fastb";
  String reads_qualb_file = READS_HEAD + ".qualb";
  String reads_pairs_file = READS_HEAD + ".pairs";

  String log_file = OUT_DIR + "/main.log";
  String cap_aligns_file = OUT_DIR + "/cap.aligns.qlt";
  String cap_minials_file = OUT_DIR + "/cap.aligns";
  String barcodes_file = OUT_DIR + "/barcodes";
  String deleted_ids_file = OUT_DIR + "/deleted.ids";
  String out_head = OUT_DIR + "/indexed_reads";

  vec<String> needed;
  needed.push_back( cap_fastb_file );

  needed.push_back( barcodesA_fastb_file );
  needed.push_back( barcodesB_fastb_file );
  needed.push_back( barcodesC_fastb_file );
  needed.push_back( barcodesD_fastb_file );
  needed.push_back( barcodesA_ids_file );
  needed.push_back( barcodesB_ids_file );
  needed.push_back( barcodesC_ids_file );
  needed.push_back( barcodesD_ids_file );

  needed.push_back( reads_fastb_file );
  needed.push_back( reads_qualb_file );
  needed.push_back( reads_pairs_file );

  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  if ( FULL_RANK < 2 || FULL_RANK > 4 ) {
    cout << "Fatal error: FULL_RANK can only be set to 2, 3, or 4.\n" << endl;
    return 1;
  }

  Mkpath( OUT_DIR );
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << "\nSending log to " << log_file << "\n" << endl;
  
  // Keep track of why we remove pairs.
  vec<int> deleted;

  // Load quals (scoped for memory ), and pairs, and remove poor quality pairs.
  log << Date( ) << ": loading pairs... " << flush;
  PairsManager pairs( reads_pairs_file );
  const size_t n_pairs = pairs.nPairs( );
  log << ToStringAddCommas( n_pairs ) << " pairs found" << endl;

  {
    deleted.resize( n_pairs, NOT_DELETED );
    
    log << Date( ) << ": loading quals" << endl;
    vecqvec quals( reads_qualb_file );

    log << Date( ) << ": removing low qual pairs (. = "
	<< ToStringAddCommas( DOTTER ) << " pairs)"
	<< endl;
    for (size_t ii=0; ii<n_pairs; ii++) {
      if ( ii % DOTTER == 0 ) Dot( log, ii / DOTTER );

      size_t id1 = pairs.ID1( ii );
      size_t id2 = pairs.ID2( ii );
      if ( IsLowQuality( LOW_QUAL_THRESHOLD, LOW_QUAL, quals[id1] ) ||
	   IsLowQuality( LOW_QUAL_THRESHOLD, LOW_QUAL, quals[id2] ) )
	deleted[ii] = CAP_LOW_QUAL;
    }
    log << endl;
  }

  // Load reads, and place cap on reads.
  log << Date( ) << ": loading reads" << endl; 
  vecbvec cap_bases( cap_fastb_file );
  vecbvec reads( reads_fastb_file );
  const size_t n_reads = reads.size( );
  ForceAssertEq( (int)cap_bases.size( ), 1 );
  
  if ( FORCE || ! IsRegularFile( cap_minials_file ) ) {
    log << Date( ) << ": running FindVector to place cap" << endl;
    vec<look_align> hits;
    FindVector( cap_aligns_file, cap_bases, reads,
		MAX_MISMATCHES, MAX_INDELS, MAX_ER, False, &log, &hits );
    
    log << Date( ) << ": converting look_aligns to minials" << endl;
    int n_improper = 0;
    vec<minial> minis;
    minis.reserve( hits.size( ) );
    for (size_t ii=0; ii<hits.size( ); ii++)
      if ( hits[ii].IsProper( ) ) minis.push_back( minial( hits[ii] ) );
      else n_improper++;
    size_t nminis = minis.size( );
    String str_nminis = ToStringAddCommas( nminis );

    log << Date( ) << ": sorting " << str_nminis << " minials" << endl;
    if ( ! is_sorted( minis.begin( ), minis.end( ) ) )
      sort( minis.begin( ), minis.end( ) );
    
    log << Date( ) << ": removing spurious minials... " << flush;
    minis.erase( unique( minis.begin( ), minis.end( ) ), minis.end( ) );
    size_t ndeleted = nminis - minis.size( );
    log << ndeleted << " removed" << endl;

    if ( n_improper > 0 )
      log << "\nWARNING: " << ToString( n_improper )
	  << " improper aligns were discarded\n" << endl;

    log << Date( ) << ": saving minials" << endl;
    WRITE( cap_minials_file, minis );
  }

  log << Date( ) << ": loading cap minials" << endl;
  READ( cap_minials_file, vec<minial>, minis );
  if ( ! is_sorted( minis.begin( ), minis.end( ) ) )
    sort( minis.begin( ), minis.end( ) );

  const size_t n_minis = minis.size( );
  vec<int64_t> firstmini( n_reads, -1 );
  for (int64_t ii=n_minis-1; ii>=0; ii--)
    firstmini[minis[ii].ReadId( )] = ii;

  // Load barcodes, and create KmerAligners.
  log << Date( ) << ": building KmerAligners" << endl;

  vec< vecbvec > fw_barcodes( 4 );
  vec< vecbvec > rc_barcodes( 4 );
  vec< KmerAligner<8> > fw_aligners( 4 );
  vec< KmerAligner<8> > rc_aligners( 4 );
  for (int ii=0; ii<4; ii++) {
    const String *in_file = 0;
    if ( ii == 0 ) in_file = &barcodesA_fastb_file;
    else if ( ii == 1 ) in_file = &barcodesB_fastb_file;
    else if ( ii == 2 ) in_file = &barcodesC_fastb_file;
    else in_file = &barcodesD_fastb_file;

    fw_barcodes[ii].ReadAll( *in_file );
    for (size_t jj=0; jj<fw_barcodes[ii].size( ); jj++)
      ForceAssertEq( fw_barcodes[ii][jj].size( ), BARCODES_LEN );
    
    rc_barcodes[ii] = fw_barcodes[ii];
    for (size_t jj=0; jj<rc_barcodes[ii].size( ); jj++)
      rc_barcodes[ii][jj].ReverseComplement( );
    
    fw_aligners[ii].SetKmerStep( KA_STEP );
    fw_aligners[ii].SetBases( fw_barcodes[ii] );
    
    rc_aligners[ii].SetKmerStep( KA_STEP );
    rc_aligners[ii].SetBases( rc_barcodes[ii] );
  }
  
  // Raw barcodes and expected ranks (in sync with pair ids).
  VecIntVec raw_barcodes( n_pairs, IntVec( 4, -1 ) );
  IntVec exp_ranks( n_pairs, -1 );
  
  // Loop over pairs.
  log << Date( ) << ": aligning barcodes (. = "
      << ToStringAddCommas( DOTTER ) << " pairs)"
      << endl;
  for (size_t pid=0; pid<n_pairs; pid++) {
    if ( pid % DOTTER == 0 ) Dot( log, pid / DOTTER );

    if ( deleted[pid] != NOT_DELETED ) continue;

    int64_t id1 = pairs.ID1( pid );
    int64_t id2 = pairs.ID2( pid );
    int n1 = HitsCount( minis, firstmini, id1 );
    int n2 = HitsCount( minis, firstmini, id2 );

    if ( n1 < 1 && n2 < 1 ) {
      deleted[pid] = CAP_UNPLACED;
      continue;
    }
    
    if ( n1 > 1 || n2 > 1 ) {
      deleted[pid] = CAP_MULT;
      continue;
    }

    const minial *al1 = n1 > 0 ? &minis[ firstmini[id1] ] : 0;
    const minial *al2 = n2 > 0 ? &minis[ firstmini[id2] ] : 0;
    if ( ( al1 && al2 ) && ( al1->Rc( ) == al2->Rc( ) ) ) {
      deleted[pid] = CAP_ILLOGICAL;
      continue;
    }
    
    IntVec index1;
    int exprank1;
    FindIndexesCore( BARCODES_LEN, PRIMER_LEN, PRIMER_SLACK,
		     MAX_MISMATCHES, MAX_INDELS, MIN_ALIGN_LEN,
		     fw_barcodes, rc_barcodes, fw_aligners, rc_aligners,
		     reads[id1], al1, index1, exprank1 );
    
    IntVec index2;
    int exprank2;
    FindIndexesCore( BARCODES_LEN, PRIMER_LEN, PRIMER_SLACK,
		     MAX_MISMATCHES, MAX_INDELS, MIN_ALIGN_LEN,
		     fw_barcodes, rc_barcodes, fw_aligners, rc_aligners,
		     reads[id2], al2, index2, exprank2 );
    
    if ( ! IntVecsConsistent( index1, index2 ) ) {
      deleted[pid] = IDX_MISMATCH;
      continue;
    }

    int rank1 = IntVecRank( index1 );
    int rank2 = IntVecRank( index2 );
    if ( rank2 > rank1 ) {
      raw_barcodes[pid] = index2;
      exp_ranks[pid] = exprank2;
    }
    else {
      raw_barcodes[pid] = index1;
      exp_ranks[pid] = exprank1;
    }
  }
  log << endl;

  // Package into a vec of CBarcodes.
  log << Date( ) << ": building CBarcodes" << endl;
  vec<CBarcodes> barcodes;
  BuildVecBarcodes( raw_barcodes, exp_ranks, barcodes );

  // Discard not full-rank barcodes.
  log << Date( ) << ": discarding non null-rank pairs" << endl;
  DiscardNonFullRank( FULL_RANK, raw_barcodes, exp_ranks, deleted );


    vec<short> vCapEnds(reads.size(),0);
    if(PRIORITIZE_BARCODE_GDNA){
        SetCapEnds(reads,cap_bases[0],minis,vCapEnds);
    }
        
  // Discard all pairs for which the genomic portion is too small.
  log << Date( ) << ": discarding pairs with short genomic lengths" << endl;
  vec<int> klens;
  DiscardShortGenomic( MIN_GENLEN, pairs, reads, minis, barcodes, vCapEnds,
		       klens, deleted );

  // Load selected quals.
  log << Date( ) << ": loading quals" << endl;
  vecqvec quals;
  {
    vec<size_t> select;
    select.reserve( reads.size( ) );
    
    for (size_t barid=0; barid<barcodes.size( ); barid++) {
      const CBarcodes &bar = barcodes[barid];
      const int n_pairs = bar.NPairs( );
      
      for (int ii=0; ii<n_pairs; ii++) {
	const size_t pid = bar.PairId( ii );
	if ( deleted[pid] ) continue;
	
	const size_t id1 = pairs.ID1( pid );
	const size_t id2 = pairs.ID2( pid );
	ForceAssert( vCapEnds[id1] != 0 || vCapEnds[id2] != 0 || klens[id1] > 0 || klens[id2] > 0 );
	
	if ( vCapEnds[id1] != 0 ) select.push_back( id1 );
	else if ( vCapEnds[id2] != 0 ) select.push_back( id2 );
	else if ( klens[id1] > 0 ) select.push_back( id1 );
	else select.push_back( id2 );
      }
    }
    sort( select.begin( ), select.end( ) );
    
    quals.SparseRead( reads_qualb_file, select, 0 );
  }  

  // DiscardDuplicates
  log << Date( ) << ": removing duplicates" << endl;
  if(DISCARD_DUPLICATES) DiscardDuplicates( pairs, reads, quals, barcodes, vCapEnds, klens, deleted );
  
  // Save.
  log << Date( ) << ": saving barcodes" << endl;
  WriteVecBarcodes( barcodes_file, barcodes );
  
  log << Date( ) << ": saving deleted.ids file" << endl;
  {
    ofstream out( deleted_ids_file.c_str( ) );
    for( size_t ii=0; ii<deleted.size( ); ii++)
      out << (int)deleted[ii] << "\n";
    out.close( );
  }
  
  log << Date( ) << ": saving genomic chunks" << endl;
  SaveIndexedFasta( pairs, reads, quals, out_head, barcodes, vCapEnds,klens, deleted);
  
  // Build report.
  ReportDeleted( deleted, pairs, log );

  // Done.
  String date = Date( );
  cout << date << ": FindIndexes done" << endl;
  log << date << ": FindIndexes done" << endl;
  log.close( );

}
