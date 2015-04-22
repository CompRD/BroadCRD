/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "feudal/IncrementalWriter.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "PrintAlignment.h"
#include "STLExtensions.h"
#include "lookup/JumpAligner.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTab.h" // Ted's LookupTab
#include "lookup/LookupTabBuilder.h"
#include "math/Functions.h"
#include "paths/GetNexts.h"
#include "paths/KmerPath.h"
#include "util/ReadTracker.h"
#include "util/RunCommand.h"
#include "feudal/BinaryStream.h"


#define is_in(x,v) (find(v.begin(), v.end(), x) != v.end())

void gather_kmers(const jump_align& ja, const vecKmerPath& unipaths, const int K, vec<longlong>&kmers)
{
  int kmer_length = max((u_int)K, ja.trusted_length) - K + 1;
  unsigned int offset = ja.target_path[0].offset;
  for (unsigned int i = 0; i < ja.target_path.size(); ++i) {
    unsigned int u = ja.target_path[i].contig;
    KmerPath unipath = unipaths[u];
    unsigned int ulength = unipath.TotalLength();
    for (unsigned int j = offset;
	 kmer_length > 0  && j < ulength;
	 ++j, --kmer_length) {
      kmers.push_back(unipath.GetKmer(j));
    }
    offset = 0;
  }
}

void gather_bases(const jump_align& ja, vecbvec const& contigs, bvec& bv)
{
  const jump_align_path& jap = ja.target_path;
  jump_align_path_iterator japi(jap, contigs);
  for (u_int i = 0; i < ja.trusted_length; ++i, ++japi)
    bv.AppendBase(*japi);
}

/**
 * ErrorCorrectJump
 *
 * Do error correction of jumping reads by aligning them onto a set of
 * unibases (plus adjacency graph), and replacing the bases of the
 * reads with those of the underlying unibase(s). If FLIP is true, run
 * error correction by seeding (and extending) from the last K-mer in
 * the read, rather than from the first.
 *
 * Input files:
 *   <REF_HEAD>.unibases.k<K>
 *   <REF_HEAD>.unipaths.k<K>
 *   <REF_HEAD>.unipath_adjgraph.k<K>
 *   <QUERY_HEAD>.fastb
 *   <QUERY_HEAD>.qualb
 *   <QUERY_HEAD>.pairs
 *
 * Output files:
 *   <OUT_HEAD>.orig_id
 *   <OUT_HEAD>.fastb
 *   <OUT_HEAD>.paths.k<K>
 *   <OUT_HEAD>.pairs*
 *   <OUT_HEAD>.qltout     (only if SAVE_ALIGNS=True)
 *
 * REF_HEAD: base name for unibases and unipaths
 * QUERY_HEAD: base name for query reads
 * OUT_HEAD: base name for output reads
 * NUM_THREADS: how many threads to use
 * FLIP: if true, rc reads, error correct, and then rc back
 * SAVE_ALIGNS: if true, save the alignments as well
 */
int main( int argc, char **argv )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String( QUERY_HEAD );
  CommandArgument_String( OUT_HEAD );
  CommandArgument_Int_OrDefault( TARGET_READ_LEN, 0);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault( FLIP, True );
  CommandArgument_Bool_OrDefault( SAVE_ALIGNS, False );
  CommandArgument_Int_OrDefault( MAX_KMER_FREQ, 10000 );
  CommandArgument_Int_OrDefault( MIN_MATCH, 20 );
  CommandArgument_Int_OrDefault( MISMATCH_THRESHOLD, 3 );
  CommandArgument_Int_OrDefault( MISMATCH_NEIGHBORHOOD, 8 );
  CommandArgument_Int_OrDefault( MISMATCH_BACKOFF, 3 );
  CommandArgument_Int_OrDefault( SCORE_DELTA, 20 );
  CommandArgument_Int_OrDefault( SCORE_MAX, 100 );
  CommandArgument_Bool_OrDefault( TRACK_READS, True );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  CommandArgument_Bool_OrDefault( WRITE, True );
  CommandArgument_String_OrDefault( DEBUG_READS, "" );
  CommandArgument_Bool_OrDefault( CROSS_UNIPATH_BOUNDARIES, True );
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // We need to extend to at least K bases to have any kmer paths!
  TARGET_READ_LEN = max(TARGET_READ_LEN, K);

  // File names.
  String strK = ToString( K );

  String ref_bases = REF_HEAD + ".unibases.k" + strK;
  String ref_tab = REF_HEAD + ".unibases.k" + strK + ".lookup";
  String ref_unipaths = REF_HEAD + ".unipaths.k" + strK;
  String ref_adjgraph = REF_HEAD + ".unipath_adjgraph.k" + strK;

  String query_bases = QUERY_HEAD + ".fastb";
  String query_quals = QUERY_HEAD + ".qualb";
  String query_pairs = QUERY_HEAD + ".pairs";

  String out_source = OUT_HEAD + ".orig_id";
  String out_bases = OUT_HEAD + ".fastb";
  String out_kpaths = OUT_HEAD + ".paths.k" + strK;
  String out_qlt = OUT_HEAD + ".qltout";

  // Generate a LookupTab from the unibases and write it to file.
  // (Skip this if the file already exists.)
  if ( !IsRegularFile( ref_tab ) ) {
    cout << Date( ) << ": Creating LookupTab file" << endl;
    LookupTabBuilder lookup_tab_builder( 12 );
    lookup_tab_builder.addFastb( ref_bases.c_str( ) );
    lookup_tab_builder.write( ref_tab.c_str( ) );
  }

  // Load.
  cout << Date( ) << ": loading reference bases" << endl;
  vecbvec ubases( ref_bases );

  cout << Date( ) << ": loading reference unipaths" << endl;
  vecKmerPath unipaths( ref_unipaths );

  cout << Date( ) << ": loading unipath adjency graph" << endl;
  digraph ugraph; 
  BinaryReader::readFile( ref_adjgraph, &ugraph );

  cout << Date( ) << ": loading query bases" << endl;
  vecbvec query( query_bases );
  size_t n_reads = query.size( );
  if ( FLIP ) {
    cout << Date( ) << ": rc-ing query bases" << endl;
    for (size_t ii=0; ii<n_reads; ii++)
      query[ii].ReverseComplement( );
  }

  cout << Date( ) << ": loading query quals" << endl;
  vecqvec quals( query_quals );
  if ( FLIP ) {
    cout << Date( ) << ": rc-ing query quals" << endl;
    for (size_t ii=0; ii<n_reads; ii++)
      reverse( quals[ii].begin( ), quals[ii].end( ) );
  }
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs;
  size_t n_pairs = pairs.nPairs( );

  // Get next unibases.
  vec< vec<int> > nexts;
  GetNexts( K, ubases, nexts );

  // Set up filtering options for call to FirstLookup.
  JumpAlignFilter lookup_filter;
  lookup_filter.orientation = FLIP ?
    JumpAlignFilter::RC_ONLY :
    JumpAlignFilter::FW_ONLY;
  lookup_filter.min_size = 20;
  lookup_filter.max_kmer_freq = MAX_KMER_FREQ;
  lookup_filter.score_delta = SCORE_DELTA;
  lookup_filter.score_max = SCORE_MAX;
  lookup_filter.min_match = MIN_MATCH;
  lookup_filter.mismatch_threshhold = MISMATCH_THRESHOLD;
  lookup_filter.mismatch_neighborhood = MISMATCH_NEIGHBORHOOD;
  lookup_filter.mismatch_backoff = MISMATCH_BACKOFF;
  lookup_filter.cross_unipath_boundaries = CROSS_UNIPATH_BOUNDARIES;
  lookup_filter.target_length = TARGET_READ_LEN;
  lookup_filter.max_placements = 2; // only care about single vs multiple
  ParseIntSet(DEBUG_READS, lookup_filter.debug_reads, False);

  // Run FirstLookup.
  vec<jump_align> hits;
  {
    // Load the LookupTab.
    cout << Date( ) << ": Loading LookupTab from file " << ref_tab << endl;
    LookupTab lookup_tab( ref_tab.c_str( ) );

    cout << Date( ) << ": running JumpAligner" << endl;
    JumpAligner jump_aligner( lookup_filter, lookup_tab, ubases, ugraph, K);
    jump_aligner.getAllAlignments( query, quals, hits, NUM_THREADS );
  }

  //Sort( hits, order_lookalign_QueryIdFullTargetMutations( ) );

  // Save alignments.
  //if ( SAVE_ALIGNS ) {
  //  cout << Date( ) << ": saving " << hits.size( ) << " hits" << endl;
  //  WriteLookAligns( out_qlt, hits );
  //}

  // Generate placement map.
  int discarded = -2;
  int unplaced = -1;
  vec<int> placement( n_reads, unplaced );
  for (int ii=0; ii<hits.isize( ); ii++) {
    int qid = hits[ii].query_id;

    // Already discarded.
    if ( placement[qid] != discarded ) {

      // SANTEMP - for now discard multiply placed reads - this is bad!
      if (placement[qid] >= 0 ) {
        placement[qid] = discarded;
      }
      else {  // Ok, set placement.
        placement[qid] = ii;
      }
    }
  }

  size_t unplaced_count = placement.CountValue(unplaced);
  size_t discarded_count = placement.CountValue(discarded);

  cout << "Alignments: " << hits.size() << endl << endl;
  cout << "Read count: " << n_reads << endl;
  cout << "Aligned uniquely: " << n_reads - unplaced_count - discarded_count << endl;
  cout << "Aligned non-uniquely : " <<  discarded_count << endl;
  cout << "Unaligned: " << unplaced_count << endl;
  cout << "Reads discarded : " << unplaced_count + discarded_count << endl;

  // Output.
  IncrementalWriter<bvec> b_out( out_bases.c_str( ) );
  IncrementalWriter<KmerPath> k_out( out_kpaths.c_str( ) );
  PairsManager pairs_out( n_reads );
  ofstream orig_out( out_source.c_str( ) );
  longlong n_good_pairs = 0, n_extended = 0;
  ReadTracker rt;
  unsigned int rt_source = 0;

  if (TRACK_READS) {
    rt_source = rt.AddSource(QUERY_HEAD);
  }

  // Parse all pairs.
  int dotter = 1e6;

  cout << Date( ) << ": parsing "
       << n_pairs << " pairs (.="
       << dotter << " pairs):"
       << endl;


  for (size_t pair_id = 0; pair_id < n_pairs; pair_id++) {
    if ( pair_id % dotter == 0 )
      Dot( cout, pair_id / dotter );

    // Only accept pairs for which both ends align uniquely.
    longlong id1 = pairs.ID1( pair_id );
    longlong id2 = pairs.ID2( pair_id );

    Bool debug = is_in(id1, lookup_filter.debug_reads) || is_in(id2, lookup_filter.debug_reads);

    if ( placement[id1] >= 0 && placement[id2] >= 0 ) {
      const jump_align &hit1 = hits[ placement[id1] ];
      const jump_align &hit2 = hits[ placement[id2] ];

//       if (hit1.extended) {
// 	PRINT4(id1, n_good_pairs*2, hit1.trusted_length, hit1.score);
//       }
//       if (hit2.extended) {
// 	PRINT4(id2, n_good_pairs*2+1, hit2.trusted_length, hit2.score);
//       }



      // sanity
      ForceAssertEq( hit1.query_id, id1 );
      ForceAssertEq( hit2.query_id, id2 );

      if (hit1.extended) ++n_extended;
      if (hit2.extended) ++n_extended;

      orig_out << hit1.query_id << "\n" << hit2.query_id << "\n";
      if (TRACK_READS) {
	rt.AddRead(rt_source, hit1.query_id);
	rt.AddRead(rt_source, hit2.query_id);
      }

      basevector b1, b2;
      gather_bases(hit1, ubases, b1);
      gather_bases(hit2, ubases, b2);
      b_out.add(b1);
      b_out.add(b2);

      vec<longlong> kmers1, kmers2;
      gather_kmers(hit1, unipaths, K, kmers1);
      gather_kmers(hit2, unipaths, K, kmers2);
      //ForceAssertEq(kmers1.size(), hit1.query_length-K+1);
      //ForceAssertEq(kmers2.size(), hit2.query_length-K+1);

      KmerPath kp1(kmers1), kp2(kmers2);
      k_out.add(kp1);
      k_out.add(kp2);

      if (VERBOSE) PRINT2( pair_id, pairs_out.nPairs( ) );
      pairs_out.addPair( 2*n_good_pairs, 2*n_good_pairs + 1,
			 pairs.sep( pair_id ), pairs.sd( pair_id ),
			 pairs.libraryName( pair_id ), True );
      n_good_pairs++;
    }
  }

  orig_out.close( );
  cout << endl;

  if (TRACK_READS) rt.Dump(OUT_HEAD);

  PRINT2( n_good_pairs, n_extended);
  if ( !WRITE ) exit(0);

  // Save.
  cout << Date( ) << ": merging output fastb and KmerPaths" << endl;
  b_out.close( );
  k_out.close( );

  cout << Date( ) << ": saving pairing info" << endl;
  // TEMP: For now, save the read pairs in the old format ( vec<read_pairing> )
  // rather than as a PairsManager.
  //vec<read_pairing> pairs_out_2 = pairs_out.convert_to_read_pairings( );
  //WritePairs( "", pairs_out_2, 2 * n_good_pairs, False, OUT_HEAD );
  pairs_out.Write( OUT_HEAD + ".pairs" );

  // Done.
  cout << Date( ) << ": Done with ErrorCorrectJump!" << endl;
  return 0;
}


