///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Fastavector.h"
#include "Histogram.h"
#include "Qualvector.h"
#include "PrintAlignment.h" 
#include "btl/ClusterOffsets.h"
#include "btl/SanitizeSmithWatBandedA2.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/KmerAligner.h"
#include "pairwise_aligners/PerfectAlignment.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "util/RunCommand.h"
  
/**
 * SelectReadsVNTR
 *
 * Select PacBio reads spanning a VNTR region. The reference is used
 * to infer the unique flanking regions.
 *
 * INPUT:
 *   <IN_DIR>/reads.{fastb,qualb,ids}
 *   <REF_HEAD>.fasta
 *
 * OUTPUT:
 *   <OUT_DIR>.anchors.fastb         bases of anchoring segments
 *   <OUT_DIR>.anchors_to_reads.qlt  aligns of anchors on reads
 *   <OUT_DIR>.reads_to_ref.qlt      aligns of reads on reference
 *   <OUT_DIR>.spanning.qlt          matching aligns of spanning reads' anchors
 *   <OUT_DIR>.spanning.names        names of spanning reads
 *   <OUT_DIR>.spanning.visual       visual aligns of reads on reference
 *   ...                             a few other log/info files
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( IN_DIR );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String( REF_HEAD );

  // Build left and right anchors out of ref.
  CommandArgument_Int_OrDefault( LEFT_ANCHOR, 700 );
  CommandArgument_Int_OrDefault( RIGHT_ANCHOR, 700 );
  CommandArgument_Int_OrDefault( MIN_VNTR_LEN, 2000 );
  CommandArgument_Int_OrDefault( MAX_VNTR_LEN, 3000 );

  // KmerAligner and Smith-Waterman args.
  CommandArgument_Int_OrDefault( STEP, 1 );
  CommandArgument_Int_OrDefault( MIN_SEEDS, 40 );
  CommandArgument_Int_OrDefault( BAND, 100 );
  CommandArgument_Int_OrDefault( MIN_ALIGN_LEN, 100 );
  CommandArgument_Double_OrDefault( MAX_ER, 0.35 );

  // Specify a bin size for the eps histogram plot.
  CommandArgument_Int_OrDefault( BIN_SIZE, 10 );
  
  // Align reads to reference, and print aligns.
  CommandArgument_Bool_OrDefault( ALIGN_READS_TO_REF, True );
  CommandArgument_Bool_OrDefault( PRINT_READS_TO_REF, False );
  
  // Do not use cached data.
  CommandArgument_Bool_OrDefault( FORCE, True );
  
  EndCommandArguments;
  
  // Dir and file names.
  String ref_fasta_file = REF_HEAD + ".fasta";

  String reads_fastb_file = IN_DIR + "/reads.fastb";
  String reads_qualb_file = IN_DIR + "/reads.qualb";
  String reads_ids_file = IN_DIR + "/reads.ids";

  String log_file = OUT_DIR + "/SelectReadsVNTR.log";
  String ref_fastb_file = OUT_DIR + "/reference.fastb";
  String anchors_fastb_file = OUT_DIR + "/anchors.fastb";
  String anchors_to_reads_qlt_file = OUT_DIR + "/anchors_to_reads.qlt";
  String reads_to_ref_qlt_file = OUT_DIR + "/reads_to_ref.qlt";
  String spanning_qlt_file = OUT_DIR + "/spanning.qlt";
  String spanning_names_file = OUT_DIR + "/spanning.names";
  String spanning_visual_file = OUT_DIR + "/spanning.visual";
  String histogram_file = OUT_DIR + "/histogram.eps";
  
  vec<String> needed;
  needed.push_back( ref_fasta_file );
  needed.push_back( reads_fastb_file );
  needed.push_back( reads_qualb_file );
  needed.push_back( reads_ids_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  Mkpath( OUT_DIR );

  ofstream log( log_file.c_str( ) );
  cout << "Sending log to " << log_file << "\n" << endl;
  PrintCommandPretty( log );
  
  // Load reference, create left and right anchors.
  vecbvec ref_bases;
  vecbvec anchors_bases;
  int right_anchor_start = 0;
  if ( FORCE || ! IsRegularFile( ref_fastb_file ) ) {
    log << Date( ) << ": building ref fastb" << endl;
    vec<fastavector> ref_fasta;
    LoadFromFastaFile( ref_fasta_file, ref_fasta );
    if ( ref_fasta.size( ) != 1 ) {
      log << "Fatal error: ref fasta must have exactly one entry.\n" << endl;
      return 1;
    }
    int ref_len = ref_fasta[0].size( );
    if ( ref_len < MIN_VNTR_LEN ) {
      log << "Fatal error: ref is too short.\n" << endl;
      return 1;
    }
    ref_bases.push_back( ref_fasta[0].ToBasevector( ) );
    right_anchor_start = ref_len - RIGHT_ANCHOR;
    bvec left_chunk( ref_bases[0], 0, LEFT_ANCHOR );
    bvec right_chunk( ref_bases[0], right_anchor_start, RIGHT_ANCHOR );
    anchors_bases.push_back( left_chunk );
    anchors_bases.push_back( right_chunk );
    ref_bases.WriteAll( ref_fastb_file );
    anchors_bases.WriteAll( anchors_fastb_file );
  }
  else {
    log << Date( ) << ": loading reference" << endl;
    ref_bases.ReadAll( ref_fastb_file );
    ForceAssertEq( (int)ref_bases.size( ), 1 );
    right_anchor_start = (int)ref_bases[0].size( ) - RIGHT_ANCHOR;
    anchors_bases.ReadAll( anchors_fastb_file );
  }
  
  vecbvec anchors_bases_rc;
  anchors_bases_rc.reserve( anchors_bases.size( ) );
  for (size_t ii=0; ii<anchors_bases.size( ); ii++) {
    bvec rcb = anchors_bases[ii];
    rcb.ReverseComplement( );
    anchors_bases_rc.push_back( rcb );
  }
  
  vecqvec anchors_quals;
  anchors_quals.reserve( anchors_bases.size( ) );
  for (size_t ii=0; ii<anchors_bases.size( ); ii++)
    anchors_quals.push_back( qvec( anchors_bases[ii].size( ), 50 ) );

  // Load reads.
  log << Date( ) << ": loading reads... " << flush;
  vecbvec reads_bases( reads_fastb_file );
  vecqvec reads_quals( reads_qualb_file );
  vecString reads_ids( reads_ids_file );
  const int n_reads = reads_bases.size( );
  log << ToStringAddCommas( n_reads ) << " found" << endl;

  // Align anchors to reads.
  if ( FORCE || ! IsRegularFile( anchors_to_reads_qlt_file ) ) {

    // Instantiate a KmerAligner.
    Assembly empty_assembly;
    KmerAligner<12> aligner12( empty_assembly );
    aligner12.SetBases( anchors_bases );
    aligner12.SetKmerStep( STEP );
    
    // Output stream for aligns of anchors to reads.
    ofstream a2rout( anchors_to_reads_qlt_file.c_str( ) );
    
    // Sometimes SmithWatBandedA2 returns invalid aligns.
    int n_invalid_als = 0;

    // Loop over all reads.
    int dotter = 1000;
    log << Date( ) << ": aligning anchors (. = " << dotter << " reads)\n";
    for (int rid=0; rid<n_reads; rid++) {
      if ( rid % dotter == 0 ) Dot( log, rid / dotter );
      
      // First align fw read, then tc read.
      for (int orient=0; orient<2; orient++) {
	
	// The read (called b2, as we'll align it to b1, aka the anchors).
	bvec rc_readb;
	qvec rc_readq;
	if ( orient == 1 ) {
	  rc_readb = reads_bases[rid];
	  rc_readq = reads_quals[rid];
	  rc_readb.ReverseComplement( );
	  rc_readq.ReverseMe( );
	}
	const bvec &b2 = ( orient == 1 ) ? rc_readb : reads_bases[rid];
	const qvec &q2 = ( orient == 1 ) ? rc_readq : reads_quals[rid];
	
	// Collect seeds.
	vec<int> ids;
	vec<int> starts;
	aligner12.FindPossibleAlignments( b2, starts, ids );
	
	// Two passes: align to left and right anchor.
	for (int pass=0; pass<2; pass++) {
	  const bvec &b1 = anchors_bases[pass];
	  const qvec &q1 = anchors_quals[pass];
	  const bvec &b1rc = anchors_bases_rc[pass];
	  const int b1len = b1.size( );
	  const int b2len = b2.size( );
	  
	  // Cluster offsets.
	  vec<int> offsets;
	  for (int ii=0; ii<ids.isize( ); ii++)
	    if ( ids[ii] == pass )
	      offsets.push_back( starts[ii] );
	  sort( offsets.begin( ), offsets.end( ) );
	  
	  vec<int> clusters;
	  ClusterOffsets( MIN_SEEDS, offsets, clusters );
	  
	  // Align anchor to read (1=anchor, 2=read).
	  for (int ii=0; ii<clusters.isize( ); ii++) {
	    int off = clusters[ii];
	    align al;
	    int err = 0;
	    SmithWatBandedA2<uint>( b1, b2, off, BAND, al, err, 0, 1, 1, 1 );
	    
	    // Sanitize align (walk around known bug in SmitWatBandedA2).
	    if ( ! SanitizeSmithWatBandedA2( al ) ) {
	      n_invalid_als++;
	      continue;
	    }

	    // Short align.
	    int al_len = al.Pos1( ) - al.pos1( );
	    if ( al_len < MIN_ALIGN_LEN ) continue;
	    
	    // Too many errors.
	    float score = 0;
	    PerfectAlignment( b1, q1, b2, q2, al, score );
	    err = al.Errors( b1, b2 );
	    if ( (float)err / (float)al_len > MAX_ER ) continue;
	    
	    // Align found, dump align.
	    bool rc1 = ( orient == 1 );

	    if ( rc1 ) {
	      // NB: ResetFromAligns will not rc b1, even if rc1 = true.
	      al.ReverseThis( b1.size( ), b2.size( ) );
	      look_align hit( pass, rid, b1len, b2len, rc1, al, 0, 0, 0 );
	      hit.ResetFromAlign( al, b1rc, reads_bases[rid] );
	      hit.PrintParseable( a2rout, &b1, &( reads_bases[rid] ) );
	    }
	    else {
	      look_align hit( pass, rid, b1len, b2len, rc1, al, 0, 0, 0 );
	      hit.ResetFromAlign( al, b1, b2 );
	      hit.PrintParseable( a2rout, &b1, &b2 );
	    }
	    
	  } // loop over all clusters of shared seeding k-mers
	} // loop over left/righ anchor
      } // loop over fw/rc of read
    } // loop over all reads
    log << endl;
    
    // Report invalid aligns.
    if ( n_invalid_als > 0 ) {
      log << "\nWARNING: found and discarded " << n_invalid_als
	  << " invalid align(s).  This is due to a bug\n"
	  << "   in SmithWatBandedA2, and it is not a big problem (unless\n"
	  << "   this number is too big).\n"
	  << endl;
    }
    
    // Close stream.
    a2rout.close( );

  } // Align anchors to reads
  
  // Load and sort aligns of anchors to reads.
  log << Date( ) << ": loading aligns of anchors to reads" << endl;
  vec<look_align> ahits;
  LoadLookAligns( anchors_to_reads_qlt_file, ahits );
  
  order_lookalign_TargetBegin sorter;
  sort( ahits.begin( ), ahits.end( ), sorter );
  vec<int> first_ids( n_reads, -1 );
  for (int ii=ahits.isize( )-1; ii>=0; ii--)
    first_ids[ ahits[ii].target_id ] = ii;
  
  // Various streams (some are used only if ALIGN_READS_TO_REF is true).
  int dotter = 1000;
  log << Date( ) << ": looking for spanning reads (. = "
      << dotter << " reads)\n";
  ofstream span_qlt_out( spanning_qlt_file.c_str( ) );
  ofstream span_names_out( spanning_names_file.c_str( ) );
  ofstream span_visual_out( spanning_visual_file.c_str( ) );
  ofstream reads_to_ref_out( reads_to_ref_qlt_file.c_str( ) );

  // List of observed separations between anchors (for the histogram).
  vec<int> seps;

  // Find reads spanning the whole VNTR.
  for (int rid=0; rid<n_reads; rid++) {
    if ( rid % dotter == 0 ) Dot( log, rid / dotter );
      
    // Unmapped anchors.
    if ( first_ids[rid] < 0 ) continue;

    // We need to see the two anchors just once.
    int n_seen = 0;
    for (int ii=first_ids[rid]; ii<ahits.isize( ); ii++) {
      if ( ahits[ii].target_id != rid ) break;
      n_seen++;
    }
    if ( n_seen != 2 ) continue;

    // Check orientation of anchors.
    const look_align &first_la = ahits[ first_ids[rid] ];
    const look_align &second_la = ahits[ first_ids[rid] + 1 ];
    if ( first_la.rc1 != second_la.rc1 ) continue;

    // Both anchors must be present.
    if ( first_la.query_id == second_la.query_id ) continue;

    // Check relative position of anchors.
    if ( first_la.rc1 && first_la.query_id == 0 ) continue;
    if ( ( ! first_la.rc1 ) && first_la.query_id == 1 ) continue;

    // Check separation.
    int sep = second_la.pos2( ) - first_la.Pos2( );
    if ( sep < MIN_VNTR_LEN || sep > MAX_VNTR_LEN ) continue;

    // Candidate found, save info.
    bool is_rc = first_la.rc1;
    bool is_LR = first_la.query_id == 0;
    const align &al1 = first_la.a;
    const align &al2 = second_la.a;

    first_la.PrintParseable( span_qlt_out );
    second_la.PrintParseable( span_qlt_out );
    span_qlt_out << "\n";

    span_names_out << reads_ids[rid] << "\t"
		   << rid << "\t"
		   << ( is_LR ? "L" : "R" ) << ( is_rc ? "[-]" : "[+]" ) << " "
		   << "[" << al1.pos2( ) << ", " << al1.Pos2( ) << ")\t"
		   << ( is_LR ? "R" : "L" ) << ( is_rc ? "[-]" : "[+]" ) << " "
		   << "[" << al2.pos2( ) << ", " << al2.Pos2( ) << ")\t"
		   << "read_len: " << reads_bases[rid].size( ) << "\t"
		   << "VNTR_len: " << sep << "\n";
    
    // Add sep to list.
    seps.push_back( sep );

    // Continue loop, if ALIGN_READS_TO_REF is false.
    if ( ! ALIGN_READS_TO_REF ) continue;
    
    // Align reads (1) to reference (2), and print align.
    int offset_left = al1.pos2( );
    int offset_right = al2.pos2( ) - right_anchor_start;
    
    int offset = ( offset_left + offset_right ) / 2;
    int max_off = Max( offset_left, offset_right );
    int min_off = Min( offset_left, offset_right );
    int sw_band = BAND + ( max_off - min_off );

    bvec rc_read;
    qvec rc_qual;
    if ( is_rc ) {
      rc_read = reads_bases[rid];
      rc_qual = reads_quals[rid];
      rc_read.ReverseComplement( );
      rc_qual.ReverseMe( );
    }
    const bvec &b1 = is_rc ? rc_read : reads_bases[rid];
    const qvec &q1 = is_rc ? rc_qual : reads_quals[rid];
    const bvec &b2 = ref_bases[0];
    
    align al;
    int err = 0;
    SmithWatBandedA2<uint>( b1, b2, offset, sw_band, al, err, 0, 1, 1, 1 );
    
    // Sanitize align (walk around known bug in SmitWatBandedA2).
    if ( ! SanitizeSmithWatBandedA2( al ) ) {
      log << "warning: SanitizeSmithWatBandedA2 failed on read " << rid << endl;
      continue;      
    }
    
    // Save align.
    int tid = 0;
    int rlen = b1.size( );
    int tlen = b2.size( );

    // NB: ResetFromAligns will not rc reads_bases[rid], even if is_rc is true.
    look_align hit( rid, tid, rlen, tlen, is_rc, al, 0, 0, 0 );
    hit.ResetFromAlign( al, b1, b2 );
    hit.PrintParseable( reads_to_ref_out );

    // Continue loop, if PRINT_READS_TO_REF is false.
    if ( ! PRINT_READS_TO_REF ) continue;
    
    // Show aligns (in its own log file).
    span_visual_out << "Alignment between " << ( is_rc ? " rc of " : "" )
     		    << rid << " and reference\n"
		    << "Using   offset = " << offset
		    << "   band = " << sw_band << "\n";
    PrintVisualAlignment( True, span_visual_out, b1, b2, al, q1 );
  }
  log << endl;

  // Close various streams.
  span_visual_out.close( );
  span_qlt_out.close( );
  span_names_out.close( );
  reads_to_ref_out.close( );

  // Generate histogram of separations.
  log << Date( ) << ": saving histogram of estimated region length" << endl;
  histogram<int> histo;
  vec<int> bins( 1, MIN_VNTR_LEN );
  while ( bins.back( ) < MAX_VNTR_LEN )
    bins.push_back( bins.back( ) + BIN_SIZE );
  histo.AddBins( bins );
  histo.AddData( seps.begin( ), seps.end( ) );
  {
    using namespace ns_psplot;
    
    vec<ns_psplot::freetext> labels( 1 );
    labels[0] = freetext( "Estimated length of VNTR region", black, 16 );
    ofstream eps_out( histogram_file.c_str( ) );
    histo.PrintAsEps( eps_out, labels, 0 );
    eps_out.close( );
  }

  // Done.
  log << "\n" << Date( ) << ": SelectReadsVNTR done" << endl;
  log.close( );
  
}
