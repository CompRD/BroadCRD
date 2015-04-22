// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "JustifyPads.h"
#include "Loader.h"
#include "Nobbit.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "ReadLocation.h"
#include "ScoreAlignment.h"
#include "TaskTimer.h"
#include "pairwise_aligners/MakeAlignsStatic.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "tiled/AlignsOnConsensus.h"
#include "tiled/TAlign.h"
#include "random/Random.h"

/**
 * AlignsOnConsensus
 */
void AlignsOnConsensus( vec<t_align> &the_taligns,
			int contig_id,
			const vecbasevector &consensus_bases,
			const vecqualvector &consensus_quals,
			const vecbasevector &read_bases,
			const vecqualvector &read_quals,
			const vec<int> &first_locs,
			const vec<read_location> &locs,
			bool trust_locs,
			int max_discrep,
			int log_level,
			bool fast,
			int *max_gap,
			ostream *out,
			bool checkPerformance,
			int numRandomReads )
{
  ofstream devnull( "dev/null" );
  ostream &log = out ? *out : devnull;

  // Clear output.
  the_taligns.clear( );

  // Read ids.
  int first_loc = first_locs[contig_id];
  
  vec<int> ids;
  map<int, int> id_to_cgpos;
  for (int ii=first_loc; ii<(int)locs.size( ); ii++) {
    if ( locs[ii].Contig( ) != contig_id )
      break;
    ids.push_back( locs[ii].ReadId( ) );
    id_to_cgpos[locs[ii].ReadId( )] = ii - first_loc;
  }

  if ( checkPerformance )
  {
    vec<int> dummy_ids;
    for ( int i=0; i<numRandomReads; ++i )
    {
      int n = randomx( ) % ids.size();
      dummy_ids.push_back( locs[first_loc+n].ReadId() );
    }       
    sort(dummy_ids.begin(), dummy_ids.end());  
    dummy_ids.erase( unique( dummy_ids.begin( ), dummy_ids.end( ) ), dummy_ids.end( ) );

    ids = dummy_ids;
  }


  // Prepare vecbasevector to be used by MakeAligns.
  const basevector &master_bases = consensus_bases[contig_id];
  const qualvector &master_quals = consensus_quals[contig_id];
  vec<bvec const*> all_bases_b;
  all_bases_b.reserve(ids.size()+1);
  all_bases_b.push_back(&master_bases);
  for (size_t ii=0; ii<ids.size( ); ii++)
    all_bases_b.push_back( &read_bases[ids[ii]] );

  // Some logging.
  longlong contig_length = master_bases.size( );
  if ( contig_length < 1 ) {
    log << "Warning: contig_" << contig_id << " has zero length" << endl;
    return;
  }
  longlong tot_read_length = 0;
  for (int ii=0; ii<(int)ids.size( ); ii++)
    tot_read_length += read_bases[ids[ii]].size( );
  if ( tot_read_length < 1 ) {
    log << "Warning: contig_" << contig_id << " has no reads" << endl;
    return;
  }
  float ratio = SafeQuotient( tot_read_length, contig_length );
  String str_ratio = ToString( ratio, 2 ) + "%";

  log << "Aligning " << ids.size( )
      << " reads to contig_" << contig_id
      << " (len=" << contig_length
      << " cov=" << str_ratio
      << ")" << flush;
  
  // MakeAligns commmon part.
  static vec<align_plus> answer;
  static vec<int> answer_counts;
  int PASSES = 10;
  int ACTUAL_PASSES = 10;
  int MAXCLIQ = fast ? 100 : 1000;
  int MAX_ALIGNS = fast ? 100 : 10000;
  int MIN_MUTMER = 0;
  int MAKEALIGNS_MAXERRS = 200;
  int MAKEALIGNS_END_STRETCH = 128;
  to_compare what_to_compare( FIRST_VS_SECOND, 1 );
  Bool AVOID_PROMISCUOUS_KMERS = False;

  // SW gap method specific.
  int MIN_MAX_MUTMER = min<int>( 200, contig_length/2 );
  int MAX_OFFSET_DIFF = 1000;
  int MAX_GAP = 1000;
  int MIN_PROG_LENGTH = 0;
  float MIN_PROG_RATIO = 0.6;
  float MIN_OVERLAP_FRAC_LENGTH = 0.2;
  Bool AFFINE_PENALTY = True;
  
  makealigns_sw_gap_method *swAffineMethodPtr = new makealigns_sw_gap_method;
  
  swAffineMethodPtr->SetMaxErrs( MAKEALIGNS_MAXERRS );
  swAffineMethodPtr->SetEndStretch( MAKEALIGNS_END_STRETCH );
  swAffineMethodPtr->SetMinMaxMutmerLength( MIN_MAX_MUTMER );
  swAffineMethodPtr->SetMaxMutmerOffsetDiff( MAX_OFFSET_DIFF );
  swAffineMethodPtr->SetMaxGap( MAX_GAP );
  swAffineMethodPtr->SetMinProgressionLength( MIN_PROG_LENGTH );
  swAffineMethodPtr->SetMinProgressionRatio( MIN_PROG_RATIO );
  swAffineMethodPtr->SetMinOverlapFraction( MIN_OVERLAP_FRAC_LENGTH );
  swAffineMethodPtr->SetAffinePenalties( AFFINE_PENALTY );

  // Orig method specific.
  makealigns_orig_method *origMethodPtr = new makealigns_orig_method;

  int MAKEALIGNS_MAXERRS_LOCAL = 200;
  int MAKEALIGNS_MAXERRS_LOCAL_DONE = 0;
  int MAKEALIGNS_NSTRETCH = 10;
  int MAKEALIGNS_STRETCH = 20;
  int MAKEALIGNS_CL = 20;
  int MAKEALIGNS_CTOR_CALLS = 1000000;

  origMethodPtr->SetMaxBadness( MAKEALIGNS_MAXERRS );
  origMethodPtr->SetMaxErrs( MAKEALIGNS_MAXERRS );
  origMethodPtr->SetLocalMaxErrs( MAKEALIGNS_MAXERRS_LOCAL );
  origMethodPtr->SetLocalMaxErrsDone( MAKEALIGNS_MAXERRS_LOCAL_DONE );
  origMethodPtr->SetStretch( MAKEALIGNS_STRETCH );
  origMethodPtr->SetEndStretch( MAKEALIGNS_END_STRETCH );
  origMethodPtr->SetNStretch( MAKEALIGNS_NSTRETCH );
  origMethodPtr->SetCl( MAKEALIGNS_CL );
  origMethodPtr->SetMaxAlignCtorCalls( MAKEALIGNS_CTOR_CALLS );

  // Storage for aligns found by MakeAligns (with their scores).
  vec< pair<t_align, float> > align_score;
  align_score.reserve( ids.size( ) );

  // Two runs loop. The first run is fast but may leave out some good
  //  alignments (although most of the reads do get aligned). The second
  //  run is slower (not a problem since very few reads, if any, get to
  //  this stage), but is more aggressive, and aligns most that is left.
  for (int run_level=0; run_level<2; run_level++) {
    
    // All was alinged in run_level 0.
    if ( align_score.size( ) == ids.size( ) )
      break;
  
    // Remove from all_bases_b the reads that have been aligned successfully.
    bvec empty;
    if ( run_level == 1 ) {
      vec<int> placed_ids;
      placed_ids.reserve( align_score.size( ) );
      for (int ii=0; ii<(int)align_score.size( ); ii++)
	placed_ids.push_back( align_score[ii].first.id_ );
      sort( placed_ids.begin( ), placed_ids.end( ) );
      for (int ii=1; ii<(int)all_bases_b.size( ); ii++) {
	int id = ids[ii-1];
	if ( binary_search( placed_ids.begin( ), placed_ids.end( ), id ) ) {
	  all_bases_b[ii] = &empty;
	}
      }
    }

    // Align.
    TaskTimer align_time;
    align_time.Start( );

    if ( run_level == 0 ) {
      makealigns_method *methodPtr = swAffineMethodPtr;
      MakeAlignsStatic<2, 48, 50>( answer,
				   answer_counts,
				   PASSES,
				   ACTUAL_PASSES,
				   all_bases_b,
				   what_to_compare, 
				   methodPtr,
				   MAXCLIQ, 
				   MAX_ALIGNS,
				   MIN_MUTMER,
				   AVOID_PROMISCUOUS_KMERS );
    }
    else {
      makealigns_method *methodPtr = origMethodPtr;
      MakeAlignsStatic<2, 24, 50>( answer,
				   answer_counts,
				   PASSES,
				   ACTUAL_PASSES,
				   all_bases_b,
				   what_to_compare, 
				   methodPtr,
				   MAXCLIQ, 
				   MAX_ALIGNS,
				   MIN_MUTMER,
				   AVOID_PROMISCUOUS_KMERS );
    }
    
    align_time.Stop( );

    // Get rid of improper aligns, and save alignment scores.
    vec<Bool> keepers( answer.size( ), False );
    vec<float> scores( answer.size( ), -1 );

    int answer_ptr = 0;
    for ( int ii=0; ii<(int)answer_counts.size( ); ii++) {
      for (int jj=0; jj<(int)answer_counts[ii]; jj++) {
	align_plus &ap = answer[answer_ptr];
	int id = ids[ap.id2 - 1];
	const basevector &b1 = master_bases;
	const qualvector &q1 = master_quals;
	const basevector &b2 = read_bases[id];
	const qualvector &q2 = read_quals[id];
	
	// Move pointer.
	answer_ptr++;
	
	// Improper align.
	if ( !Proper( ap.a, b1.size( ), b2.size( ) ) )
	  continue;

	// Too large gaps.
	if ( max_gap ) {
	  int al_large_gap = 0;
	  const align &al = ap.a;
	  for (int tt=0; tt<al.Nblocks( ); tt++)
	    al_large_gap = Max( al_large_gap, Abs( al.Gaps( tt ) ) );
	  if ( al_large_gap > *max_gap )
	    continue;
	}
	
	// Wrong orientation.
	int loc_id = -1;
	if ( trust_locs ) {
	  map<int,int>::iterator iter;
	  iter = id_to_cgpos.find( id );
	  ForceAssert( iter != id_to_cgpos.end( ) );
	  loc_id = first_loc + iter->second;
	  int rc2 = ( ReverseOr == locs[loc_id].OrientationOnContig( ) );
	  if ( ap.RC != rc2 )
	    continue;
	}
	
	// This align would place read too far away from loc_expected site.
	if ( trust_locs ) {
	  int loc_start = locs[loc_id].StartOnContig( );
	  int al_start = ap.a.pos1( );
	  if ( Abs( loc_start - al_start ) > max_discrep ) 
	    continue;
	}	

	// Not enough of the read aligns consensus.
	if ( !trust_locs ) {
	  int read_begin = ap.a.pos1( );
	  int read_end = ap.a.Pos1( );
	  int align_len = read_end - read_begin;
	  int read_len = read_bases[id].size( );
	  if ( align_len < read_len / 2 )
	    continue;
	}
	
	// A keeper.
	float score = ScoreAlignment( ap.RC, ap.a, b1, q1, b2, q2 );
	keepers[answer_ptr-1] = True;
	scores[answer_ptr-1] = score;
      }
    }
    
    // For each read, select best alignment.
    vec<Bool> best( answer.size( ), False );

    // Winner align is the one placing read closer to initial locs estimate.
    if ( trust_locs ) {
      vec< pair< int, pair<int, int> > > chubber;
      chubber.reserve( answer.size( ) );
      for (int ii=0; ii<(int)answer.size( ); ii++) {
	if ( !keepers[ii] )
	  continue;
	int chub1 = answer[ii].id2 - 1;
	int read_id = ids[answer[ii].id2 - 1];
	map<int,int>::iterator iter = id_to_cgpos.find( read_id );
	ForceAssert( iter != id_to_cgpos.end( ) );
	int loc_id = first_loc + iter->second;
	int loc_start = locs[loc_id].StartOnContig( );
	int distance = Abs( answer[ii].a.pos1( ) - loc_start );
	pair<int,int> chub2 = make_pair( distance, ii );
	chubber.push_back( make_pair( chub1, chub2 ) );
      }
      sort( chubber.begin( ), chubber.end( ) );
      
      int last_best = -1;
      for (int ii=0; ii<(int)chubber.size( ); ii++) {
	int answer_id = chubber[ii].second.second;
	int pos_id = chubber[ii].first;
	if ( !keepers[answer_id] )
	  continue;      
	if ( pos_id != last_best ) {
	  best[answer_id] = True;
	  last_best = pos_id;
	}
      }      
    }
    
    // Winner align is lowest score.
    else {
      vec< pair< int, pair<float, int> > > chubber;
      chubber.reserve( answer.size( ) );
      for (int ii=0; ii<(int)answer.size( ); ii++) {
	if ( !keepers[ii] )
	  continue;
	int chub1 = answer[ii].id2 - 1;
	pair<float,int> chub2 = make_pair( scores[ii], ii );
	chubber.push_back( make_pair( chub1, chub2 ) );
      }
      sort( chubber.begin( ), chubber.end( ) );
      
      int last_best = -1;
      for (int ii=0; ii<(int)chubber.size( ); ii++) {
	int answer_id = chubber[ii].second.second;
	int pos_id = chubber[ii].first;
	if ( !keepers[answer_id] )
	  continue;      
	if ( pos_id != last_best ) {
	  best[answer_id] = True;
	  last_best = pos_id;
	}
      }
    }
    
    // Store alignments.
    int n_found = 0;
    for (int ii=0; ii<(int)answer.size( ); ii++) {
      if ( !best[ii] ) continue;

      // The alignment.
      int read_id = ids[answer[ii].id2 - 1];
      float the_score = scores[ii];
      t_align the_tal( read_id, answer[ii].RC, answer[ii].a );
      the_tal.al_.Flip( );

#ifdef SKIP // Leave it out for now, this slows down the code too much.
      // Run affine Smith-Waterman (clean up the align).
      bool isRC = answer[ii].RC;
      float score;
      basevector rc_read;
      if ( isRC ) {
	rc_read = read_bases[read_id];
	rc_read.ReverseComplement( );
      }
      const basevector &rb = isRC ? rc_read : read_bases[read_id];

      const basevector &cborig = consensus_bases[contig_id];
      int start = the_tal.al_.pos2( );
      int end = the_tal.al_.Pos2( );
      int len = end - start;
      basevector cb;
      cb.SetToSubOf( cborig, start, len );

      alignment theAlignment;
      theAlignment.Set( the_tal.al_ );
      SmithWatAffine( rb, cb, theAlignment );
      int pos1;
      int pos2;
      int errors;
      avector<int> gaps;
      avector<int> lens;
      theAlignment.Unpack( pos1, pos2, errors, gaps, lens );

      // Correct align only if start was not changed.
      if ( pos2 == 0 )
	the_tal.al_.Set( pos1, start, gaps, lens );
#endif
      
      // Add align.
      align_score.push_back( make_pair( the_tal, the_score ) );
      n_found++;      
    }
    
    // Log percent aligned and run time.
    float ratio = 100.0 * SafeQuotient( n_found, (int)ids.size( ) );
    String percent_found = ToString( ratio, 2 ) + "%";
    String run_time = ToString( align_time.Elapsed( ), 2 ) + "s";

    log << " [run_" << run_level
 	<< ": " << n_found
 	<< " aligned (" << percent_found
	<< ") in " << run_time
 	<< "]" << flush;
  }
  
  // Fill output.
  the_taligns.reserve( align_score.size( ) );
  for (int ii=0; ii<(int)align_score.size( ); ii++)
    the_taligns.push_back( align_score[ii].first );

  // Final base logging.
  if ( the_taligns.size( ) < ids.size( ) ) {
    int miss = ids.size( ) - the_taligns.size( );
    log << " missing " << miss << " reads";
  }
  
  log << endl;

  // Verbose logging now.
  if ( log_level > 1 ) {
    vec<String> printables( ids.size( ), "" );
    vec<int> alscore_pos( ids.size( ), -1 );
    
    // printables for each read (either "" or a full line);
    String str_n_reads_tot = ToString( (int)ids.size( ) - 1 );
    map<int, int>::const_iterator iter;
    for (int ii=0; ii<(int)align_score.size( ); ii++) {
      int read_id = align_score[ii].first.id_;
      float al_score = align_score[ii].second;
      iter = id_to_cgpos.find( read_id );
      ForceAssert( iter != id_to_cgpos.end( ) );
      int cgpos = iter->second;
      String str_orient = align_score[ii].first.isRC_ ? "_rc" : "";

      alscore_pos[cgpos] = ii;

      printables[cgpos]
	= " " + ToString( cgpos )
	+ "." + str_n_reads_tot
	+ ": c" + ToString( contig_id )
	+ " vs r" + ToString( read_id ) + str_orient
	+ " (len=" + ToString( read_bases[read_id].size( ) )
	+ ", score=" + ToString( al_score, 2 )
	+ ")";
    }

    // replace "" lines in printables with meaningful entries.
    for (int ii=0; ii<(int)printables.size( ); ii++ ) {
      if ( printables[ii] != "" )
	continue;

      int read_id = ids[ii];

      int q20 = 0;
      for (int jj=0; jj<(int)read_quals[read_id].size( ); jj++)
	if ( int( read_quals[read_id][jj] ) >= 20 )
	  q20++;
      
      printables[ii]
	= " " + ToString( ii )
	+ "." + str_n_reads_tot
	+ ": c" + ToString( contig_id )
	+ " vs r" + ToString( read_id )
	+ " align_not_found (len=" + ToString( read_bases[read_id].size( ) )
	+ " q20=" + ToString( q20 )
	+ ")";
    }
    
    // Print alignments with high score.
    for (int ii=0; ii<(int)printables.size( ); ii++) {
      log << printables[ii] << "\n";
      
      // Align not found.
      if ( alscore_pos[ii] < 0 )
	continue;
      
      // Do not print alignments.
      if ( log_level < 3 )
	continue;
      
      // Good align found.
      if ( log_level == 3 && align_score[alscore_pos[ii]].second < 100.0 )
	continue;

      // Print high-score align.
      align flip_al = align_score[alscore_pos[ii]].first.al_;
      flip_al.Flip( );
      bool alRC = align_score[alscore_pos[ii]].first.isRC_;
      int read_id = align_score[alscore_pos[ii]].first.id_;
      const basevector &b1 = consensus_bases[contig_id];
      const qualvector &q1 = consensus_quals[contig_id];
      const basevector &b2 = read_bases[read_id];
      const qualvector &q2 = read_quals[read_id];
      
      PrintVisualAlignment( alRC, True, log, b1, b2, flip_al, q1, q2 );
    }
  }
  
  // Clear memory.
  delete( swAffineMethodPtr );
  delete( origMethodPtr );

}

/**
 * AlignsOnConsensus
 */
void AlignsOnConsensus( vec<t_align> &the_taligns,
			int contig_id,
			const vecbasevector &cbases,
			const vecbasevector &rbases,
			const lhandler &locs,
			int max_discrep,
			float max_error,
			ostream *log,
			ostream *fastalog)
{
  // A realigner.
  CRealign aligner( cbases, rbases, locs, log, fastalog );
  aligner.SetVerbose( false );
  aligner.SetContig( contig_id );
  aligner.SetMaxDiscrep( max_discrep );
  aligner.SetMaxError( max_error );
  
  // Empty contig.
  int floc = locs.FirstLoc( contig_id );
  if ( floc < 0 )
    return;

  // Loop over locs in contig.
  t_align tal;
  for (int ii=floc; ii<locs.Size( ); ii++) {
    if ( locs[ii].Contig( ) != contig_id )
      break;
    if ( ! aligner.PlaceLoc( ii, tal ) ) {
      if ( log ) *log << "ERROR - read not placed.\n\n";
      continue;
    }
    the_taligns.push_back( tal );
  }
  if ( log ) *log << endl;

}

