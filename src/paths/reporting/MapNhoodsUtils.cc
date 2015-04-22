///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "lookup/LookAlign.h"
#include "paths/reporting/MapNhoodsUtils.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"

/**
 * EstimatePlacement
 */
triple<int,int,int> EstimatePlacement( const HyperBasevector &hyper,
				       const look_align &la,
				       const vec<int> &alens )
{
  triple<int,int,int> result;

  longlong total_nhood_len = hyper.TotalEdgeLength( );
  result.first  = la.target_id;
  result.second = Max( (longlong)0, la.a.pos2( ) - total_nhood_len );
  result.third  = Min( (longlong)alens[la.target_id] , la.a.Pos2( ) + total_nhood_len );
  
  return result;
}

/**
 * RefinePlacement
 */
void RefinePlacement( triple<int,int,int> &placement,
		      vec<look_align> &aligns,
		      const int K,
		      const int nhood_id,
		      const String &assembly_head,
		      const String &outdir_head,
		      const vec<fastavector> &assembly_fasta,
		      const vec<HyperBasevector> &hypers,
		      const Bool FORCE )
{
  // Dir and file names (common part).
  String out_dir = outdir_head + "/" + ToString( nhood_id / 1000 );
  
  String assembly_lookup_file = assembly_head + ".assembly.lookup";

  String out_head = out_dir + "/" + ToString( nhood_id );
  String fastb_file = out_head + ".fastb";
  String aligns_file = out_head + ".qlt";
  String chunk_head = out_head + ".chunk";
  String chunk_range_file = chunk_head + ".range";
  String chunk_fastb_file = chunk_head + ".fastb";
  String chunk_lookup_file = chunk_head + ".lookup";
  
  Mkpath( out_dir );
  
  // Various and sundry.
  const int target_id = placement.first;
  const int guess_beg = placement.second;
  const int guess_end = placement.third;
  const int guess_len = guess_end - guess_beg;
  const HyperBasevector &hyper = hypers[nhood_id];
  
  // Generate nhood fastb and range files.
  if ( FORCE || ! IsRegularFile( fastb_file ) ) {
    vecbvec bases;
    hyper.GenerateVecbasevector( bases );
    bases.WriteAll( fastb_file );

    ofstream rout( chunk_range_file.c_str( ) );
    rout << target_id << ":" << guess_beg << "-" << guess_end << "\n";
    rout.close( );
  }
  
  // Lookup file.
  if ( ( FORCE || ! IsRegularFile( aligns_file ) ) && target_id > -1 ) {
    vec<int> select( 1, target_id );
    fastavector fsel;
    fsel.SetToSubOf( assembly_fasta[target_id], guess_beg, guess_len );
    vecbvec chunk;
    chunk.push_back( fsel.ToBasevector( ) );
    chunk.WriteAll( chunk_fastb_file );
    
    String theCommand
      = "MakeLookupTable LOOKUP_ONLY=True SOURCE=" + chunk_fastb_file
      + " OUT_HEAD=" + chunk_head;
    RunCommandWithLog( theCommand, "/dev/null" );
  }
  String lookup_file = target_id < 0 ? assembly_lookup_file : chunk_lookup_file;

  // Generate (or load) aligns.
  GetAlignsFast( K, fastb_file, lookup_file, aligns_file,
		 aligns, ! FORCE, out_dir );
  
  // May have to adjust aligns.
  if ( target_id > -1 ) {
    for (int ii=0; ii<aligns.isize( ); ii++) {
      aligns[ii].SetTargetId( target_id );
      aligns[ii].SetStartOnTarget( aligns[ii].StartOnTarget( ) + guess_beg );
    }
  }

  // Generate reftigs.
  digraph agraph;
  vec< pair<int,ho_interval> > reftigs;
  HyperToReftigsCore( K, hyper, aligns, reftigs, &agraph );
  
  // Tune up placement.
  int new_target_id = -1;
  int new_beg = -1;
  int new_end = -1;
  for (int ii=0; ii<reftigs.isize( ); ii++) {
    new_target_id = reftigs[ii].first;
    if ( target_id > -1 )
      ForceAssertEq( reftigs[ii].first, target_id );
    if ( new_beg < 0 || reftigs[ii].second.Start( ) < new_beg )
      new_beg = reftigs[ii].second.Start( );
    if ( new_end < 0 || reftigs[ii].second.Stop( ) > new_end )
      new_end = reftigs[ii].second.Stop( );
  }

  // Update placement.
  if ( new_beg < 0 ) placement = triple<int,int,int>( -1, -1, -1 );
  else placement = triple<int,int,int> ( new_target_id, new_beg, new_end );
  
}

/**
 * MapRange
 */
void MapRange( ostream &log,
	       const String RANGE,
	       const Bool SKIP_UNMAPPED,
	       const Bool FORCE,
	       const int K,
	       const int n_seeds,
	       const String &out_dir,
	       const String &assembly_head,
	       const vec<fastavector> &assembly_fasta,
	       const vec<int> &assembly_lengths, 
	       const vec<int> &seeds_ids,
	       const vec<int> &aligns_index,
	       const vec<look_align> &seeds_aligns,
	       const vec<HyperBasevector> &hypers )
{
  // Parse RANGE.
  int range_id = -1;
  int range_beg = -1;
  int range_end = -1;
  if ( RANGE.Contains( ":" ) && RANGE.Contains( "-" ) ) {
    String str_id = RANGE.Before( ":" );
    String str_beg = RANGE.After( ":" ).Before( "-" );
    String str_end = RANGE.After( "-" );
    if ( str_id.IsInt( ) && str_beg.IsInt( ) && str_end.IsInt( ) ) {
      range_id = str_id.Int( );
      range_beg = str_beg.Int( );
      range_end = str_end.Int( );
    }
  }
  if ( range_id < 0 ) {
    log << "Invalid RANGE (q to quit).\n" << endl;
    return;
  }
  
  // Initial selection of nhoods.
  log << "\nCandidate nhoods for RANGE = " << RANGE << "\n" << endl;
  int n_selected = 0;
  int n_unmapped = 0;
  vec<bool> selected( seeds_ids.size( ), false );
  for (int nhood_id=0; nhood_id<n_seeds; nhood_id++) {
    if ( aligns_index[nhood_id] < 0 ) {
      if ( ! SKIP_UNMAPPED ) {
	n_selected++;
	selected[nhood_id] = true;
	log << "nhood " << nhood_id
	    << "\tunmapped (or multiply mapped) seed\n";
      }
      n_unmapped++;
      continue;
    }
    
    triple<int,int,int> placement( -1, -1, -1 );
    const HyperBasevector &hyper = hypers[nhood_id];
    const look_align &al = seeds_aligns[ aligns_index[nhood_id] ];
    placement = EstimatePlacement( hyper, al, assembly_lengths );
    int target = placement.first;
    int beg = placement.second;
    int end = placement.third;
    if ( beg < 0 ) {
      n_selected++;
      selected[nhood_id] = true;
      log << "nhood " << nhood_id
	  << "\tu_" << seeds_ids[nhood_id]
	  << "\ton t_" << target
	  << " seed mapped, but EstimatePlacement failed\n";
      continue;
    }
    
    if ( target != range_id ) continue;
    if ( end <= range_beg || range_end <= beg ) continue;
    
    n_selected++;
    selected[nhood_id] = true;
    log << "nhood " << nhood_id
	<< "\tu_" << seeds_ids[nhood_id]
	<< "\ton t_" << target
	<< " [" << beg
	<< ", " << end
	<< ")_" << assembly_lengths[target]
	<< "\n";
  }

  // Report unmapped nhoods.
  if ( SKIP_UNMAPPED ) {
    double ratio = 100.0 * SafeQuotient( n_unmapped, seeds_ids.isize( ) );
    log << "\nWarning: skipping " << n_unmapped
	<< " unmapped nhoods (out of " << seeds_ids.size( ) 
	<< " nhoods - " << ToString( ratio, 1 )
	<< "% of the total)\n";
  }

  // Tuning up alignments.
  log << "\nAligning " << n_selected << " candidate nhoods\n" << endl;
  int running = 0;
  for (int nhood_id=0; nhood_id<n_seeds; nhood_id++) {
    if ( ! selected[nhood_id] ) continue;
    running++;

    log << running << "." << n_selected << flush;

    triple<int,int,int> placement( -1, -1, -1 );
    const HyperBasevector &hyper = hypers[nhood_id];
    const look_align &al = seeds_aligns[ aligns_index[nhood_id] ];
    placement = EstimatePlacement( hyper, al, assembly_lengths );

    vec<look_align> aligns;
    RefinePlacement( placement, aligns, K, nhood_id, assembly_head,
		     out_dir, assembly_fasta, hypers, FORCE );
    
    int target = placement.first;
    int beg = placement.second;
    int end = placement.third;
    if ( target != range_id || end <= range_beg || range_end <= beg ) {
      log << "\n";
      continue;
    }

    log << "\tnhood " << nhood_id
	<< "\tu_" << seeds_ids[nhood_id]
	<< "\ton t_" << target
	<< " [" << placement.second
	<< ", " << placement.third
	<< ")_" << assembly_lengths[target]
	<< endl;
  }
  cout << endl;
  
}
