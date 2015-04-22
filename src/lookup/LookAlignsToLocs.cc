/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "ReadLocation.h"
#include "lookup/LookAlign.h"
#include "lookup/LookAlignsToLocs.h"

/**
 * LookAlignsToLocs
 */
void LookAlignsToLocs( vec<look_align_plus> &hits,
		       vec<read_location> &locs,
		       double *max_cov,
		       int *min_reads,
		       int *max_gap,
		       ostream *log )
{
  locs.clear( );

  // Sort hits if needed.
  order_lookalign_TargetLookAlignOffset sorterTO;
  sort( hits.begin( ), hits.end( ), sorterTO );
 
  // Reserve space.
  locs.reserve( hits.size( ) );

  // Generate contigs.
  int start_hit_id = 0;
  int c_id = 0;
  while( start_hit_id < (int)hits.size( ) ) {

    // Check if some gap between reads is allowed.
    int slack = max_gap ? *max_gap : 0;

    // Detect reads in contig.
    const look_align_plus &start_hit = hits[start_hit_id];
    int begin_base = LookAlignOffset( start_hit );
    int end_base = begin_base + start_hit.query_length;
    int end_hit_id = start_hit_id + 1;
    for (int ii=start_hit_id+1; ii<(int)hits.size( ); ii++) {
      if ( hits[ii].target_id != start_hit.target_id )
	break;
      if ( LookAlignOffset( hits[ii] ) > end_base + slack )
	break;
      int begin_hit = LookAlignOffset( hits[ii] );
      int end_hit = begin_hit + hits[ii].query_length;
      begin_base = Min( begin_base, begin_hit );
      end_base = Max( end_base, end_hit );
      end_hit_id++;
    }

    int c_len = end_base - begin_base;
    ForceAssert( c_len > 0 );

    // Skip contig?
    bool skip = false;

    int tot_len = 0;
    for (int ii=start_hit_id; ii<end_hit_id; ii++)
      tot_len += hits[ii].query_length;
    float cov = SafeQuotient( tot_len, c_len );
    if ( max_cov && cov > *max_cov )
      skip = true;
    
    int n_reads = end_hit_id - start_hit_id;
    if ( min_reads && n_reads < *min_reads )
      skip = true;
    
    if ( log )
      *log << "CONTIG_" << c_id
	   << "  len: " << c_len
	   << "  reads: " << n_reads
	   << "  coverage: " << ToString( cov, 2 )
	   << ( skip ? " (skip contig)" : "" )
	   << "\n";
      
    // Add reads to contig and log read info.
    int contig_id = ( locs.size( ) < 1 ) ? 0 : 1 + locs.back( ).Contig( );
    for (int hit_id=start_hit_id; hit_id<end_hit_id; hit_id++) {
      const look_align_plus &hit = hits[hit_id];
      int offset_hit = LookAlignOffset( hit );
      int t_id = hit.target_id;
      int r_id = hit.query_id;
      int r_len = hit.query_length;
      int c_start = offset_hit - begin_base;
      orientation orient = hit.rc1 ? ReverseOr : ForwardOr;
      String str_orient = hit.rc1 ? "rc" : "fw";
      
      if ( log )
	*log << " c_" << contig_id
	     << "_" << hit_id - start_hit_id
	     << "." << end_hit_id - start_hit_id
	     << "  r_" << r_id
	     << "  t_" << t_id
	     << " [" << offset_hit
	     << ", " << offset_hit + r_len 
	     << ") " << (int)orient
	     << ( skip ? " skipped" : "" )
	     << "\n";

      if ( !skip ) {
	read_location loc( r_id, r_len, contig_id, c_start, orient, c_len );
	locs.push_back ( loc );
      }
    }    
    
    // Next.
    start_hit_id = end_hit_id;
    c_id++;
  }

  // Sort locations.
  sort( locs.begin( ), locs.end( ) );

}

/**
 * LookAlignsToGapLocs
 */
void LookAlignsToGapLocs( const vec<int> &tlens,
			  vec<look_align_plus> &hits,
			  vec<read_location> &locs,
			  ostream *log )
{
  locs.clear( );
  
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  // Sort hits if needed.
  order_lookalign_TargetLookAlignOffset sorterTO;
  sort( hits.begin( ), hits.end( ), sorterTO );
 
  // Reserve, create fhits.
  locs.reserve( hits.size( ) );

  vec<int> fhits( tlens.size( ), -1 );
  for (int ii=hits.size( )-1; ii>=0; ii--)
    fhits[ hits[ii].target_id ] = ii;

  // Loop over all targets.
  for (int target_id=0; target_id<(int)tlens.size( ); target_id++) {
    int fhit = fhits[target_id];
    int count = 0;
    if ( fhit > -1 ) {
      for (int ii=fhit; ii<(int)hits.size( ); ii++) {
	if ( hits[ii].target_id != target_id ) break;
	count++;
      }
    }

    out << "TARGET_" << target_id
	<< " (" << tlens[target_id]
	<< " bp, " << count
	<< " hits)\n";
    
    // No hits on this target.
    if ( count < 1 ) {
      out << "\n";
      continue;
    }

    // Loop over hits on target.
    for (int hit_id=fhit; hit_id<(int)hits.size( ); hit_id++) {
      if ( hits[hit_id].target_id != target_id ) break;

      const look_align_plus &hit = hits[hit_id];
      int c_start = LookAlignOffset( hit );
      int c_id = hit.target_id;
      int r_id = hit.query_id;
      int r_len = hit.query_length;
      int c_len = tlens[target_id];
      orientation orient = hit.rc1 ? ReverseOr : ForwardOr;
      String str_orient = hit.rc1 ? "rc" : "fw";
      
      read_location loc( r_id, r_len, c_id, c_start, orient, c_len );
      locs.push_back ( loc );

      out << " c_" << c_id
	  << "_" << hit_id - fhit
	  << "." << count - 1
	  << "  r_" << r_id
	  << " [" << c_start
	  << ", " << c_start + r_len 
	  << ") " << (int)orient
	  << "\n";
    }

  }
  
}
