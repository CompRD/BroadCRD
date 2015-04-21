/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "LocsHandler.h"
#include "SeqInterval.h"
#include "Vec.h"
#include "util/RemoveLocsOverWindows.h"

/**
 * RemoveLocsOverWindows
 */
void RemoveLocsOverWindows( const int algo_select,
			    const vec<seq_interval> &wins,
			    const lhandler &locs_in,
			    vec<read_location> &locs_out,
			    ostream *plog )
{
  ForceAssert( algo_select == 1 );
  ForceAssert( is_sorted( wins.begin( ), wins.end( ) ) );

  ofstream devnull ( "/dev/null" );
  ostream &log = plog ? *plog : devnull;

  int n_contigs = locs_in.FirstLocs( ).size( );
  vec<int> fwins( n_contigs, -1 );
  for (int ii=wins.size( )-1; ii>=0; ii--)
    fwins[ wins[ii].SeqId( ) ] = ii;
  
  // Select locs to be removed.
  vec<Bool> to_remove( locs_in.Size( ), False );
  
  // Loop over contigs.
  for (int contig_id=0; contig_id<n_contigs; contig_id++) {
    int fwin = fwins[contig_id];
    if ( fwin < 0 )
      continue;
    
    // Build the set of bad intervals.
    vec< pair<int,int> > bad;
    for (int ii=fwin; ii<(int)wins.size( ); ii++) {
      if ( wins[ii].SeqId( ) != contig_id ) break;
      bad.push_back( make_pair( wins[ii].Begin( ), wins[ii].End( ) ) );
    }
    
    // Loop over locs in contig.
    int floc = locs_in.FirstLoc( contig_id );
    if ( floc < 0 ) continue;
    for (int ii=floc; ii<locs_in.Size( ); ii++) {
      if ( locs_in[ii].Contig( ) != contig_id ) break;
      for (int jj=0; jj<(int)bad.size( ); jj++) {
	int right = Min( locs_in[ii].StopOnContig( ) + 1, bad[jj].second );
	int left = Max( locs_in[ii].StartOnContig( ), bad[jj].first );

	// Remove loc.
	if ( algo_select == 1 )
	  if ( left < right )
	    to_remove[ii] = True;
	
	// Report removal.
	if ( to_remove[ii] ) {
	  log << "removing c_" << contig_id
	      << "/r_" << locs_in[ii].ReadId( )
	      << " [" << locs_in[ii].StartOnContig( )
	      << "," << locs_in[ii].StopOnContig( ) + 1
	      << "): it overlaps [" << bad[jj].first
	      << "," << bad[jj].second
	      << ")\n";
	  break;
	}
	
      } // next bad interval
    } // next loc
  } // next contig
  log << endl;
  
  // Generate output locs.
  int n_del = 0;
  int n_tot = locs_in.Size( );
  longlong tot_length = 0;
  longlong tot_del = 0;

  locs_out.reserve( locs_in.Size( ) );
  for (int ii=0; ii<locs_in.Size( ); ii++) {
    tot_length += locs_in[ii].LengthOfRead( );
    if ( to_remove[ii] ) {
      tot_del += locs_in[ii].LengthOfRead( );
      n_del++;
      continue;
    }
    locs_out.push_back( locs_in[ii] );
  }

  // Final stats.
  float n_ratio = 100.0 * SafeQuotient( n_del, n_tot );
  float len_ratio = 100.0 * SafeQuotient( tot_del, tot_length );

  log << "SUMMARY:\n"
      << " reads deleted: " << n_del
      << " out of " << n_tot
      << " (" << ToString( n_ratio, 2 )
      << "%)\n"
      << " bases deleted: " << tot_del
      << " out of " << tot_length
      << " (" << ToString( len_ratio, 2 )
      << "%)\n"
      << endl;    
  
}

