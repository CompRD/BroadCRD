///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "STLExtensions.h"
#include "util/StroboAligner.h"

void StroboAligner( vec< triple<int64_t,int64_t,int> > &aligns,
		    const String &TARGET,
		    const String &QUERY,
		    const String &OUTDIR,
		    const int K,
		    const int STROBE,
		    const int MAX_GAP,
		    const int MIN_LENGTH,
		    const double MIN_RATIO )
{
  aligns.clear( );

  String kmers_file = OUTDIR + "/kmers.fastb";
  String log_file = OUTDIR + "/StroboAligner.log";
  
  Mkpath( OUTDIR );
  ofstream log( log_file.c_str( ) );

  log << Date( ) << ": loading query" << endl;
  vecbvec query( QUERY );

  log << Date( ) << ": selecting kmers" << endl;
  int nkmers = query.size( );
  for (int ii=0; ii<(int)query.size( ); ii++)
    nkmers += query[ii].isize( ) / STROBE;   // approximation for reserve
  
  vecbvec kmers;
  vec< pair<int,int> > kmers_map;
  kmers.reserve( nkmers );
  kmers_map.reserve( nkmers );

  vec<int> nk( query.size( ), 0 );   // exact count of kmers for each read
  for (int ii=0; ii<(int)query.size( ); ii++) {
    int pos = 0;
    int count = 0;
    while ( pos + K < (int)query[ii].size( ) ) {
      kmers.push_back( bvec( query[ii], pos, K ) );
      kmers_map.push_back( make_pair( ii, pos ) );
      pos += STROBE;
      nk[ii] += 1;
      count++;
    }
  }
  
  log << Date( ) << ": saving selected kmers" << endl;
  kmers.WriteAll( kmers_file );

  log << Date( ) << ": aligning kmers" << endl;
  vec< triple<int64_t,int64_t,int> > kmer_aligns;
  SearchFastb2( kmers_file, TARGET, K, &kmer_aligns );

  vec< vec< pair<int,int> > > placements( query.size( ) );
  
  log << Date( ) << ": parsing alignments" << endl;
  for (size_t id=0; id<kmer_aligns.size( ); id++) {
    const triple<int64_t,int64_t,int> &al = kmer_aligns[id];
    
    // Fw aligns only.
    if ( al.third < 0 ) continue;
    
    // Only care for long queries.
    int kmer_id = al.first;
    int query_id = kmers_map[kmer_id].first;
    if ( (int)query[query_id].size( ) < MIN_LENGTH ) continue;
    
    int target_id = al.second;
    int kmer_start = al.third;
    int kmer_offset = kmers_map[kmer_id].second;
    int query_start = kmer_start - kmer_offset;
    
    placements[query_id].push_back( make_pair( target_id, query_start ) );
  }
  for (size_t ii=0; ii<placements.size( ); ii++)
    sort( placements[ii].begin( ), placements[ii].end( ) );
  
  log << Date( ) << ": digesting alignments" << endl;
  for (int query_id=0; query_id<(int)query.size( ); query_id++) {
    const vec< pair<int,int> > &placs = placements[query_id];
    int query_len = query[query_id].size( );
    int query_weight = nk[query_id];
    
    // Clusters of placements.
    int start = 0;
    int pos = 1;
    while ( pos < placs.isize( ) ) {
      int end = -1;
      if ( pos == placs.isize( )-1 ) end = pos;
      else if ( placs[pos-1].first != placs[pos].first ) end = pos;
      else if ( placs[pos].second - placs[pos-1].second > MAX_GAP ) end = pos;
      if ( end > -1 ) {
	int weight = 1 + end - start;
	if ( (double)weight / (double)query_weight > MIN_RATIO ) {
	  triple<int64_t,int64_t,int> 
	    newal( query_id, placs[start].first, placs[start].second );
	  aligns.push_back( newal );
	  log << query_id << " (l=" << query_len << "): t"
	      << placs[start].first << " @"
	      << placs[start].second << "   w="
	      << weight << " ("
	      << query_weight << " strobes total)\n";
	}
	start = pos;
      }
      pos++;
    }
  }
  
  log << Date( ) << ": StroboAligner done" << endl;
  log.close( );

}
