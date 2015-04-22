///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef KMERS__KMERS_TAGGER_CORE__H
#define KMERS__KMERS_TAGGER_CORE__H

#include "Basevector.h"
#include "Bitvector.h"
#include "Qualvector.h"
#include "kmers/ReadPather.h"
#include "kmers/ReadPatherDefs.h"

/**
 * KmersTaggerCore
 *
 * Error correction based on kmer frequency. Tag regions of reads with
 * frequency >= MIN_FREQ, and generate in output a vecbvec of valid
 * sequences, defined as all stretches of tagged kmers (discard
 * sequences < MIN_READLEN).
 *
 * MIN_FREQ: discard kmers with frequency < MIN_FREQ
 * MIN_READLEN: do not save short segments (kmer length)
 * VERBOSE: toggled verbose mode
 */
template <unsigned int K>
void KmersTaggerCore( const int MIN_FREQ,
		      const int MIN_READLEN,
		      const bool VERBOSE,
		      const KmerDict<K> &dict,
		      VirtualMasterVec<bvec> &bases,
		      vecbvec &all_segments,
		      ostream &log )
{
  all_segments.clear( );
  size_t n_reads = bases.size( );
  
  // Create a vecbitvector matching reads.
  vecbitvector is_valid;
  is_valid.reserve( n_reads );
  for (size_t ii=0; ii<n_reads; ii++) {
    int klen = Max( 0, (int)bases[ii].size( ) - (int)K + 1 );
    is_valid.push_back( BitVec( klen , false ) );
  }
  
  // Loop over all reads.
  String str_n_reads = ToStringAddCommas( n_reads );
  size_t n_dotter = 100000;
  size_t n_done = 0;
  log << Date( ) << ": digesting "
      << str_n_reads << " reads (. = "
      << ToStringAddCommas( n_dotter ) << " reads)"
      << endl;
  
  for (size_t rid=0; rid<bases.size( ); rid++) {
    if ( rid % n_dotter == 0 ) Dot( log, rid / n_dotter );
    
    const bvec &read = bases[rid];
    if ( read.size( ) < K ) continue;

    KMer<K> theKmer( read.begin( ) );
    bvec::const_iterator itr( read.begin( ) + K );
    while ( itr != read.end( ) + 1 ) {
      KDef const *pKDef = dict.lookup( theKmer );
      if ( pKDef ) {
	size_t pos = distance( read.begin( ), itr - K );
	bool is_freq = ( (int)pKDef->getCount( ) >= MIN_FREQ );
	is_valid[rid].Set( pos, is_freq );
      }
      
      theKmer.toSuccessor( *itr );
      itr++;
    }
  }
  log << endl;

  // Generating list of valid segments (find count in the first pass).
  log << Date( ) << ": generating valid segments" << endl;
  if ( VERBOSE ) log << "\n";
  size_t n_segments = 0;
  for (int pass=0; pass<2; pass++) {
    if ( pass == 1 ) all_segments.reserve( n_segments );

    for (size_t rid=0; rid<bases.size( ); rid++) {
      const int klen = bases[rid].size( ) - K + 1;

      // Log info.
      if ( VERBOSE && pass == 0 ) {
	for (size_t ii=0; ii<is_valid[rid].size( ); ii++)
	  log << ( is_valid[rid][ii] ? '.' : 'x' ); 
	log << "   r" << rid;
      }
      
      // Loop over kmers in read.
      int pos = 0;
      while ( pos < klen ) {
	
	// Go to start of valid segment.
	while ( pos < klen && ( ! is_valid[rid][pos] ) ) pos++;
	if ( pos >= klen ) break;

	// Find end of valid segment.
	int end = pos + 1;
	while ( end < klen && is_valid[rid][end] ) end++;
	
	// Add segment, if long enough.
	bool keeper = ( end - pos >= MIN_READLEN );
	if ( keeper ) {
	  if ( pass == 0 ) n_segments++;
	  else all_segments.push_back( bvec( bases[rid], pos, end-pos+K-1 ) );
	}
	
	// Log info.
	if ( VERBOSE && pass == 0 )
	  log << ( keeper ? "   [" : "   (" )
	      << pos << ", " << end
	      << ( keeper ? "]" : ")" );
	
	// Move on.
	pos = end;
      }

      // Close log line.
      if ( VERBOSE && pass == 0 ) log << "\n";
    }
  }
  if ( VERBOSE ) log << "\n";
  
  // Done.
  log << Date( ) << ": KmersTaggerCore done" << endl;
  
}

#endif
