///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoverageAnalyzer.h"
#include "PairsManager.h"
#include "SeqInterval.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "paths/reporting/CHits.h"

/**
 * CHits
 * Constructor
 */
CHits::CHits( const vec<look_align> &hits,
	      const PairsManager &pairs,
	      const vec<int> &cglens,
	      ostream *log ) :
  hits_ ( hits ),
  pairs_ ( pairs ),
  cglens_ ( cglens ),
  log_ ( log )
{
  this->Cleanup( );
  this->GenerateMaps( );
  this->EstimateLibStats( );
}

/**
 * CHits
 * CoveragesAtMost
 *
 * Coverages are computed on ilens_ (whole inserts, rather than just
 * the separations between end points of reads).
 */
void CHits::CoveragesAtMost( const int lib_id,
			     const int max_cov,
			     vec<seq_interval> &regions ) const
{
  CoverageAnalyzer coverator( ilens_[lib_id], &cglens_ );
  coverator.GetCoveragesAtMost( max_cov, regions );
}
  
/**
 * CHits
 * CoveragesExactly
 *
 * As before, coverages are computed on ilens_ (whole inserts, rather
 * than just the separations between end points of reads).
 */
void CHits::CoveragesExactly( const int lib_id,
			      const int cov,
			      vec<seq_interval> &regions ) const
{
  CoverageAnalyzer coverator( ilens_[lib_id], &cglens_ );
  coverator.GetCoveragesExactly( cov, regions );
}
  
/**
 * CHits
 * CoveragesAtLeast
 *
 * As before, coverages are computed on ilens_ (whole inserts, rather
 * than just the separations between end points of reads).
 */
void CHits::CoveragesAtLeast( const int lib_id,
			      const int min_cov,
			      vec<seq_interval> &regions ) const
{
  CoverageAnalyzer coverator( ilens_[lib_id], &cglens_ );
  coverator.GetCoveragesAtLeast( min_cov, regions );
}
  
/**
 * CHits
 * UpdatePairsFile
 *
 * Change pairs file (and save a backup copy of the original) to
 * reflect lib size and stdev computed by EstimateLibStats( ).
 * Notice we need to reload the pairs, since pairs_ is const.
 *
 * NB. An existing backup file will not be overwritten!
 */
void CHits::UpdatePairsFile( const String &pairs_file, ostream *log ) const
{
  String backup_file = pairs_file + ".UpdatePairsFile.orig";

  PairsManager pairs( pairs_file );
  if ( IsRegularFile( backup_file ) ) {
    String str_mess = ": NOT saving backup pairs file (it already exists)";
    if ( log ) *log << Date( ) << str_mess << endl;
  }
  else {
    if ( log ) *log << Date( ) << ": saving backup pairs file" << endl;
    pairs.Write( backup_file );
  }

  vec< vec<String> > table;
  vec<String> line;
  if ( log ) {
    line = MkVec( String( "lib_name" ),
		  String( "orig_stats" ),
		  String( "new_stats" ) );
    table.push_back( line );
  }

  for (size_t ii=0; ii<pairs.nLibraries( ); ii++) {
    String lib_name = pairs.getLibraryName( ii );
    int old_sep = pairs.getLibrarySep( ii );
    int old_dev = pairs.getLibrarySD( ii );
    int new_sep = nds_[ii].mu_;
    int new_dev = nds_[ii].sigma_;
    pairs.changeLibrarySepSd( ii, new_sep, new_dev );
    if ( ! log ) continue;

    line.clear( );
    line.push_back( lib_name,
		    ToString( old_sep ) + " +/- " + ToString( old_dev ),
		    ToString( new_sep ) + " +/- " + ToString( new_dev ) );
    table.push_back( line );
  }

  if ( log ) {
    *log << "\n";
    PrintTabular( *log, table, 3, "rrr" );
    *log << endl;
  }

  if ( log ) *log << Date( ) << ": saving corrected pairs file" << endl;
  pairs.Write( pairs_file );
}
  
/**
 * CHits
 * Cleanup
 * private
 */
void CHits::Cleanup( )
{
  to_hit_.clear( );
  nds_.clear( );
  seps_.clear( );
  ilens_.clear( );
}

/**
 * CHits
 * GenerateMaps
 * private
 */
void CHits::GenerateMaps( )
{
  to_hit_.clear( );
  to_hit_.resize( pairs_.nReads( ), -1 );
  for (int ii=0; ii<hits_.isize( ); ii++) {
    int read_id = hits_[ii].query_id;
    if ( to_hit_[read_id] == -2 ) continue;
    else if ( to_hit_[read_id] == -1 ) to_hit_[read_id] = ii;
    else to_hit_[read_id] = -2;
  }
}

/**
 * CHits
 * FindSeparations
 * private
 *
 * If max_stretch is not null, use only valid inserts.
 */
void CHits::FindSeparations( const double *max_stretch )
{
  // Numbers of pairs.
  vec<size_t> lib_sizes = pairs_.getLibrarySizes( );
  
  // Reserve memory.
  seps_.clear( );
  seps_.resize( pairs_.nLibraries( ) );
  for (int lib_id=0; lib_id<(int)pairs_.nLibraries( ); lib_id++)
    seps_[lib_id].reserve( lib_sizes[lib_id] );
  
  if ( max_stretch ) {
    ilens_.clear( );
    ilens_.resize( pairs_.nLibraries( ) );
    for (int lib_id=0; lib_id<(int)pairs_.nLibraries( ); lib_id++)
      ilens_[lib_id].reserve( lib_sizes[lib_id] );
  }
  
  // Loop over all pairs.
  for (int pair_id=0; pair_id<(int)pairs_.nPairs( ); pair_id++) {
    int id1 = pairs_.ID1( pair_id );
    int id2 = pairs_.ID2( pair_id );
    if ( to_hit_[id1] < 0 || to_hit_[id2] < 0 ) continue;
    
    // Skip non fully embedded reads.
    if ( ! this->IsFullyEmbedded( hits_[ to_hit_[id1] ] ) ) continue;
    if ( ! this->IsFullyEmbedded( hits_[ to_hit_[id2] ] ) ) continue;
    
    // Different target.
    int t1 = hits_[ to_hit_[id1] ].target_id;
    int t2 = hits_[ to_hit_[id2] ].target_id;
    if ( t1 != t2 ) continue;
    
    // Same orientation.
    bool rc1 = hits_[ to_hit_[id1] ].Rc1( );
    bool rc2 = hits_[ to_hit_[id2] ].Rc1( );
    if ( rc1 == rc2 ) continue;
    
    // A good pair (unless max_stretch is defined, see test below).
    int hit_fw = rc2 ? to_hit_[id1] : to_hit_[id2];
    int hit_rc = rc2 ? to_hit_[id2] : to_hit_[id1];
    const look_align &alFw = hits_[hit_fw];
    const look_align &alRc = hits_[hit_rc];
    
    // This test is only run if max_stretch_ is given.
    int lib_id = pairs_.libraryID( pair_id );
    if ( max_stretch ) {
      int observed_sep = alRc.a.pos2( ) - alFw.a.Pos2( );
      double given_sep = nds_[lib_id].mu_;
      double given_dev = nds_[lib_id].sigma_;
      double stretch = double( observed_sep - given_sep ) / given_dev;
      if ( Abs( stretch ) > *max_stretch ) continue;
    }
    
    // Fill seps_.
    seq_interval sint( pair_id, t1, alFw.a.Pos2( ), alRc.a.pos2( ) );
    seps_[lib_id].push_back( sint );
    if ( ! max_stretch ) continue;
    
    // Fill ilens_.
    sint.Set( pair_id, t1, alFw.a.pos2( ), alRc.a.Pos2( ) );
    ilens_[lib_id].push_back( sint );
  }
  
}
  
/**
 * CHits
 * EstimateLibStats
 * private
 *
 * Crudely remove outliers (use guess if unable to estimate!)
 */
void CHits::EstimateLibStats( )
{
  // Heuristics!
  const int min_sep = 0;            // min allowed sep
  const int max_sep = 1000;         // max allowed sep
  const int min_to_call = 2000;     // min seps needed to estim stats
  const double guess_mu = 150;      // default guess fragment size
  const double guess_sigma = 100;   // default guess fragment stdev
  const double max_stretch = 5.0;   // max stretch allowed for final seps_
  
  // Compute seps_ (first pass).
  this->FindSeparations( );
  
  // Compute stats.
  nds_.clear( );
  nds_.reserve( pairs_.nLibraries( ) );
  for (int lib_id=0; lib_id<seps_.isize( ); lib_id++) {
    vec<float> lens;
    for (int ii=0; ii<seps_[lib_id].isize( ); ii++) {
      float len = seps_[lib_id][ii].Length( );
      if ( len < min_sep ) continue;
      if ( len > max_sep ) continue;
      lens.push_back( len );
    }
    sort( lens.begin( ), lens.end( ) );
    if ( lens.isize( ) < min_to_call ) {
      if ( log_ )
	*log_ << "lib_" << lib_id
	      << ": not enough seps! Guessing <" << ToString( guess_mu, 0 )
	      << " +/- " << ToString( guess_sigma ) << ">\n";
      nds_.push_back( NormalDistribution( guess_mu, guess_sigma ) );
      continue;
    }
    nds_.push_back( SafeMeanStdev( lens ) );
    if ( log_ )
      *log_ << "lib_" << lib_id
	    << ": <" << ToString( nds_.back( ).mu_, 0 )
	    << " +/- " << ToString( nds_.back( ).sigma_, 0 )
	    << "> (estimated on " << ToStringAddCommas( lens.size( ) )
	    << " inserts)\n";
  }
  if ( log_ ) *log_ << endl;
  
  // Recompute seps_, only allowing "reasonable" stretches.
  this->FindSeparations( &max_stretch );
}

/**
 * CHits
 * IsFullyEmbedded
 * private
 */
bool CHits::IsFullyEmbedded( const look_align &hit ) const
{
  if ( hit.a.pos1( ) != 0 ) return false;
  if ( hit.a.Pos1( ) != (int)hit.query_length ) return false;
  return true;
}

