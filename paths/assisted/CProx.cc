/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/assisted/CProx.h"

/**
 * class CProx
 * Constructor
 */
CProx::CProx( const CProxInfo *info ) :
  info_ ( info ),
  lwin_ ( make_pair( 0, 0 ) ),
  rgap_ ( INT_MAX ),
  lgap_ ( INT_MAX ),
  ldev_ ( 0 )
{ }

/**
 * class CProx
 * AddRefGap
 */
void CProx::AddRefGap( const int gap )
{
  if ( gap < info_->MinGap( ) ) return;
  if ( gap > info_->MaxGap( ) ) return;
  rgap_ = gap;
}

/**
 * class CProx
 * AddLinksGaps
 *
 * Append gaps to existing list of link based gaps.
 */
void CProx::AddLinkGaps( const int tag_id, const vec<int> &gaps )
{
  gaps_.reserve( gaps_.size( ) + gaps.size( ) );
  for (size_t ii=0; ii<gaps.size( ); ii++)
    gaps_.push_back( make_pair( gaps[ii], tag_id ) );
}

/**
 * class CProx
 * ComputeLinkGap
 */
void CProx::ComputeLinkGap( )
{
  // HEURISTICS
  const double SIGMA_CLUSTER = 1.0;   // used to define consistent clusters

  // Defaults (in case of failure).
  lwin_ = make_pair( 0, 0 );
  lgap_ = INT_MAX;
  ldev_ = 0;

  // Nothing to do.
  if ( gaps_.size( ) < 1 ) return;
  
  // Sort gaps.
  sort( gaps_.begin( ), gaps_.end( ) );

  // Find clusters.
  vec<int> cluster_starts( 1, 0 );
  for (int ii=1; ii<gaps_.isize( ); ii++) {
    int gap1 = this->Gap( ii-1 );
    int gap2 = this->Gap( ii );
    int dev1 = int( SIGMA_CLUSTER * double( this->LibDev( ii-1 ) ) );
    int dev2 = int( SIGMA_CLUSTER * double( this->LibDev( ii ) ) );
    ho_interval win_prev( gap1 - dev1, gap1 + dev1 );
    ho_interval win_curr( gap2 - dev2, gap2 + dev2 );
    if ( Overlap( win_prev, win_curr ) < 1 )
      cluster_starts.push_back( ii );
  }
  
  // Build clusters, where clusters[ii] = ( count, window_in_gaps_ ).
  vec< pair<int, pair<int,int> > > clusters;
  for (int ii=0; ii<cluster_starts.isize( ); ii++) {
    bool is_last = ( ii == cluster_starts.isize( ) - 1 );
    int beg = cluster_starts[ii];
    int end = is_last ? gaps_.isize( ) : cluster_starts[ii+1];

    // Discard invalid gap.
    vec<NormalDistribution> nds;
    nds.reserve( end - beg );
    for (int ii=beg; ii<end; ii++) {
      float mu = this->Gap( ii );
      float sigma = this->LibDev( ii );
      nds.push_back( NormalDistribution( mu, sigma ) );
    }
    NormalDistribution result = CombineNormalDistributions( nds );
    int current_gap = int( result.mu_ );
    if ( current_gap < info_->MinGap( ) ) continue;
    if ( current_gap > info_->MaxGap( ) ) continue;
    if ( end - beg < info_->MinLinks( ) ) continue;
    
    // Ok, add valid cluster.
    clusters.push_back( make_pair( end - beg, make_pair( beg, end ) ) );
  }

  // No valid cluster found.
  if ( clusters.size( ) < 1 ) return;

  // Sort cluster, and find best one.
  sort( clusters.rbegin( ), clusters.rend( ) );
  
  int best_id = 0;
  if ( clusters.size( ) > 1 ) {
    int best_count = clusters[0].first;
    int best2_count = clusters[1].first;
    if ( best2_count <= best_count )
      best_id = -1;
  }
  if ( best_id < 0 ) return;

  // Compute gap.
  lwin_ = clusters[0].second;
  vec<NormalDistribution> nds;
  nds.reserve( lwin_.second - lwin_.first );
  for (int ii=lwin_.first; ii<lwin_.second; ii++) {
    float mu = this->Gap( ii );
    float sigma = this->LibDev( ii );
    nds.push_back( NormalDistribution( mu, sigma ) );
  }
  NormalDistribution result = CombineNormalDistributions( nds );
  lgap_ = int( result.mu_ );
  ldev_ = Max( 1, int( result.sigma_ ) );
  
}

/**
 * class CProx
 * EstimatedGap
 */
pair<int,int> CProx::EstimatedGap( ) const
{
  const int MIN_DEV = 10;

  if (clb_.Weight() > 0)
    return make_pair(clb_.Sep(), clb_.Dev());

  int ref_gap = this->RefGap( );
  pair<int,int> link_gap = this->LinkGap( );

  ForceAssert( ref_gap != INT_MAX || link_gap.first != INT_MAX );
  if ( link_gap.first == INT_MAX )
    return make_pair( ref_gap, Max( MIN_DEV, ref_gap / 4 ) );

  link_gap.second = Max( MIN_DEV, link_gap.second );
  return link_gap;
}

/**
 * class CProx
 * Print
 */
void CProx::Print( const bool brief,
		   const int v,
		   const int w,
		   ostream &out ) const
{
  bool skip = ( rgap_ == INT_MAX && lwin_.first == lwin_.second );
  if ( skip ) return;   // this should not happen
  
  out << v / 2 << ( v % 2 == 0 ? "[+]" : "[-]" ) << "  =>  "
      << w / 2 << ( w % 2 == 0 ? "[+]" : "[-]" ) << "\n";

  if ( rgap_ != INT_MAX ) 
    out << "  ref-based gap: " << rgap_ << "\n";

  if ( lwin_.first != lwin_.second )
    out << "  link-based gap: " << lgap_ << " +/- " << ldev_
	<< " (" << lwin_.second - lwin_.first << " links)\n";

  int beg = brief ? lwin_.first : 0;
  int end = brief ? lwin_.second : gaps_.isize( );
  for (int ii=beg; ii<end; ii++) {
    int gap = gaps_[ii].first;
    int tag = gaps_[ii].second;
    out << "    " << ii << "." << gaps_.isize( ) - 1 << "  "
	<< gap << " +/- " << info_->LibDev( tag ) << "  "
	<< "(" << info_->LibName( tag ) << ")\n";
  }

  out << "\n";
}

/**
 * class CProx
 * writeBinary
 */
void CProx::writeBinary( BinaryWriter& writer ) const
{
  writer.write( gaps_ );
  writer.write( lwin_ );
  writer.write( rgap_ );
  writer.write( lgap_ );
  writer.write( ldev_ );
  writer.write( clb_ );
}

/**
 * class CProx
 * readBinary
 */
void CProx::readBinary( BinaryReader& reader )
{
  reader.read( &gaps_ );
  reader.read( &lwin_ );
  reader.read( &rgap_ );
  reader.read( &lgap_ );
  reader.read( &ldev_ );
  reader.read( &clb_ );
}

