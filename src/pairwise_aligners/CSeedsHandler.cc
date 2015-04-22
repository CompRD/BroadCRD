/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "ReadLocation.h"
#include "String.h"
#include "STLExtensions.h"
#include "pairwise_aligners/CSeedsHandler.h"

/**
 * Class CSeedsHandler
 * Constructor
 */
CSeedsHandler::CSeedsHandler( int read_id ) :
  read_id_ ( read_id )
{
  this->SetDefaults( );
}

/**
 * Class CSeedsHandler
 * Add
 */
void CSeedsHandler::Add( int target_id, int start, Bool rc )
{
  CSeed seed( target_id, start, rc );
  seeds_.push_back( seed );
}

/**
 * Class CSeedsHandler
 * ClusterSeeds
 */
void CSeedsHandler::ClusterSeeds( )
{
  if ( ! is_sorted( seeds_.begin( ), seeds_.end( ) ) )
    sort( seeds_.begin( ), seeds_.end( ) );
  fseeds_.clear( );
  for (int ii=seeds_.size( )-1; ii>0; ii--) {
    if ( seeds_[ii-1].TargetId( ) != seeds_[ii].TargetId( ) )
      fseeds_.push_back( ii );
    else if ( seeds_[ii-1].RC( ) != seeds_[ii].RC( ) )
      fseeds_.push_back( ii );
    else if ( seeds_[ii].Start( ) - seeds_[ii-1].Start( ) > max_diff_ )
      fseeds_.push_back( ii );
  }
  if ( seeds_.size( ) > 0 )
    fseeds_.push_back( 0 );
  for (int ii=0; ii<(int)fseeds_.size( )/2; ii++)
    swap( fseeds_[ii], fseeds_[fseeds_.size( ) - 1 - ii] );
}

/**
 * Class CSeedsHandler
 * HeavySeed
 *
 * Select cluster with more seeds (if two top clusters have the same
 * of seeds, it will pick one of the two), and package an average
 * seed. Return false if nothing is found.
 */
bool CSeedsHandler::HeavySeed( CSeed &heavy_seed ) const
{
  ForceAssert( seeds_.size( ) == 0 || fseeds_.size( ) > 0 );
    
  // No seeds.
  heavy_seed.Set( -1, 0, false );
  if ( seeds_.size( ) < 1 )
    return false;

  // Select bag with highest weight.
  int last = fseeds_.size( ) - 1;
  vec< pair<int,int> > weight2fseedpos;
  for (int ii=0; ii<(int)fseeds_.size( ); ii++) {
    int beg = fseeds_[ii];
    int end = ii == last ? (int)seeds_.size( ) : fseeds_[ii+1];
    weight2fseedpos.push_back( make_pair( end - beg, ii ) );
  }
  sort( weight2fseedpos.rbegin( ), weight2fseedpos.rend( ) );

  // Fill heavy_seed and return.
  heavy_seed = this->AverageSeed( weight2fseedpos[0].second );
  return true;
}

/**
 * Class CSeedsHandler
 * AverageSeed
 */
CSeed CSeedsHandler::AverageSeed( int fseedpos ) const
{
  ForceAssert( fseedpos < (int)fseeds_.size( ) );
    
  int last = (int)fseeds_.size( ) - 1;
  int beg = fseeds_[fseedpos];
  int end = fseedpos == last ? (int)seeds_.size( ) : fseeds_[fseedpos+1];
    
  vec<int> starts;
  for (int ii=beg; ii<end; ii++)
    starts.push_back( seeds_[ii].Start( ) );

  int target = seeds_[beg].TargetId( );
  int start = (int)Mean( starts );
  bool rc = seeds_[beg].RC( );

  CSeed result( target, start, rc );
  return result;
}

/**
 * Class CSeedsHandler
 * PrintInfo
 *
 * For every cluster it prints the average of the seeds in the cluster.
 */
void CSeedsHandler::PrintInfo( ostream &out ) const
{
  ForceAssert( seeds_.size( ) == 0 || fseeds_.size( ) > 0 );
  out << "r_" << read_id_ << " " << fseeds_.size( ) << "_clusters:";
  for (int ii=0; ii<(int)fseeds_.size( ); ii++) {
    int n_entries
      = ( ii == (int)fseeds_.size( ) - 1 )
      ? seeds_.size( ) - fseeds_[ii]
      : fseeds_[ii+1] - fseeds_[ii];
    out << " " << n_entries << "x";
    this->AverageSeed( ii ).Print( out );
  }
  out << "\n";
}

/**
 * Class CSeedsHandler
 * SetDefaults
 */
void CSeedsHandler::SetDefaults( )
{
  max_diff_ = 12;
}

