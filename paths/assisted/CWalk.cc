///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Intvector.h"
#include "paths/KmerPath.h"
#include "paths/assisted/CWalk.h"

/**
 * class CWalk
 * Constructor
 */
CWalk::CWalk( const longlong gid,
	      const KmerPath &kpath,
	      const termination_code termin )
{
  this->Extend( gid, kpath, 0, &termin );
}

/**
 * class CWalk
 * Clear
 */
void CWalk::Clear( )
{
  gids_.clear( );
  govers_.clear( );
  kpath_.Clear( );
  termin_ = NOT_TERMINATED;
}

/**
 * class CWalk
 * Reserve
 */
void CWalk::Reserve( const size_t n_extensions )
{
  gids_.reserve( n_extensions );
  govers_.reserve( n_extensions );
}

/**
 * class CWalk
 * Extend
 */
void CWalk::Extend( const longlong gid,
		    const KmerPath &ext,
		    const int *over,
		    const termination_code *termin )
{
  kpath_ = ext;
  gids_.push_back( gid );
  if ( over ) govers_.push_back( *over );
  if ( termin ) termin_ = *termin;
  ForceAssertEq( gids_.size( ) - 1, govers_.size( ) );
}

/**
 * class CWalk
 * Extend
 */
void CWalk::PrintGidsInfo( const vecKmerPath &cloud_paths, ostream &out ) const
{
  String str_termin = "";
  if ( termin_ == NOT_TERMINATED ) str_termin = "NOT_TERMINATED";
  else if ( termin_ == TARGET_FOUND ) str_termin = "TARGET_FOUND";
  else if ( termin_ == NO_EXTENSIONS ) str_termin = "NO_EXTENSION";
  else if ( termin_ == TOO_LONG ) str_termin = "TOO_LONG";
  else ForceAssert( 1 == 0 );
  
  int u1_len = cloud_paths[gids_.front( )].TotalLength( );
  int u2_len = cloud_paths[gids_.back( )].TotalLength( );
  int closure_len = kpath_.TotalLength( ) - u1_len;
  if ( termin_ == TARGET_FOUND ) closure_len += -u2_len;
  
  int n_steps = gids_.size( ) - 1;
  if ( termin_ == TARGET_FOUND ) n_steps += -1;
  
  out << closure_len << " kmers, " << n_steps << " steps";
  for (size_t ii=0; ii<gids_.size( ); ii++) {
    bool spec = ( ii==0 || ( ii==gids_.size( )-1 && termin_==TARGET_FOUND ) );
    char openbra = ( spec ) ? '{' : '[';
    char closebra = ( spec ) ? '}' : ']';
    int clen = cloud_paths[ gids_[ii] ].TotalLength( );
    
    out << "   " << openbra << gids_[ii] << closebra << "." << clen;
    if ( ii < gids_.size( ) - 1 )
      out << "   (" << govers_[ii] << ")";
  }
  
  out << "   " << str_termin << "\n";
  
}
