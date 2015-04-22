/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "paths/CKmerAlign.h"
#include "paths/KmerPath.h"

/**
 * CKmerAlign
 * Constructor
 */
CKmerAlign::CKmerAlign( ) :
  qid_ ( -1 ),
  tid_ ( -1 ),
  qstart_ ( -1 ),
  tstart_ ( -1 ),
  qrc_ ( false ),
  len_ ( 0 )
{ }

/**
 * CKmerAlign
 * Constructor
 */
CKmerAlign::CKmerAlign( const KmerPathInterval &kpi ) :
  kpi_ ( kpi ),
  qid_ ( -1 ),
  tid_ ( -1 ),
  qstart_ ( -1 ),
  tstart_ ( -1 ),
  qrc_ ( false ),
  len_ ( 0 )
{ }

/**
 * CKmerAlign
 * PrintCoreInfo
 */
void CKmerAlign::PrintCoreInfo( ostream &out ) const
{
  out << "[" << kpi_.Start( ) << "," << kpi_.Stop( ) << "]"
      << " q" << qid_
      << ( qrc_ ? "_rc" : "" )
      << ( tstart_ > 0 ? " +" + ToString( tstart_ ) : "" )
      << " @" << qstart_
      << " l=" << len_
      << "\n";
}

/**
 * CKmerAlign
 * Set
 */
void CKmerAlign::Set( const tagged_rpint &rpint, const KmerPath &path,
		      const int t_id, const int t_offset )
{
  // Given query path interval does not overlap kpi_ (leave).
  longlong alpha = Max( kpi_.Start( ), rpint.Start( ) );
  longlong beta = Min( kpi_.Stop( ), rpint.Stop( ) );
  if ( alpha > beta ) return;
  
  len_ = 1 + beta - alpha;
  tid_ = t_id;
  
  // Query sequence's id and orientation flag.
  qid_ = rpint.ReadId( );
  qrc_ = rpint.Rc( );
  
  // Start on target and query.
  tstart_ = Max( 0, int( rpint.Start( ) - kpi_.Start( ) ) );
  qstart_ = Max( 0, int( kpi_.Start( ) - rpint.Start( ) ) );
  
  // Adjust the start locations for the offsets of path intervals on their paths
  for (int ii=0; ii<rpint.PathPos( )-1; ii++)
    qstart_ += path.Length( ii );
  tstart_ += t_offset;
}

/**
 * CKmerAlign
 * Intersect
 *
 * It only looks for matching query id_'s.
 */
int Intersect( const vec<CKmerAlign> &ll, const vec<CKmerAlign> &rr )
{
  vec<int> qids;
  qids.reserve( rr.size( ) );
  for (int ii=0; ii<rr.isize( ); ii++)
    qids.push_back( rr[ii].QueryId( ) );
  sort( qids.begin( ), qids.end( ) );

  int nshared = 0;
  for (int ii=0; ii<ll.isize( ); ii++)
    if ( binary_search( qids.begin( ), qids.end( ), ll[ii].QueryId( ) ) )
      nshared++;
  
  return nshared;
}

/**
 * CKmerAlign
 * operator<
 */
bool operator< ( const CKmerAlign &a1, const CKmerAlign &a2 )
{
  if( a1.Kpi( )      != a2.Kpi( ) )      return ( a1.Kpi( )      < a2.Kpi( ) );
  if( a1.TargetId( ) != a2.TargetId( ) ) return ( a1.TargetId( ) < a2.TargetId( ) );
  if( a1.TStart( )   != a2.TStart( ) )   return ( a1.TStart( )   < a2.TStart( ) );
  if( a1.QueryId( )  != a2.QueryId( ) )  return ( a1.QueryId( )  < a2.QueryId( ) );
  if( a1.QStart( )   != a2.QStart( ) )   return ( a1.QStart( )   < a2.QStart( ) );
  if( a1.Rc( )       != a2.Rc( ) )       return ( a1.Rc( )       < a2.Rc( ) );
  
  return ( a1.Length( ) < a2.Length( ) );
}

