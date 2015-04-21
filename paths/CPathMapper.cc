/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Vec.h"
#include "paths/CKmerAlign.h"
#include "paths/CPathMapper.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"

/**
 * CPathMapper
 * Constructor
 */
CPathMapper::CPathMapper( ) :
  paths_fw_ ( 0 ),
  paths_rc_ ( 0 ),
  paths_db_ ( 0 ),
  bases_ ( 0 )
{ }

/**
 * CPathMapper
 * Constructor
 */
CPathMapper::CPathMapper( const vecKmerPath *paths_fw,
			  const vecKmerPath *paths_rc,
			  const vec<tagged_rpint> *paths_db,
			  const vecbasevector *bases ) :
  paths_fw_ ( paths_fw ),
  paths_rc_ ( paths_rc ),
  paths_db_ ( paths_db ),
  bases_ ( bases )
{ }

/**
 * CPathMapper
 * SetPaths
 */
void CPathMapper::SetPaths( const vecKmerPath *paths_fw,
			    const vecKmerPath *paths_rc,
			    const vec<tagged_rpint> *paths_db )
{
  paths_fw_ = paths_fw;
  paths_rc_ = paths_rc;
  paths_db_ = paths_db;
}

/**
 * CPathMapper
 * SetBases
 */
void CPathMapper::SetBases( const vecbasevector *bases )
{
  bases_ = bases;
}

/**
 * CPathMapper
 * Matches
 */
void CPathMapper::Matches( const KmerPathInterval &kpi,
			   vec<tagged_rpint> &answer ) const
{
  answer.clear( );
  vec<longlong> pint_ids;
  Contains( *paths_db_, kpi, pint_ids );
  answer.reserve( pint_ids.size( ) );
  for (int ii=0; ii<(int)pint_ids.size( ); ii++)
    answer.push_back( (*paths_db_)[ pint_ids[ii] ] );
}

/**
 * CPathMapper
 * Matches
 */
void CPathMapper::Matches( const KmerPathInterval &kpi,
			   vec<CKmerAlign> &answer,
			   const int t_id,
			   const int t_offset ) const
{
  answer.clear( );
  vec<tagged_rpint> rpints;
  this->Matches( kpi, rpints );

  answer.resize( rpints.size( ), CKmerAlign( kpi ) );
  
  for (int ii=0; ii<(int)rpints.size( ); ii++) {
    read_id_t id = rpints[ii].ReadId( );
    bool rc = rpints[ii].Rc( );
    const KmerPath &path = rc ? (*paths_rc_)[id] : (*paths_fw_)[id];

    answer[ii].Set( rpints[ii], path, t_id, t_offset );
  }
}

/**
 * CPathMapper
 * KmerIdToPos
 */
int CPathMapper::KmerIdToPos( int read_id, longlong kmer_id, bool rc ) const
{
  const KmerPath &path = rc ? (*paths_rc_)[read_id] : (*paths_fw_)[read_id];

  int pos = 0;
  for (int ii=0; ii<(int)path.NSegments( ); ii++) {
    if ( path.Start( ii ) <= kmer_id && kmer_id <= path.Stop( ii ) )
      return pos + ( kmer_id - path.Start( ii ) );
    else
      pos += path.Length( ii );
  }

  return -1;
}

/**
 * CPathMapper
 * PosToKmerId
 */
longlong CPathMapper::PosToKmerId( int read_id, int pos, bool rc ) const
{
  const KmerPath &path = rc ? (*paths_rc_)[read_id] : (*paths_fw_)[read_id];

  int cur = 0;
  for (int ii=0; ii<(int)path.NSegments( ); ii++) {
    if ( cur <= pos && pos <= path.Length( ii ) )
      return path.Start( ii ) + (longlong)( pos - cur );
    else
      cur += path.Length( ii );
  }

  return -1;
}

/**
 * CPathMapper
 * operator=
 */
CPathMapper &CPathMapper::operator= ( const CPathMapper &in )
{
  if ( this == &in ) return *this;
  paths_fw_ = in.paths_fw_;
  paths_rc_ = in.paths_rc_;
  paths_db_ = in.paths_db_;
  bases_ = in.bases_;

  return *this;
}

