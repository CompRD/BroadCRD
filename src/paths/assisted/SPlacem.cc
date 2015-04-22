///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "PairsManager.h"
#include "SupersHandler.h"
#include "VecUtilities.h"
#include "paths/ReadLoc.h"
#include "paths/assisted/SPlacem.h"

/**
 * SPlacem
 * Constructor
 */
SPlacem::SPlacem( ) :
  rid_( -1 ),
  cid_( -1 ),
  sid_( -1 ),
  beg_( 0 ),
  end_( 0 ),
  rc_( true )
{ }

/**
 * SPlacem
 * Constructor
 */
SPlacem::SPlacem( int64_t rid,
		  int32_t cid,
		  int sid,
		  int beg,
		  int end,
		  bool rc ) :
  rid_( rid ),
  cid_( cid ),
  sid_( sid ),
  beg_( beg ),
  end_( end ),
  rc_( rc )
{ }

/**
 * SPlacem
 * FindAlign
 */
int SPlacem::FindAlign( const vec< triple<int64_t,int64_t,int> > &aligns ) const
{ 
  triple<int64_t,int64_t,int> bogus( rid_, -1 , 0 );
  vec< triple<int64_t,int64_t,int> >::const_iterator it
    = lower_bound( aligns.begin( ), aligns.end( ), bogus );
  ForceAssert( it != aligns.end( ) );
  ForceAssert( it->first == rid_ );
  return distance ( aligns.begin( ), it );
}

/**
 * SPlacem
 * IsLogicalPair
 */
bool SPlacem::IsLogicalPair( const int MIN_SEP,
			     const int MAX_SEP,
			     const SPlacem &other ) const
{
  if ( this->sid_ != other.sid_ ) return false;
  if ( this->rc_ == other.rc_ ) return false;

  const SPlacem &fw_plac = this->rc_ ? other : *this;
  const SPlacem &rc_plac = this->rc_ ? *this : other;
  int sep = rc_plac.beg_ - fw_plac.end_;
  if ( sep < MIN_SEP || sep > MAX_SEP ) return false;

  return true;
}

/**
 * SPlacem
 * AddReadLoc
 */
void SPlacem::AddReadLoc( const vec< triple<int64_t,int64_t,int> > &aligns,
			  const vec<int> &clens,
			  const PairsManager &pairs,
			  vec<read_loc> &locs ) const
{
  const int idx = this->FindAlign( aligns );
  const triple<int64_t,int64_t,int> &al = aligns[idx];
  
  int64_t pairid = pairs.getPairID( rid_ );
  int64_t p_id =
    pairs.ID1( pairid ) == rid_
    ? pairs.ID2( pairid )
    : pairs.ID1( pairid );
  int32_t cid = cid_ < 0 ? - cid_ - 1 : cid_;
  int32_t p_cid = -1;   // unmapped partner
  uint16_t rlen = end_ - beg_;
  int initpos = al.third < 0 ? - al.third - 1 : al.third;
  int cpos = ( cid < 0 ) ? clens[cid] - initpos - rlen : initpos;
  int p_cpos = -1;   // unmapped partner
  Bool fw = ( ( al.third < 0 ) == ( cid_ < 0 ) );
  Bool p_fw = false;   // unmapped partner
  int8_t rclass = 1;   // read class = jump by default
  uint8_t lib_id = pairs.libraryID( pairid );
  uint16_t p_rlen = 0;   // unmapped partner
  uint8_t band = 0;    // bandwidth - undefined!
  uint8_t p_band = 0;    // bandwidth - undefined!
  bool p_placed = false;   // unmapped partner
  
  locs.push_back( read_loc( rid_, p_id, cid, p_cid, cpos, p_cpos,
			    fw, p_fw, rclass, lib_id, rlen, p_rlen,
			    band, p_band, p_placed ) );

}

/**
 * SPlacem
 * AddReadLocs
 */
void SPlacem::AddReadLocs( const vec< triple<int64_t,int64_t,int> > &aligns,
			   const vec<int> &clens,
			   const PairsManager &pairs,
			   const SPlacem &partner,
			   vec<read_loc> &locs ) const
{
  const int idx = this->FindAlign( aligns );
  const triple<int64_t,int64_t,int> &al = aligns[idx];
  
  const int p_idx = partner.FindAlign( aligns );
  const triple<int64_t,int64_t,int> &p_al = aligns[p_idx];
  
  int64_t pairid = pairs.getPairID( rid_ );
  int64_t p_id =
    pairs.ID1( pairid ) == rid_
    ? pairs.ID2( pairid )
    : pairs.ID1( pairid );
  int32_t cid = cid_ < 0 ? - cid_ - 1 : cid_;
  int32_t p_cid = partner.cid_< 0 ? - partner.cid_ - 1 : partner.cid_;
  uint16_t rlen = end_ - beg_;
  int initpos = al.third < 0 ? - al.third - 1 : al.third;
  int cpos = ( cid < 0 ) ? clens[cid] - initpos - rlen : initpos;
  uint16_t p_rlen = partner.end_ - partner.beg_;
  int p_initpos = p_al.third < 0 ? - p_al.third - 1 : p_al.third;
  int p_cpos = ( p_cid < 0 ) ? clens[p_cid] - p_initpos - p_rlen : p_initpos;
  Bool fw = ( ( al.third < 0 ) == ( cid_ < 0 ) );
  Bool p_fw = ( ( p_al.third < 0 ) == ( partner.cid_ < 0 ) );
  int8_t rclass = 1;   // read class = jump by default
  uint8_t lib_id = pairs.libraryID( pairid );
  uint8_t band = 0;    // bandwidth - undefined!
  uint8_t p_band = 0;    // bandwidth - undefined!
  bool p_placed = true;
  
  locs.push_back( read_loc( rid_, p_id, cid, p_cid, cpos, p_cpos,
			    fw, p_fw, rclass, lib_id, rlen, p_rlen,
			    band, p_band, p_placed ) );
  
  locs.push_back( read_loc( p_id, rid_, p_cid, cid, p_cpos, cpos,
			    p_fw, fw, rclass, lib_id, p_rlen, rlen,
			    p_band, band, p_placed ) );
  
}

