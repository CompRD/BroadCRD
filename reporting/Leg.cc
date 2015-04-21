// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "LocsHandler.h"
#include "PairsHandler.h"
#include "String.h"
#include "SupersHandler.h"
#include "reporting/Leg.h"

/*
 * leg
 * Set
 */
void leg::Set( int loc_id,
	       const shandler &supers,
	       const phandler &pairs,
	       const lhandler &locs )
{
  // Default.
  loc_id_ = loc_id;
  obs_length_ = numeric_limits<int>::quiet_NaN( );
  stretch_ = numeric_limits<float>::quiet_NaN( );
  status_ = uninitialized;

  // Partner missing.
  int read_id = locs[loc_id].ReadId( );
  int part_id = pairs.GetPartnerId( read_id );
  if ( part_id < 0 ) {
    status_ = partner_missing;
    return;
  }
  
  // Partner unassembled or multiply assembled.
  vec<int> all_part_locs = locs.GetAllLocationIds( part_id );
  if ( all_part_locs.size( ) < 1 ) {
    status_ = partner_unassembled;
    return;
  }
  if ( all_part_locs.size( ) > 1 ) {
    status_ = partner_multiply_assembled;
    return;
  }

  // Partner in other super.
  int read_contig = locs[loc_id].Contig( );
  int read_super = supers.ToSuper( read_contig );

  int part_loc_id = all_part_locs[0];
  int part_contig = locs[part_loc_id].Contig( );
  int part_super = supers.ToSuper( part_contig );

  if ( read_super != part_super ) {
    status_ = partner_in_other_super;
    return;
  }
  
  // Read and partner illogically oriented.
  orientation read_or = locs[loc_id].OrientationOnContig( );
  orientation part_or = locs[part_loc_id].OrientationOnContig( );
  if ( read_or == part_or ) {
    status_ = partner_illogical;
    return;
  }

  // Ok, this is a logical link: find stretch.
  status_ = partner_logical;
  
  const read_location &locA = locs[loc_id];
  const read_location &locB = locs[part_loc_id];
  const read_location &loc_fw = read_or == ForwardOr ? locA : locB;
  const read_location &loc_rc = read_or == ForwardOr ? locB : locA;
  
  int fw_cg = loc_fw.Contig( );
  int fw_cg_start = supers.StartOnSuper( fw_cg );
  int fw_end = fw_cg_start + loc_fw.StartOnContig( ) + loc_fw.LengthOfRead( );

  int rc_cg = loc_rc.Contig( );
  int rc_cg_start = supers.StartOnSuper( rc_cg );
  int rc_beg = rc_cg_start + loc_rc.StartOnContig( );
  
  int observ_len = rc_beg - fw_end;

  int pair_id = pairs.GetPairId( read_id );
  int given_len = pairs[pair_id].sep;
  int given_stdev = pairs[pair_id].sd;

  int len1 = loc_fw.LengthOfRead( );
  int len2 = loc_rc.LengthOfRead( );

  obs_length_ = observ_len + len1 + len2;
  stretch_ = SafeQuotient( ( observ_len - given_len ), given_stdev );

}
  
/*
 * leg
 * InsertType
 *
 * It will return the insert length closer to one of the given types. For
 * example, say inserts = (4, 10, 40) and actual insert length = 3750; then
 * InsertType will return "4Kb". Notice "inserts" has to be given in Kb.
 */
String leg::InsertType( const lhandler &locs,
			const phandler &pairs,
			const vec<int> inserts ) const
{
  ForceAssert( inserts.size( ) > 0 );
  
  String str_type;

  if ( uninitialized == status_ )
    str_type = "uninitialized";
  else { 
    int read_id = locs[loc_id_].ReadId( );
    int pair_id = pairs.GetPairId( read_id );
    if ( pair_id < 0 )
      str_type = "na";
    else {
      int size = pairs[pair_id].sep;
      int cont = 0;
      int dist = Abs( size - 1000 * inserts[0] );
      for (int ii=1; ii<(int)inserts.size( ); ii++) {
	int newdist = Abs( size - 1000 * inserts[ii] );
	if ( newdist < dist ) {
	  cont = ii;
	  dist = newdist;
	}
      }
      str_type = ToString( inserts[cont] ) + "Kb";
    }
  }

  return str_type;
}

/*
 * leg
 * StatusAsString
 *
 * The returned String refers to the partner read of the leg (eg
 * missing is short for partner_missing). If max_stretch is not null,
 * then a logical insert is defined valid iff Abs( stretch ) <=
 * *max_stretch.
 */
String leg::StatusAsString( const float *max_stretch ) const
{
  String str_status;

  switch( status_ ) {
  case uninitialized :
    str_status = "uninitialized";
    break;
  case partner_missing :
    str_status = "missing";
    break;
  case partner_unassembled :
    str_status = "unassembled";
    break;
  case partner_multiply_assembled :
    str_status = "multiply_assembled";
    break;
  case partner_in_other_super :
    str_status = "in_other_super";
    break;
  case partner_illogical :
    str_status = "illogical";
    break;
  case partner_logical :
    if ( max_stretch )
      str_status = Abs( stretch_ ) <= *max_stretch ? "valid" : "invalid";
    else 
      str_status = "logical (stretch=" + ToString( stretch_, 2 ) + ")";
    break;
  default :
    str_status = "error";
  }

  return str_status;
}
