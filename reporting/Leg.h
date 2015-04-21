// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Describe the status of one (assembled) leg of an insert.

#ifndef LEG_H
#define LEG_H

#include "LocsHandler.h"
#include "PairsHandler.h"
#include "String.h"
#include "SupersHandler.h"

// Set of valid stati.
enum leg_status {
  uninitialized,
  partner_missing,
  partner_unassembled,
  partner_multiply_assembled,
  partner_in_other_super,
  partner_illogical,
  partner_logical
};

// Class containing leg status for one assembled read.
class leg {

public:
  
  leg( ) :
    loc_id_ ( -1 ), stretch_ ( 0.0 ), status_ ( uninitialized ) { }

  leg( int loc_id, float stretch, leg_status status ) :
    loc_id_ ( loc_id ), stretch_ ( stretch ), status_ ( status ) { }

  leg( int loc_id,
       const shandler &supers,
       const phandler &pairs,
       const lhandler &locs ) { this->Set( loc_id, supers, pairs, locs ); }
  
  void Set( int loc_id,
	    const shandler &supers,
	    const phandler &pairs,
	    const lhandler &locs );
  
  int LocId( ) const { return loc_id_; }
  
  int ObservedInsertLength( ) const { return obs_length_; }

  float Stretch( ) const { return stretch_; }

  leg_status Status( ) const { return status_; }

  String InsertType( const lhandler &locs,
		     const phandler &pairs,
		     const vec<int> inserts ) const;

  String StatusAsString( const float *max_stretch = 0 ) const;
  
  
private:

  int loc_id_;          // id of loc
  int obs_length_;      // observed insert length (iff status = logical)
  float stretch_ ;      // stretch (iff status = logical)
  leg_status status_;   // leg status

};

#endif
