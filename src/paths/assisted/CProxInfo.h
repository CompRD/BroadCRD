/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_PROX_INFO__H
#define PATHS__ASSISTED__C_PROX_INFO__H

#include "PairsManager.h"
#include "String.h"

/**
 * class CProxInfo
 *
 * Ancillary metadata and heuristics for gaps informations stored as
 * CProx objects.
 */
class CProxInfo {

public:
  
  CProxInfo( ) { this->SetDefaults( ); }

  CProxInfo( const int min_gap, const int max_gap, const int min_links ) :
    MIN_GAP_ ( min_gap ), MAX_GAP_ ( max_gap ), MIN_LINKS_ ( min_links )
  { }

  void SetDefaults( ) { MIN_GAP_ = -2500; MAX_GAP_ = 5000; MIN_LINKS_ = 2; }
  
  void SetMinGap( int min_gap ) { MIN_GAP_ = min_gap; }
  void SetMaxGap( int max_gap ) { MAX_GAP_ = max_gap; }
  void SetMinLinks ( int min_links ) { MIN_LINKS_ = min_links; }

  void SetLibTags( const PairsManager *jumps = 0,
		   const PairsManager *Jumps = 0 ) {
    if ( jumps ) {
      for (size_t ii=0; ii<jumps->nLibraries( ); ii++) {
	lib_names_.push_back( jumps->getLibraryName( ii ) );
	lib_types_.push_back( 0 );
	lib_ids_.push_back( int( ii ) );
	lib_sizes_.push_back( jumps->getLibrarySep( ii ) );
	lib_devs_.push_back( jumps->getLibrarySD( ii ) );
      }
    }
    if ( Jumps ) {
      for (size_t ii=0; ii<Jumps->nLibraries( ); ii++) {
	lib_names_.push_back( Jumps->getLibraryName( ii ) );
	lib_types_.push_back( 1 );
	lib_ids_.push_back( int( ii ) );
	lib_sizes_.push_back( Jumps->getLibrarySep( ii ) );
	lib_devs_.push_back( Jumps->getLibrarySD( ii ) );
      }
    }
  }
  
  int MinGap( ) const { return MIN_GAP_; }
  int MaxGap( ) const { return MAX_GAP_; }
  int MinLinks( ) const { return MIN_LINKS_; }

  // Returns -1 on error.
  int TagId( const int lib_type, const int lib_id ) const {
    int tag_id = -1;
    for (int ii=0; ii<lib_names_.isize( ); ii++) {
      if ( lib_types_[ii] != lib_type ) continue;
      if ( lib_ids_[ii] != lib_id ) continue;
      tag_id = ii;
      break;
    }
    return tag_id;
  }

  String LibName( const int tag ) const { return lib_names_[tag]; }
  int LibType( const int tag ) const { return lib_types_[tag]; }
  int LibId( const int tag ) const { return lib_ids_[tag]; }
  int LibSize( const int tag ) const { return lib_sizes_[tag]; }
  int LibDev( const int tag ) const { return lib_devs_[tag]; }
  
  
private:
  
  int MIN_GAP_;           // min allowed gap size
  int MAX_GAP_;           // max allowed gap size
  int MIN_LINKS_;         // min number of links to join two contigs
  
  vec<String> lib_names_; // name of library, as from the pairs file
  vec<int> lib_types_;    // 0 (jump) or 1 (long jump)
  vec<int> lib_ids_;      // id of library in the pairs file
  vec<int> lib_sizes_;    // lab-given estimate for library size...
  vec<int> lib_devs_;     // ... and its dev
  
};

#endif
