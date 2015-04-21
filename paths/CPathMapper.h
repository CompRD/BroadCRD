/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_PATH_MAPPER_H
#define C_PATH_MAPPER_H

#include "Basevector.h"
#include "Vec.h"
#include "paths/CKmerAlign.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"

/**
 * class CPathMapper
 *
 * Class to facilitate the mapping of segments or kmers onto a path
 * databse.
 */
class CPathMapper {

public:

  CPathMapper( );

  CPathMapper( const vecKmerPath *paths_fw,
	       const vecKmerPath *paths_rc,
	       const vec<tagged_rpint> *paths_db,
	       const vecbasevector *bases = 0 );

  void SetPaths( const vecKmerPath *paths_fw,
		 const vecKmerPath *paths_rc,
		 const vec<tagged_rpint> *paths_db );

  void SetBases( const vecbasevector *bases );

  // All segments (whole segments) matching kpi
  void Matches( const KmerPathInterval &kpi, vec<tagged_rpint> &answer ) const;
  
  // Find alignments of reads onto kpi
  void Matches( const KmerPathInterval &kpi, vec<CKmerAlign> &answer,
		const int t_id = 0, const int t_offset = 0 ) const;
    
  // Returns -1 if kmer_id is not in the read
  int KmerIdToPos( int read_id, longlong kmer_id, bool rc = false ) const;
  
  // Returns -1 if pos is out of range
  longlong PosToKmerId( int read_id, int pos, bool rc = false ) const;
  
  // operator=
  CPathMapper& operator= ( const CPathMapper &in );
  
  
private:
  
  const vecKmerPath *paths_fw_;       // paths...
  const vecKmerPath *paths_rc_;       // ...and rc paths
  const vec<tagged_rpint> *paths_db_; // the read path interval database
  const vecbasevector *bases_;        // bases (match paths_fw_, may be null)

};

#endif
