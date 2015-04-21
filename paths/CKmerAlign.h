/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_KMER_ALIGN_H
#define C_KMER_ALIGN_H

#include "Vec.h"
#include "paths/KmerPath.h"

/**
 * class CKmerAlign
 *
 * Class to deal with alignments seeded by overlapping KmerPathIntervals.
 */
class CKmerAlign {

public:

  CKmerAlign( );
  
  CKmerAlign( const KmerPathInterval &kpi );
  
  KmerPathInterval Kpi( ) const { return kpi_; }

  int Length( ) const { return len_; }

  int QueryId( ) const { return qid_; }

  int TargetId( ) const { return tid_; }

  int QStart( ) const { return qstart_; } // align loc on query KmerPath

  int TStart( ) const { return tstart_; } // align loc on target KmerPath

  bool Rc( ) const { return qrc_; }
  
  // Begin, End are provided to be used with the Channeler class.
  int Begin( ) const { return tstart_; }
  
  int End  ( ) const { return tstart_ + len_; }
  
  // Print a one-liner with basic info
  void PrintCoreInfo( ostream &out ) const;

  // Set from read path interval, and (properly oriented) path
  void Set( const tagged_rpint &rpint, const KmerPath &path,
	    const int t_id = 0, const int t_offset = 0 );
  
  // Count how many aligns in ll and rr share same query_id_'s.
  friend int Intersect( const vec<CKmerAlign> &ll, const vec<CKmerAlign> &rr );
  
  // operator<
  friend bool operator< ( const CKmerAlign &a1, const CKmerAlign &a2 );
  
  
private:

  KmerPathInterval kpi_; // the target path interval
  int qid_;              // id of query path
  int tid_;              // id of target path
  int qstart_;           // where align starts on query path
  int tstart_;           // where align starts on target path
  bool qrc_;             // orientation of query on target
  int len_;              // length of alignment

};

#endif
