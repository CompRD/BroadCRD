/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header file: AlignUtilsTemplate.h

   Utilities for working with alignments,
   especially generic utils for working with align classes that
   model the Align concept.

   @file
 */

#ifndef __INCLUDE_lookup_AlignUtils_h
#define __INCLUDE_lookup_AlignUtils_h

#include "lookup/LookAlign.h"

template <class Align>
void MakeAlignIndex( const vec< Align >& aligns,
		     int nqueries,
		     vec< align_id_t >& alignsIdx ) {
  ForceAssertLe( 0, nqueries );
  alignsIdx.resize( nqueries, -1 );
  for ( align_id_t alignId = 0; alignId < aligns.isize(); alignId++ ) {
    const Align& a = aligns[ alignId ];
    ForceAssertLe( 0, a.QueryId() );
    ForceAssertLt( a.QueryId(), nqueries );
    alignsIdx[ a.QueryId() ] = alignId;
  }
}

template <class Align>
Bool AlignIsValid( const Align& a ) {
  return a.QueryId() >= 0 &&
    a.TargetId() >= 0 &&
    a.StartOnQuery() >= 0 &&
    a.EndOnQuery() >= 0 &&
    a.StartOnTarget() >= 0 &&
    a.EndOnTarget() >= 0;
}

#define DEFINE_ALIGN_UTILS(T) \
  template \
    void MakeAlignIndex<T>( const vec< T >& aligns,	\
	  		   int nqueries, \
			   vec< align_id_t >& alignsIdx ); \
  template \
    Bool AlignIsValid<T>( const T& a )
   

#endif
// #ifndef __INCLUDE_lookup_AlignUtils_h
