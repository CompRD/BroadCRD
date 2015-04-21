/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header file: AlignUtils.h

   Utilities for working with alignments,
   especially generic utils for working with align classes that
   model the Align concept.

   @file
 */

#ifndef __INCLUDE_lookup_AlignUtils_h
#define __INCLUDE_lookup_AlignUtils_h

/// Sort by: target_id - begin_on_target (begin of aligned portion).
template <class Align>
struct order_align_TargetBegin
  : public binary_function<const Align&, const Align&, bool>
{
 public:
  bool operator() ( const Align &left, const Align &right ) const {
    return left.TargetId() == right.TargetId() ?
      left.StartOnTarget() < right.EndOnTarget( ) :
      left.target_id < right.target_id;
  }
};

template <class Align>
void MakeAlignIndex( const vec< Align >& aligns,
		     int nqueries,
		     vec< align_id_t >& alignsIdx );

template <class Align>
Bool AlignIsValid( const Align& align );


#endif
// #ifndef __INCLUDE_lookup_AlignUtils_h
