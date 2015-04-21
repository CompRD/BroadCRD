/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// StartBias.  Print out a bunch of stats regarding bias of read start points.
// This is at present a poorly-documented hodge-podge.

#ifndef START_BIAS_H
#define START_BIAS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "lookup/LookAlign.h"

String StartBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec<look_align>& aligns, 
     const Bool brief = False, String marked_reference_file = "",
     const Bool summary_form = False, const Bool silent = False );

String StartBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec< vec<int> >& starts, 
     const Bool brief = False, String marked_reference_file = "", 
     const Bool summary_form = False, const Bool silent = False );

String StartBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec<placement_mark>& places, 
     const Bool brief = False, String marked_reference_file = "", 
     const Bool summary_form = False, const Bool silent = False );

String StartBiasRelRaw( const vecbasevector& ref, const vecbasevector& rref, 
     const vec<placement_mark>& places1, const vec<placement_mark>& places2 );

#endif
