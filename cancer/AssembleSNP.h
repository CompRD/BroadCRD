/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// AssembleSNP.  Given a putative SNP, attempt to assemble a tiny region around
// it, to see if it really looks like a SNP.  Return True if so.
//
// The key problem is to define the region that is to be assembled.  We want
// this to be as small as possible.

#ifndef ASSEMBLE_SNP_H
#define ASSEMBLE_SNP_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"

Bool AssembleSNP( const int maxread, const vecbasevector& bases, 
     const vecqualvector& quals, const vec<look_align>& aligns, 
     const vecbasevector& ref, const int tig, const int pos, const char altbase, 
     const Bool verbose = False );

#endif
