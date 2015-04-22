/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef VEC_T_ALIGNS_H
#define VEC_T_ALIGNS_H

#include "PackAlign.h"
#include "String.h"
#include "tiled/TAlign.h"

// Core bin write (it assumes all aligns belong to the same contig id).
void BinWriteTAligns( ostream &out, int id, const vec<t_align> &t_aligns );

// Core bin read (it loads all aligns for one contig, saving contig_id in id).
void BinReadTAligns( istream &in, int &id, vec<t_align> &t_aligns );

// Load from file all TAlignPluses.
void LoadAllTAlignPluses( const String &in_file, vec<t_align_plus> &t_aplus );

#endif
