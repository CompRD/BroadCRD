///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EXCISETRANSPOSONS
#define EXCISETRANSPOSONS

#include <fstream>
#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "Vec.h"
#include "VecString.h"

void ExciseTransposons( vecbasevector& reads, vecqualvector& Q, 
     const vecString& ids,
     String transposon_file, vec< vec<int>* >& killed_reads, vec<int>& start, 
     vec<int>& stop, ostream& log, String run_dir, int ml, 
     const vec<int>& min_overlap, const vec<int>& min_perfect_overlap,
     const vec<Bool>& kill_all, const Bool quiet = False );

#endif
