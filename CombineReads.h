// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef COMBINE_READS
#define COMBINE_READS

#include <map>

#include "Basevector.h"
#include "Qualvector.h"
#include "system/Types.h"
#include "Vec.h"

Bool CombineReads( const vecbasevector& reads, 
     vecqualvector& reads_q, vec<int>& keep,
     vecbasevector& contigs, vecqualvector& contigs_q,
     map<int,int>& reads_to_contigs, vec<int>& reads_to_pos2,
     vec<int>& reads_to_Pos2, vec<int>& reads_RC, ostream& log,
     int k = 24, Bool ok_to_call_phrap = True, int max_badness = 100,
     int max_clique = 1000, int max_badness_for_makealigns = 1000,
     int max_errors_for_makealigns = 200 );

#endif
