// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
// 

#ifndef PRINT_LOOK_ALIGN_H
#define PRINT_LOOK_ALIGN_H

#include "Basevector.h"
#include "lookup/LookAlign.h"

// Print the given look_align_plus
void PrintLookAlignPlus( ostream &out,
			 const vecbasevector &target_bases,
			 const vecbasevector &query_bases,
			 const look_align_plus &look_alplus,
			 bool abbreviate = true );

#endif
