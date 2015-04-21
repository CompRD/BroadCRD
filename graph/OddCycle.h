// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef ODD_CYCLE_H
#define ODD_CYCLE_H

#include "CoreTools.h"

void MinimalOddCycle( const vec< vec<Bool> >& edge, vec<int>& odd_cycle );

void OddCycleVertices( const vec< vec<Bool> >& edge, vec<int>& to_delete );

#endif
