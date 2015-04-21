// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "CoreTools.h"
#include "ReadLocation.h"
#include "testing/TraceTest2.h"

int EvalVec( const vec<int>& v, int n )
{    return v[n];    }

int EvalVec2( const vec<int>& v, int n )
{    read_location r;
     r.SetStartOnContig( v[n] );
     return -1;    }
