// Copyright (c) 2003 Whitehead Institute for Biomedical Research
//

#ifndef LIN_MIN_SOLVE_H
#define LIN_MIN_SOLVE_H

#include "CoreTools.h"
#include "math/Matrix.h"

float DiscrepToSolveNonneg( const matrix<float>& A, const vec<float>& b,
     Bool no_special = False );

#endif
