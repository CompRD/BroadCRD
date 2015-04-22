// Copyright (c) 2000, 2001 Whitehead Institute of Biomedical Research
// 


#ifndef PADDEDSMITHWAT
#define PADDEDSMITHWAT

#include "Alignment.h"
#include "Basevector.h"
#include "PaddedBasevector.h"

float PaddedSmithWat( const vec<char>& S, const vec<char>& T, 
                      int offset, int bandwidth, alignment& a, 
                      const int mismatch_penalty = 3, 
                      const int gap_penalty = 6,
                      const float divider = 4 );

#endif
