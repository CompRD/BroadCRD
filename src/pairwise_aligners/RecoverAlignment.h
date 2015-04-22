// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef RECOVERALIGNMENT
#define RECOVERALIGNMENT

#include "Alignment.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "system/Types.h"
#include "Vec.h"

Bool RecoverAlignment( int id1, int id2, const vecbasevector& EE,
     const vecqualvector& Q, Bool rc, alignment_plus& ap,
     int kmer_size = 8, int min_offset = -1000000000, int max_offset = 1000000000 );

#endif
