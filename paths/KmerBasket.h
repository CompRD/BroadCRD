///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// KmerBasket.  Assemble a bunch of reads.

#ifndef KMER_BASKET_H
#define KMER_BASKET_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"

void KmerBasket( const vecbasevector& R, const vecqualvector& Q, 
     HyperBasevector& hb, const int K, const int max_mismatches, 
     const int verbosity, ostream& out );

#endif
