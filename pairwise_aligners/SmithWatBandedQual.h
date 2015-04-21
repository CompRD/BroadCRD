///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SMITH_WAT_BANDED_QUAL
#define SMITH_WAT_BANDED_QUAL

// This is the same as SmithWatBandedA, except that we use quality scores here.

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"

// Warning.  An object of type X contains the error count at a given position in
// the dynamic programming matrix.  Below, we use X = unsigned char, which could
// easily overflow.  Use of a larger type is recommended but carries with it a
// performance penalty.

// Note also the extremely ugly inclusion of the dummy argument x, which is needed,
// for unknown reasons.

template<class X> float SmithWatBandedQual2( const basevector& S, 
     const basevector& T, const qualvector& SQ, int offset, int bandwidth, 
     align& a, int& errors, ostream* log = 0, int mis=2, int ins=3, int del=3 );

inline float SmithWatBandedQual( const basevector& S, const basevector& T, 
     const qualvector& SQ, int offset, int bandwidth, align& a, int& errors, 
     ostream* log = 0, int mis=2, int gap=3  )
{    
     return SmithWatBandedQual2<unsigned char>( 
          S, T, SQ, offset, bandwidth, a, errors, log, mis, gap, gap  );    }

inline float SmithWatBandedQual( const basevector& S, const basevector& T, 
     const qualvector& SQ, int offset, int bandwidth, alignment& a, 
     ostream *log = 0 )
{    align temp; int errors;
     float result 
          = SmithWatBandedQual( S, T, SQ, offset, bandwidth, temp, errors, log );
     a.Set( temp, errors );
     return result;    }

#endif
