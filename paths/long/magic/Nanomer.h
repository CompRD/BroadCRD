///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Jan 26, 2015 - <crdhelp@broadinstitute.org>
//

// The one million + 1th class to hold a kmer as two bits per base
// and handle the various bit twiddling associated with conversion
// to basevector, etc.

// We explicitly don't spend cycles checking for overflow of the requested
// underlying datatype


#ifndef NANOMER_H
#define NANOMER_H
#include "feudal/BaseVec.h"

template <unsigned int K, class STYPE = unsigned short>
class Nanomer {
public:
    Nanomer(STYPE val) : mBits(val) {};

    BaseVec ToBaseVec() {
        STYPE val = mBits;
        BaseVec kmer(K);
        for ( size_t i = 0; i < K; ++i ) {
            kmer.Set(K - i - 1, val & BASEMASK );
            val >>= BASEBITS;
        }
        return kmer;
    }

    String ToString() { return this->ToBaseVec().ToString(); }

private:
    static const unsigned int BASEBITS = 2u;
    static const STYPE BASEMASK = 0x3;
    STYPE mBits;
};

#endif  // NANOMER_H
