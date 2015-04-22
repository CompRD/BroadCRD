///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * KMerHasherTests.h
 *
 *  Created on: Dec 10, 2013
 *      Author: tsharpe
 */

#ifndef KMERHASHERTESTS_H_
#define KMERHASHERTESTS_H_

#include <cxxtest/TestSuite.h>
#include "Basevector.h"
#include "kmers/KMerHasher.h"

class KMerHasherTests : public CxxTest::TestSuite
{
public:
    void testStepping()
    {
        bvec bv("ACGTACGTACGT");
        KMerHasher<5> oddHasher;
        auto itr = bv.begin();
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.hash(itr));
        ++itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepF(itr));
        ++itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepF(itr));
        ++itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepF(itr));
        ++itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepF(itr));
        ++itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepF(itr));
        ++itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepF(itr));
        ++itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepF(itr));
        --itr;
        TS_ASSERT_EQUALS(oddHasher(itr),oddHasher.stepR(itr));

        itr = bv.begin();
        BuzHasher<5> bh1;
        BuzHasher<5> bh2;
        uint64_t hVal = bh1.hash(itr);
        ++itr;
        TS_ASSERT_EQUALS(bh2.hash(itr),bh1.step(itr,hVal));
    }

    void testRC()
    {
        bvec bv1("ACGTACGTACG");
        bvec bv2("CGTACGTACGT");
        KMerHasher<11> h;
        TS_ASSERT_EQUALS(h(bv1.begin()),h(bv2.begin()));

        BuzHasher<11> b;
        TS_ASSERT_DIFFERS(b.hash(bv1.begin()),b.hash(bv2.begin()));
    }
};



#endif /* KMERHASHERTESTS_H_ */
