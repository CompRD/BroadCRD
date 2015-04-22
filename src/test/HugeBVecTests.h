///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file HugeBVecTests.h
 * \author tsharpe
 * \date Mar 14, 2012
 *
 * \brief
 */
#ifndef HUGEBVECTESTS_H_
#define HUGEBVECTESTS_H_

#include <cxxtest/TestSuite.h>
#include "feudal/HugeBVec.h"
#include "Basevector.h"

class HugeBVecTests : public CxxTest::TestSuite
{
public:
    void test_append()
    {
        HugeBVec hbv;
        bvec bv1("A");
        bvec bv2("AC");
        bvec bv3("ACG");
        bvec bv4("ACGT");
        hbv.append(bv1.begin(),bv1.end());
        hbv.append(bv2.begin(),bv2.end());
        hbv.append(bv3.begin(),bv3.end());
        hbv.append(bv4.begin(),bv4.end());
        bvec bv5("AACACGACGT");
        HugeBVec hbv2;
        hbv2.resize(10);
        for ( unsigned idx = 0; idx < 10; ++idx )
            hbv2.set(idx,bv5[idx]);
        TS_ASSERT_EQUALS(hbv,hbv2);
    }
};

#endif /* HUGEBVECTESTS_H_ */
