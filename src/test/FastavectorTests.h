////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file FastavectorTests.h
 * \author tsharpe
 * \date Oct 22, 2009
 *
 * \brief
 */
#ifndef TEST_FASTAVECTORTESTS_H_
#define TEST_FASTAVECTORTESTS_H_

#include <cxxtest/TestSuite.h>
#include "Fastavector.h"

class FastavectorTests : public CxxTest::TestSuite
{
public:
    void test_AllBasevectors()
    {
        fastavector fv1("A");
        vecbvec vbv1; vbv1.push_back(bvec("A"));
        TS_ASSERT(vbv1 == fv1.AllBasevectors(1));

        fastavector fv2("ANA");
        vecbvec vbv2;
        vbv2.push_back(bvec("AAA")).push_back(bvec("ACA"));
        vbv2.push_back(bvec("AGA")).push_back(bvec("ATA"));
        TS_ASSERT(vbv2 == fv2.AllBasevectors(4));
        TS_ASSERT(vecbvec() == fv2.AllBasevectors(3));

        fastavector fv3("SW");
        vecbvec vbv3;
        vbv3.push_back(bvec("CA")).push_back(bvec("GA"));
        vbv3.push_back(bvec("CT")).push_back(bvec("GT"));
        TS_ASSERT(vbv3 == fv3.AllBasevectors(4));
    }

    void test_setToSubOfSelf()
    {
        fastavector fv("XYZZY");
        fv.SetToSubOf(fv,1,-1);
        TS_ASSERT(fv == fastavector("YZZY"));
        fv.SetToSubOf(fv,2,1);
        TS_ASSERT(fv == fastavector("Z"));
    }
};

#endif /* TEST_FASTAVECTORTESTS_H_ */
