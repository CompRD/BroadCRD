///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file AlignmentCalculatorTest.h
 * \author tsharpe
 * \date Jun 21, 2010
 *
 * \brief
 */
#ifndef ALIGNMENTCALCULATORTEST_H_
#define ALIGNMENTCALCULATORTEST_H_

#include "system/AlignmentCalculator.h"
#include <cxxtest/TestSuite.h>

class AlignmentCalculatorTests : public CxxTest::TestSuite
{
    class S1_1 {};
    class S1_0 { char mV; };
    class S2_0 { short mV; };
    class S4_0 { int mV; };
    class S8_7 { long mV1; char mV2[1]; };
    class S8_6 { long mV1; char mV2[2]; };
    class S8_5 { long mV1; char mV2[3]; };
    class S8_4 { long mV1; char mV2[4]; };
    class S8_3 { long mV1; char mV2[5]; };
    class S8_2 { long mV1; char mV2[6]; };
    class S8_1 { long mV1; char mV2[7]; };
    class S8_0 { long mV1; char mV2[8]; };
    class S4_3 { short mV1; int mV2; char mV3; };

    template <class T> void check( size_t alignment, size_t padding )
    { TS_ASSERT_EQUALS(AlignmentCalculator<T>::getAlignment(),alignment);
      TS_ASSERT_EQUALS(AlignmentCalculator<T>::getTailPadding(),padding); }

public:
    void test_all()
    {
        check<S1_1>(1,1);
        check<S1_0>(1,0);
        check<S2_0>(2,0);
        check<S4_0>(4,0);
        check<S8_0>(8,0);
        check<S8_1>(8,1);
        check<S8_2>(8,2);
        check<S8_3>(8,3);
        check<S8_4>(8,4);
        check<S8_5>(8,5);
        check<S8_6>(8,6);
        check<S8_7>(8,7);
        check<S4_3>(4,3);
    }
};

#endif /* ALIGNMENTCALCULATORTEST_H_ */
