///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file KmerTests.h
 * \author tsharpe
 * \date Apr 5, 2012
 *
 * \brief
 */
#ifndef KMERTESTS_H_
#define KMERTESTS_H_

#include <cxxtest/TestSuite.h>
#include "kmers/KMer.h"
#include "Basevector.h"

class KMerTests : public CxxTest::TestSuite
{
public:
    void test_rc()
    {
      bvec b1("ACGGCTCTTGCATATAACGGTCTTA");
      KMer<25> k1(b1.begin());
      KMer<25> k2(b1.rcbegin());
      TS_ASSERT_EQUALS(k1.rc(),k2);
    }

};

#endif /* KMERTESTS_H_ */
