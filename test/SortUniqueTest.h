////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file SortUniqueTest.h
 * \author tsharpe
 * \date Oct 30, 2009
 *
 * \brief
 */
#ifndef SORTUNIQUETEST_H_
#define SORTUNIQUETEST_H_

#include <cxxtest/TestSuite.h>
#include "STLExtensions.h"
#include <vector>

class SortUniqueTest : public CxxTest::TestSuite
{
public:
    void test_sort_unique()
    { using std::vector;
      vector<int> v1;
      v1.push_back(3);
      v1.push_back(1);
      v1.push_back(4);
      v1.push_back(1);
      v1.push_back(5);
      v1.push_back(9);
      sort_unique(v1);
      vector <int> v2;
      v2.push_back(1);
      v2.push_back(3);
      v2.push_back(4);
      v2.push_back(5);
      v2.push_back(9);
      TS_ASSERT_EQUALS(v1,v2);
    }

};

#endif /* SORTUNIQUETEST_H_ */
