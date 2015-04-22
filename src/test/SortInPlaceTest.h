///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SortInPlaceTest.h
 * \author tsharpe
 * \date Feb 24, 2010
 *
 * \brief
 */
#ifndef SORTINPLACETEST_H_
#define SORTINPLACETEST_H_

#include <cxxtest/TestSuite.h>
#include "Intvector.h"
#include "TaskTimer.h"
#include "VecString.h"
#include "random/Random.h"
#include "system/SortInPlace.h"
#include <algorithm>
#include <vector>

class SortInPlaceTest : public CxxTest::TestSuite
{
public:
    void testOrder( ULongVec::iterator itr1,
                    ULongVec::iterator itr2 )
    {
        while ( ++itr1 != itr2 )
            TS_ASSERT_LESS_THAN_EQUALS(itr1[-1],*itr1);
    }

    void test_SortInPlaceLotsOfEqualElements()
    {
        size_t const nElts = 5000000;
        size_t const maxInstances = 500;
        ULongVec vec;
        vec.reserve(nElts);
        while ( vec.size() < vec.capacity() )
        {
            size_t val = randomx()%nElts;
            size_t nnn = randomx()%maxInstances + 1;
            while ( vec.size() < vec.capacity() && nnn-- )
                vec.push_back(val);
        }
        size_t const nPermutations = nElts;
        for ( size_t jjj = 0; jjj < nPermutations; ++jjj )
        {
            size_t idx1 = randomx()%nElts;
            size_t idx2 = randomx()%nElts;
            using std::swap; swap(vec[idx1],vec[idx2]);
        }
        sortInPlace(vec.begin(),vec.end());
        testOrder(vec.begin(),vec.end());
    }

    void test_SortInPlaceRandomUnique()
    {
        size_t const nElts = 10000;
        ULongVec vec;
        vec.reserve(nElts);
        for ( size_t iii = 0; iii < nElts; ++iii )
            vec.push_back(iii);

        size_t const nTests = 1000;
        for ( size_t iii = 0; iii < nTests; ++iii )
        {
            size_t const nPermutations = nElts;
            for ( size_t jjj = 0; jjj < nPermutations; ++jjj )
            {
                size_t idx1 = randomx()%nElts;
                size_t idx2 = randomx()%nElts;
                using std::swap; swap(vec[idx1],vec[idx2]);
            }
            sortInPlace(vec.begin(),vec.end());
            testOrder(vec.begin(),vec.end());
        }
    }

    void test_SortInPlaceLotsOfEqualElementsP()
    {
        size_t const nElts = 5000000;
        size_t const maxInstances = 500;
        ULongVec vec;
        vec.reserve(nElts);
        while ( vec.size() < vec.capacity() )
        {
            size_t val = randomx()%nElts;
            size_t nnn = randomx()%maxInstances + 1;
            while ( vec.size() < vec.capacity() && nnn-- )
                vec.push_back(val);
        }
        size_t const nPermutations = nElts;
        for ( size_t jjj = 0; jjj < nPermutations; ++jjj )
        {
            size_t idx1 = randomx()%nElts;
            size_t idx2 = randomx()%nElts;
            using std::swap; swap(vec[idx1],vec[idx2]);
        }
        sortInPlaceParallel(vec.begin(),vec.end());
        testOrder(vec.begin(),vec.end());
    }

    void test_SortInPlaceRandomUniqueP()
    {
        size_t const nElts = 10000;
        ULongVec vec;
        vec.reserve(nElts);
        for ( size_t iii = 0; iii < nElts; ++iii )
            vec.push_back(iii);

        size_t const nTests = 1000;
        for ( size_t iii = 0; iii < nTests; ++iii )
        {
            size_t const nPermutations = nElts;
            for ( size_t jjj = 0; jjj < nPermutations; ++jjj )
            {
                size_t idx1 = randomx()%nElts;
                size_t idx2 = randomx()%nElts;
                using std::swap; swap(vec[idx1],vec[idx2]);
            }
            sortInPlaceParallel(vec.begin(),vec.end());
            testOrder(vec.begin(),vec.end());
        }
    }
};


#endif /* SORTINPLACETEST_H_ */
