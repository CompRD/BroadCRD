///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PriorityQueueTests.h
 * \author tsharpe
 * \date Jul 22, 2011
 *
 * \brief
 */
#ifndef PRIORITYQUEUETESTS_H_
#define PRIORITYQUEUETESTS_H_

#include <cxxtest/TestSuite.h>
#include "feudal/PriorityQueue.h"

class PriorityQueueTests : public CxxTest::TestSuite
{
public:
    void test_stuff()
    {
        PriorityQueue<int> pq(5);
        TS_ASSERT_EQUALS(pq.empty(),true);
        pq.push(30);
        pq.push(40);
        pq.push(10);
        pq.push(0);
        pq.push(20);
        TS_ASSERT_EQUALS(pq.size(),5ul);
        TS_ASSERT_EQUALS(pq.full(),true);
        TS_ASSERT_EQUALS(pq.top(),0);
        pq.pop();
        TS_ASSERT_EQUALS(pq.size(),4ul);
        TS_ASSERT_EQUALS(pq.top(),10);
        pq.replaceTop(35);
        TS_ASSERT_EQUALS(pq.top(),20);
        pq.pop();
        TS_ASSERT_EQUALS(pq.top(),30);
        pq.pop();
        TS_ASSERT_EQUALS(pq.top(),35);
        pq.pop();
        TS_ASSERT_EQUALS(pq.top(),40);
        pq.pop();
        TS_ASSERT_EQUALS(pq.empty(),true);
    }
};

#endif /* PRIORITYQUEUETESTS_H_ */
