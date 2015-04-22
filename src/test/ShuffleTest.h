#ifndef _SHUFFLE_TEST_H
#define _SHUFFLE_TEST_H

#include <cxxtest/TestSuite.h>
#include "random/Shuffle.h"
#include "Vec.h"

// this tests that shuffle produces invariants based on a seed.
//
// if this would change, all hell would break loose.


typedef int integral_type;
typedef uint64_t integral64_type;

class ShuffleTest : public CxxTest::TestSuite
{
public:
    void verify(vec<integral_type> &res, integral_type *answerKey)
    {
        for(size_t idx = 0; idx < res.size(); ++idx)
        {
            TS_ASSERT_EQUALS(res[idx], answerKey[idx]);
        }
    }

    void run_shuffle(integral_type len, integral_type seed, integral_type *answerKey)
    {
        vec<integral_type> res;
        Shuffle(len, res, seed);
        /*
        for(size_t idx = 0; idx < res.size(); ++idx)
            printf("%d,", res[idx]);
        printf("\n");
        */
        verify(res, answerKey);
    }

    void verify64(vec<integral64_type> &res, integral64_type *answerKey)
    {
        for(size_t idx = 0; idx < res.size(); ++idx)
        {
            TS_ASSERT_EQUALS(res[idx], answerKey[idx]);
        }
    }


    void run_shuffle64(integral64_type len, integral64_type seed, integral64_type *answerKey)
    {
        vec<integral64_type> res;
        Shuffle64(len, res, seed);
        /*
        for(size_t idx = 0; idx < res.size(); ++idx)
            printf("%lu,", res[idx]);
        printf("\n");
        */
        verify64(res, answerKey);
    }

    void test_shuffle64_1()
    {
        integral64_type answerKey[10] = {3,8,1,4,0,2,5,7,9,6};
        run_shuffle64(10UL, 101UL, answerKey);
    }

    void test_shuffle64_2()
    {
        integral64_type answerKey[10] = {9,3,5,4,0,2,6,8,7,1};
        run_shuffle64(10UL, 28286UL, answerKey);
    }

    void test_shuffle64_3()
    {
        integral64_type answerKey[10] = {4,0,5,3,1,7,2,8,9,6};
        run_shuffle64(10UL, 8268UL, answerKey);
    }

    void test_shuffle1()
    {
        integral_type answerKey[10] = {2,3,7,5,0,8,6,4,1,9};
        run_shuffle(10, 101, answerKey);
    }

    void test_shuffle2()
    {
        integral_type answerKey[10] = {6,3,8,0,5,2,1,7,4,9};
        run_shuffle(10, 28286, answerKey);
    }

    void test_shuffle3()
    {
        integral_type answerKey[10] = {1,8,0,9,5,6,3,7,4,2};
        run_shuffle(10, 8262, answerKey);
    }
};

#endif // __MEMPOOLTESTS_H
