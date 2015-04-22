#ifndef __MEMPOOLTESTS_H
#define __MEMPOOLTESTS_H

#include <cxxtest/TestSuite.h>
#include "feudal/Mempool.h"

class MempoolTests : public CxxTest::TestSuite
{
public:
    void test_allocatePool()
    {
        MempoolFinder& mpf = MempoolFinder::getInstance();
        unsigned short poolID = mpf.allocatePool();
        TS_ASSERT(poolID);
        mpf.freePool(poolID);
    }

    void test_freePool()
    {
        MempoolFinder& mpf = MempoolFinder::getInstance();
        unsigned short poolID = mpf.allocatePool();
        mpf.freePool(poolID);
        TS_ASSERT(mpf.resolvePool(poolID)->isUnused());
    }

    // not sure if this really tests anything
    void test_resolvePool()
    {
        MempoolFinder& mpf = MempoolFinder::getInstance();
        unsigned short poolID1 = mpf.allocatePool();
        unsigned short poolID2 = mpf.allocatePool();
        // should be the same ptr
        TS_ASSERT_EQUALS(mpf.resolvePool(poolID1), mpf.resolvePool(poolID1));
        // pointers should be different
        TS_ASSERT_DIFFERS(mpf.resolvePool(poolID1), mpf.resolvePool(poolID2));
        mpf.freePool(poolID1);
        mpf.freePool(poolID2);
    }
    
    void test_ref_deref()
    {
        MempoolFinder& mpf = MempoolFinder::getInstance();
        unsigned short poolID = mpf.allocatePool();
        Mempool* pPool = mpf.resolvePool(poolID);
        // since we allocated the pool, we have one reference
        TS_ASSERT_EQUALS(pPool->ref(), 1u);
        TS_ASSERT_EQUALS(pPool->deref(), 0u);
        mpf.freePool(poolID);
    }

    void test_allocate_and_free()
    {
        const int N = 1024;
        void* Ptrs[N];
        MempoolFinder& mpf = MempoolFinder::getInstance();
        unsigned short poolID = mpf.allocatePool();
        Mempool* pPool = mpf.resolvePool(poolID);

        for(int x = 0; x < N; x++)
        {
            int sz = 2 * x * x;
            void* p = pPool->allocate(sz,2);
            Ptrs[x] = p;
        }

        for(int x = 0; x < N; x++)
        {
            void* p = Ptrs[x];
            int sz = 2 * x * x;
            pPool->free(p, sz);
        }

        mpf.freePool(poolID);
    }

    void test_alignment()
    {
        MempoolFinder& mpf = MempoolFinder::getInstance();
        unsigned short poolID = mpf.allocatePool();
        Mempool* pPool = mpf.resolvePool(poolID);
        char* p1 = (char*)pPool->allocate(1,1);
        char* p2 = (char*)pPool->allocate(2,2);
        TS_ASSERT(p2-p1 == 2);
        char* p3 = (char*)pPool->allocate(1,1);
        TS_ASSERT(p3-p1 == 4);
        char* p4 = (char*)pPool->allocate(4,4);
        TS_ASSERT(p4-p1 == 8);
        char* p5 = (char*)pPool->allocate(8,8);
        TS_ASSERT(p5-p1 == 16);
        char* p6 = (char*)pPool->allocate(1,1);
        TS_ASSERT(p6-p1 == 24);
        char* p7 = (char*)pPool->allocate(7,1);
        TS_ASSERT(p7-p1 == 25);
        pPool->free(p7,7);
        pPool->free(p6,1);
        pPool->free(p5,8);
        pPool->free(p4,4);
        pPool->free(p3,1);
        pPool->free(p2,2);
        pPool->free(p1,1);
        mpf.freePool(poolID);
    }
};

#endif // __MEMPOOLTESTS_H
