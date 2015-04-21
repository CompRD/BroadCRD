#ifndef __OUTERVECTESTS_H
#define __OUTERVECTESTS_H

#include <cxxtest/TestSuite.h>
#include "feudal/Mempool.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

typedef SmallVec< int, MempoolAllocator<int> > SVec;
template class SmallVec< int, MempoolAllocator<int> >;

typedef OuterVec<SVec> OVec;
template class OuterVec<SVec>;
typedef OVec::size_type idx_type;

#define ELEMS 0xFFu

class OuterVecTests : public CxxTest::TestSuite
{
public:
    void setUp()
    {
        // vec3 is a vector, with 1, 2, 3 as its items
        vec3 = build_vec();
    }

    void tearDown()
    {
        delete vec3;
    }

    // Construct a vector of N length for testing purposes
    OVec* build_vec()
    {
        OVec* vec = new OVec();
        vec->resize(ELEMS);
        for(SVec::size_type idx = 0; idx < ELEMS; ++idx)
        {
            (*vec)[idx] = SVec(idx, static_cast<SVec::value_type>(idx));
        }
        return vec;
    }


    void verify_vec(const OVec *vec_ptr)
    {
        const OVec &vec = *vec_ptr;
        TS_ASSERT(!vec.empty());
        TS_ASSERT_EQUALS(vec.size(), ELEMS);

        for(SVec::value_type o_idx = 0; 
            o_idx < static_cast<SVec::value_type>(ELEMS); ++o_idx)
        {
            TS_ASSERT_EQUALS(vec[o_idx].size(), 
                                static_cast<SVec::size_type>(o_idx));
            for(SVec::value_type i_idx = 0; i_idx < o_idx; ++i_idx)
            {
                TS_ASSERT_EQUALS(vec[o_idx][i_idx], o_idx);
                TS_ASSERT_EQUALS(vec.at(o_idx).at(i_idx), o_idx);
            }
        }
    }


    // test assertions of an empty vec
    void test_OuterVec_empty()
    {
        OVec vec;
        TS_ASSERT(vec.empty());
        TS_ASSERT_EQUALS(vec.size(), 0u);
    }

    // test the size and capacity methods
    void test_OuterVec_capacity()
    {
        OVec vec;

        const OVec::size_type old_size     = vec.size();
        const OVec::size_type old_capacity = vec.capacity();
        const OVec::size_type new_capacity = old_capacity + 10;

        // Assert max_size returns a positive value
        TS_ASSERT(vec.max_size());

        vec.reserve(new_capacity);

        // we should of allocated at least new_capacity
        TS_ASSERT(new_capacity <= vec.capacity());
        // the size should not have changed
        TS_ASSERT(old_size == vec.size());
    }

    // test access of vec
    void test_OuterVec_access()
    {
        verify_vec(vec3);
    }

    // test the copy operator
    void test_OuterVec_copy()
    {
        OVec vec = *(vec3);
        verify_vec(&vec);
    }

    // the the ability to swap to vectors
    void test_OuterVec_swap()
    {
        OVec target_vec;
        OVec& source_vec = *(vec3);
        source_vec.swap(target_vec);
        verify_vec(&target_vec);
        TS_ASSERT_EQUALS(source_vec.size(), 0U);
    }

    // test the ability to swap with a vector
    // in a different memory pool
    void test_OuterVec_swap2()
    { 
        OVec* vPtr = new OVec;
        OVec& target_vec = *(vPtr);
        OVec& source_vec = *(vec3);

        /* target vector should be empty */
        TS_ASSERT(target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), 0u); 

        source_vec.swap(target_vec);

        /* target vector should now have content */
        verify_vec(&target_vec);

        /* source vector should now be empty */
        TS_ASSERT(source_vec.empty()); 
        TS_ASSERT_EQUALS(source_vec.size(), 0u); 

        delete vPtr;
    } 

    void test_OuterVec_swap3()
    {
        OVec a, b;

        a.push_back(SVec(0U,1));
        a.push_back(SVec(1U,2));
        a.push_back(SVec(2U,3));

        b.push_back(SVec(3U,4));
        b.push_back(SVec(4U,5));

        a.swap(b);
        TS_ASSERT( a.size()==2 && b.size() == 3 );
    }

    void test_OuterVec_construct_sized()
    { 
        OVec target_vec = OVec(ELEMS, SVec(ELEMS, 1));
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(ELEMS, 1));
        }
    } 

    void test_OuterVec_construct_iter()
    { 
        OVec target_vec = OVec(vec3->begin(), vec3->end());
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_OuterVec_construct_copy()
    { 
        OVec target_vec = OVec(*vec3);
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_OuterVec_assign()
    { 
        OVec& target_vec = *(vec3);
        target_vec.assign(ELEMS, SVec(ELEMS, 1));

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(ELEMS, 1));
        }
    } 

    void test_OuterVec_assign_iters()
    { 
        OVec target_vec;
        target_vec.assign(vec3->begin(), vec3->end());
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_OuterVec_erase()
    { 
        OVec& target_vec = *(vec3);
        OVec::iterator itr = target_vec.begin();
        target_vec.erase(itr);

        TS_ASSERT_EQUALS(target_vec.size(), ELEMS - 1); 
        TS_ASSERT_EQUALS(target_vec.front(), SVec(1u, 1));
        TS_ASSERT_EQUALS(target_vec.back(), SVec(ELEMS - 1, static_cast<SVec::value_type>(ELEMS - 1)));
    } 

    void test_OuterVec_erase_iters()
    { 
        OVec& target_vec = *(vec3);
        OVec::iterator beg = target_vec.begin() + 1;
        OVec::iterator end = target_vec.end() - 1;
        target_vec.erase(beg, end);

        TS_ASSERT_EQUALS(target_vec.size(), 2u); 
        TS_ASSERT_EQUALS(target_vec.front(), SVec());
        TS_ASSERT_EQUALS(target_vec.back(), SVec(ELEMS - 1, static_cast<SVec::value_type>(ELEMS - 1)));
    } 

    void test_OuterVec_pop_back()
    { 
        OVec& target_vec = *(vec3);
        target_vec.pop_back();

        TS_ASSERT_EQUALS(target_vec.size(), ELEMS - 1); 
        TS_ASSERT_EQUALS(target_vec.back(), SVec(ELEMS - 2, static_cast<SVec::value_type>(ELEMS - 2)));
    } 

    void test_OuterVec_append()
    { 
        OVec& target_vec = *(vec3);
        target_vec.append(ELEMS, SVec(ELEMS, 1));

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS * 2);
        unsigned x = 0;
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(static_cast<SVec::size_type>(x), static_cast<SVec::value_type>(x)));
        }
        for(; x < (ELEMS * 2); ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(ELEMS, 1));
        }
    } 

    void test_OuterVec_append_iters()
    { 
        OVec *append_vec = build_vec();
        OVec& target_vec = *(vec3);
        target_vec.append(append_vec->begin(), append_vec->end());

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS * 2);
        unsigned x = 0;
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(static_cast<SVec::size_type>(x), static_cast<SVec::value_type>(x)));
        }
        for(; x < (ELEMS * 2); ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(static_cast<SVec::size_type>(x - ELEMS), static_cast<SVec::value_type>(x - ELEMS)));
        }
        delete append_vec;
    } 

    void test_OuterVec_eraseIf()
    { 
        OVec target_vec;
        SVec evens_vec;

        for(unsigned int x = 0; x < ELEMS; x++)
        {
            target_vec.push_back(SVec(x, static_cast<SVec::value_type>(x)));
            evens_vec.push_back(x & 1);
        }

        target_vec.eraseIf(evens_vec.begin(), evens_vec.end());
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS / 2 + ELEMS % 2);
        for(unsigned int x = 0; x < target_vec.size(); x++)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(x*2, static_cast<SVec::value_type>(x*2)));
        }
    } 

    void test_OuterVec_eraseUnless()
    { 
        OVec target_vec;
        SVec evens_vec;

        for(unsigned int x = 0; x < ELEMS; x++)
        {
            target_vec.push_back(SVec(x, static_cast<SVec::value_type>(x)));
            evens_vec.push_back(!(x & 1));
        }

        target_vec.eraseUnless(evens_vec.begin(), evens_vec.end());
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS / 2 + ELEMS % 2);
        for(unsigned int x = 0; x < target_vec.size(); x++)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(x*2, static_cast<SVec::value_type>(x*2)));
        }
    } 
   
    void test_OuterVec_eraseEntries()
    { 
        OVec target_vec;
        SVec evens_vec;

        for(unsigned int x = 0; x < ELEMS; x++)
        {
            target_vec.push_back(SVec(x, static_cast<SVec::value_type>(x)));
            if(x & 1)
            {
                evens_vec.push_back(x);
            }
        }

        target_vec.eraseEntries(evens_vec.begin(), evens_vec.end());
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS - evens_vec.size());
        for(unsigned int x = 0; x < target_vec.size(); x++)
        {
            TS_ASSERT_EQUALS(target_vec[x], SVec(x*2, static_cast<SVec::value_type>(x*2)));
        }
    } 
   

    // test the operators
    void test_OuterVec_operators()
    {
        OVec smaller_vec;
        OVec& bigger_vec = *(vec3);

        TS_ASSERT(bigger_vec != smaller_vec);
        TS_ASSERT(!(bigger_vec == smaller_vec));
        TS_ASSERT(bigger_vec > smaller_vec);
        TS_ASSERT(bigger_vec >= smaller_vec);
        TS_ASSERT(smaller_vec < bigger_vec);
        TS_ASSERT(smaller_vec <= bigger_vec);
    }

    // test iterator access
    void test_OuterVec_iterators()
    {
        OVec& vec = *(vec3);
        OVec::iterator itr;
        OVec::reverse_iterator ritr;

        // test forward access via dereference
        OVec::size_type x = 0;
        for(itr = vec.begin(); itr != vec.end(); ++itr)
        {
            TS_ASSERT_EQUALS(*itr, SVec(x, static_cast<SVec::value_type>(x)));
            ++x;
        }

        // test forward access via index
        itr = vec.begin();
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(itr[x], SVec(x, static_cast<SVec::value_type>(x)));
        }

        // test reverse access via dereference
        x = ELEMS - 1;
        for(ritr = vec.rbegin(); ritr != vec.rend(); ++ritr)
        {
            TS_ASSERT_EQUALS(*ritr, SVec(x, static_cast<SVec::value_type>(x)));
            --x;
        }

        // test reverse access via index
        ritr = vec.rbegin();
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(ritr[x], SVec(ELEMS - 1 - x,  static_cast<SVec::value_type>(ELEMS - 1 - x)));
        }
    }

    // test preallocation
    void test_OuterVec_preallocation()
    {
        SVec sv(1024U,47);
        OVec ov1(1024U,sv);
        OVec ov2(ov1);
        TS_ASSERT( ov1 == ov2 );
        ov2.assign(ov1.begin(),ov1.end());
        TS_ASSERT( ov1 == ov2 );
        ov2.append(ov1.begin(),ov1.end());
    }

private:
    OVec*               vec3;
};

#endif // __OUTERVECTESTS_H
