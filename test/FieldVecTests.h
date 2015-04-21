#ifndef __FIELDVECTESTS_H
#define __FIELDVECTESTS_H

#include <cxxtest/TestSuite.h>
#include "feudal/Mempool.h"
#include "feudal/FieldVecDefs.h"

typedef FieldVec<1, MempoolAllocator<unsigned char> > bitvec;
typedef FieldVec<2, MempoolAllocator<unsigned char> > basevec;
typedef FieldVec<4, MempoolAllocator<unsigned char> > nibblevec;
typedef FieldVec<4, std::allocator<unsigned char> > nibblevec_stdAlloc;

#define ELEMS 0xFFu

class FieldVecTests : public CxxTest::TestSuite
{
public:
    void setUp()
    { 
        _vec1 = build_vec<1>();
        _vec2 = build_vec<2>();
        _vec4 = build_vec<4>();
    }

    void tearDown()
    {
        delete _vec1;
        delete _vec2;
        delete _vec4;
    }

    template <int N, class A1>
    void verify_vec(FieldVec<N, A1> *vecPtr)
    { 
        TS_ASSERT_EQUALS(vecPtr->size(), ELEMS);
        unsigned char maxval = static_cast<unsigned char>(pow(2, N));

        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(vecPtr->at(x), (x % maxval));
        }
    }

    template <int N>
    FieldVec<N, MempoolAllocator<unsigned char> >* build_vec()
    {
        typedef FieldVec<N, MempoolAllocator<unsigned char> > FV;
        unsigned char maxval = static_cast<unsigned char>(pow(2, N));
        FV* ret = new FV;

        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            ret->push_back(x % maxval);
        }
        return ret;
    }

    void test_FieldVec_empty()
    { 
        bitvec vec;
        TS_ASSERT(vec.empty()); 
    }

    void test_FieldVec_access()
    { 
        verify_vec(_vec1);
        verify_vec(_vec2);
        verify_vec(_vec4);
    }
    
    // test the size and capacity methods
    void test_FieldVec_capacity()
    { 
        bitvec vec;

        const bitvec::size_type old_size     = vec.size();
        const bitvec::size_type old_capacity = vec.capacity();
        const bitvec::size_type new_capacity = old_capacity + 10;

        // Assert max_size returns a positive value
        TS_ASSERT(vec.max_size());

        vec.reserve(new_capacity);

        // we should of allocated at least new_capacity
        TS_ASSERT(new_capacity <= vec.capacity());
        // the size should not have changed
        TS_ASSERT(old_size == vec.size());
    } 

    void test_FieldVec_copy()
    { 
        nibblevec target_vec;
        target_vec = *(_vec4);
        verify_vec(&target_vec);
    } 

    void test_FieldVec_construct_sized()
    { 
        nibblevec target_vec = nibblevec(ELEMS, 1);
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1);
        }
    } 

    void test_FieldVec_construct_iter()
    { 
        nibblevec target_vec = nibblevec(_vec4->begin(), _vec4->end());
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_FieldVec_construct_copy()
    { 
        nibblevec target_vec = nibblevec(*_vec4);
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 
    
    // test the operators of FieldVec
    void test_FieldVec_operators()
    { 
        nibblevec that_vec;
        nibblevec& this_vec = *(_vec4);

        that_vec = *(build_vec<4>());
        TS_ASSERT(this_vec == that_vec);
        TS_ASSERT(that_vec <= this_vec);
        TS_ASSERT(this_vec >= that_vec);

        // this vec is now one element larger
        that_vec.pop_back();

        TS_ASSERT(that_vec < this_vec);
        TS_ASSERT(that_vec <= this_vec);
        TS_ASSERT(this_vec >= that_vec);
        TS_ASSERT(this_vec > that_vec);
        TS_ASSERT(this_vec != that_vec);
        TS_ASSERT(!(this_vec == that_vec));
    } 

    // test iterator access
    void test_FieldVec_iterators()
    { 
        nibblevec& vec = *(_vec4);
        nibblevec::iterator itr;
        nibblevec::reverse_iterator ritr;
        nibblevec::lvalue_iterator litr;

        // test forward access 
        int x = 0;
        for(itr = vec.begin(); itr != vec.end(); ++itr)
        {
            TS_ASSERT_EQUALS(*itr, (x++ % 16));
        }

        // test reverse access 
        x = ELEMS - 1;
        for(ritr = vec.rbegin(); ritr != vec.rend(); ++ritr)
        {
            TS_ASSERT_EQUALS(*ritr, (x-- % 16));
        }

        litr = vec.lbegin();
        *litr = 3;
        itr = vec.begin();
        TS_ASSERT_EQUALS(*litr, *itr);
        TS_ASSERT_EQUALS(3, *itr);
        TS_ASSERT_EQUALS(3, vec[0]);
    }

    void test_FieldVec_assign()
    { 
        nibblevec& target_vec = *(_vec4);
        target_vec.assign(ELEMS, 1);

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1);
        }
    } 

    void test_FieldVec_assign_iters()
    { 
        nibblevec target_vec;
        target_vec.assign(_vec4->begin(), _vec4->end());
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_FieldVec_erase()
    { 
        nibblevec& target_vec = *(_vec4);
        nibblevec::iterator itr = target_vec.begin();
        target_vec.erase(itr);

        TS_ASSERT_EQUALS(target_vec.size(), ELEMS - 1); 
        TS_ASSERT_EQUALS(target_vec.front(), 1); 
        TS_ASSERT_EQUALS(target_vec.back(), 14); 
    } 

    void test_FieldVec_erase_iters()
    { 
        nibblevec& target_vec = *(_vec4);
        nibblevec::iterator beg = target_vec.begin() + 1;
        nibblevec::iterator end = target_vec.end() - 1;
        target_vec.erase(beg, end);

        TS_ASSERT_EQUALS(target_vec.size(), 2u); 
        TS_ASSERT_EQUALS(target_vec.front(), 0); 
        TS_ASSERT_EQUALS(target_vec.back(), 14); 
    } 

    void test_FieldVec_pop_back()
    { 
        nibblevec& target_vec = *(_vec4);
        target_vec.pop_back();

        TS_ASSERT_EQUALS(target_vec.size(), ELEMS - 1); 
        TS_ASSERT_EQUALS(target_vec.back(), 13); 
    } 

    void test_FieldVec_append()
    { 
        nibblevec& target_vec = *(_vec4);
        target_vec.append(ELEMS, 1);

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS * 2);
        unsigned x = 0;
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], (x % 16));
        }
        for(; x < (ELEMS * 2); ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1);
        }
    } 

    void test_FieldVec_append_iters()
    { 
        nibblevec append_vec(ELEMS, 1);
        nibblevec& target_vec = *(_vec4);
        target_vec.append(append_vec.begin(), append_vec.end());

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS * 2);
        unsigned x = 0;
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], (x % 16));
        }
        for(; x < (ELEMS * 2); ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1);
        }
    } 
   
   
    void test_FieldVec_swap()
    { 
        nibblevec* vPtr = new nibblevec;
        nibblevec& source_vec = *(_vec4);
        nibblevec& target_vec = *(vPtr);

        /* target vector should be empty */
        TS_ASSERT(target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), 0u); 

        source_vec.swap(target_vec);

        /* target vector should now have content */
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], (x % 16));
        }

        /* source vector should now be empty */
        TS_ASSERT(source_vec.empty()); 
        TS_ASSERT_EQUALS(source_vec.size(), 0u); 

        delete vPtr;
    } 
    
    // test the ability to swap with a vector
    // in a different memory pool
    void test_FieldVec_swap2()
    {
        MempoolOwner<unsigned char> owner;
        nibblevec* vPtr = new nibblevec(owner);
        nibblevec& source_vec = *(_vec4);
        nibblevec& target_vec = *(vPtr);

        /* target vector should be empty */
        TS_ASSERT(target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), 0u); 

        source_vec.swap(target_vec);

        /* target vector should now have content */
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], (x % 16));
        }

        /* source vector should now be empty */
        TS_ASSERT(source_vec.empty()); 
        TS_ASSERT_EQUALS(source_vec.size(), 0u); 

        delete vPtr;
    } 

    bitvec* _vec1;
    basevec* _vec2;
    nibblevec* _vec4;
};

#endif // __FIELDVECTESTS_H
