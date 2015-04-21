#ifndef __SMALLVECTESTS_H
#define __SMALLVECTESTS_H

#include <cxxtest/TestSuite.h>
#include "STLExtensions.h"
#include "feudal/Mempool.h"
#include "feudal/SmallVecDefs.h"

typedef unsigned int value_type;
typedef SmallVec<value_type, MempoolAllocator<value_type> > SVec;
template class SmallVec< value_type, MempoolAllocator<value_type> >;

#define ELEMS 0xFFu

class SmallVecTests : public CxxTest::TestSuite
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
    SVec* build_vec()
    {
        SVec* vec = new SVec(mPool);
        vec->resize(ELEMS); 
        for(unsigned int x = 0; x < ELEMS; x++)
        {
            vec->at(x) = x; 
        }
        return vec;
    }

    SVec _build_vec(value_type* array, size_t len)
    { //SVec vec(MempoolAllocator<unsigned int>(poolID));
      SVec vec;
      for(size_t x = 0; x < len; x++)
        vec.push_back(array[x]);
      return vec; }

    void _verify_vec(SVec& target, value_type* array, size_t len)
    { TS_ASSERT_EQUALS(target.size(), len);
      for(size_t x = 0; x < len; x++)
        TS_ASSERT_EQUALS(target[x], array[x]); }
    
    void verify_vec(SVec *vec)
    {
        TS_ASSERT(!vec->empty()); 
        TS_ASSERT_EQUALS(vec->size(), ELEMS);
        for(unsigned int x = 0; x < ELEMS; x++)
        {
            TS_ASSERT_EQUALS(vec->at(x), x); 
            TS_ASSERT_EQUALS((*vec)[x], x); 
        }
    }

    // test assertions of an empty vec
    void test_SmallVec_empty()
    { 
        SVec vec;
        TS_ASSERT(vec.empty()); 
        TS_ASSERT_EQUALS(vec.size(), 0u); 
    }

    void test_SmallVec_Insert()
    { 
        { value_type source[] = {9,8,7,6};
          value_type target[] = {9,8,1,7,6};
          SVec vec = _build_vec(source, 4);
          vec.insert(vec.begin(2), 1);
          _verify_vec(vec, target, 5); }
        { value_type source[] = {9,8,7,6,5,4,3,2};
          value_type target[] = {9,8,1,1,1,7,6,5,4,3,2};
          SVec vec = _build_vec(source, 8);
          // the "3u" differs the types.  otherwise, this call
          // will match the iterator pairs.
          vec.insert(vec.begin(2), 3u, 1);
          _verify_vec(vec, target, 11); }
        { value_type source[] = {9,8,7,6,5,4,3,2};
          value_type iarray[] = {1,2,3};
          value_type target[] = {9,8,1,2,3,7,6,5,4,3,2};
          SVec vec = _build_vec(source, 8);
          SVec ivec = _build_vec(iarray, 3);
          vec.insert(vec.begin(2), ivec.begin(), ivec.end());
          _verify_vec(vec, target, 11); }
    }

    // test access of vec 
    void test_SmallVec_access()
    {
        verify_vec(vec3);
    } 

    // test the copy operator
    void test_SmallVec_copy()
    { 
        SVec target_vec(mPool);
        // start a copy
        target_vec = *(vec3);
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 

        for(unsigned int x = 0; x < ELEMS; x++)
        {
            TS_ASSERT_EQUALS(target_vec.at(x), x); 
            TS_ASSERT_EQUALS(target_vec[x], x); 
        }
    } 

    // test the size and capacity methods
    void test_SmallVec_capacity()
    { 
        SVec vec(mPool);

        const SVec::size_type old_size     = vec.size();
        const SVec::size_type old_capacity = vec.capacity();
        const SVec::size_type new_capacity = old_capacity + 10;

        // Assert max_size returns a positive value
        TS_ASSERT(vec.max_size());

        vec.reserve(new_capacity);

        // we should of allocated at least new_capacity
        TS_ASSERT(new_capacity <= vec.capacity());
        // the size should not have changed
        TS_ASSERT(old_size == vec.size());
    } 

    void test_SmallVec_construct_sized()
    { 
        SVec target_vec = SVec(ELEMS, 1u);
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1u);
        }
    } 

    void test_SmallVec_construct_iter()
    { 
        SVec target_vec = SVec(vec3->begin(), vec3->end());
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_SmallVec_construct_copy()
    { 
        SVec target_vec = SVec(*vec3);
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_SmallVec_assign()
    { 
        SVec& target_vec = *(vec3);
        target_vec.assign(ELEMS, 1u);

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        for(unsigned int x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1u);
        }
    } 

    void test_SmallVec_assign_iters()
    { 
        SVec target_vec;
        target_vec.assign(vec3->begin(), vec3->end());
        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS); 
        verify_vec(&target_vec);
    } 

    void test_SmallVec_erase()
    { 
        SVec& target_vec = *(vec3);
        SVec::iterator itr = target_vec.begin();
        target_vec.erase(itr);

        TS_ASSERT_EQUALS(target_vec.size(), ELEMS - 1); 
        TS_ASSERT_EQUALS(target_vec.front(), 1u); 
        TS_ASSERT_EQUALS(target_vec.back(), ELEMS - 1); 
    } 

    void test_SmallVec_erase_iters()
    { 
        SVec& target_vec = *(vec3);
        SVec::iterator beg = target_vec.begin() + 1;
        SVec::iterator end = target_vec.end() - 1;
        target_vec.erase(beg, end);

        TS_ASSERT_EQUALS(target_vec.size(), 2u); 
        TS_ASSERT_EQUALS(target_vec.front(), 0u); 
        TS_ASSERT_EQUALS(target_vec.back(), ELEMS - 1); 
    } 

    void test_SmallVec_pop_back()
    { 
        SVec& target_vec = *(vec3);
        target_vec.pop_back();

        TS_ASSERT_EQUALS(target_vec.size(), ELEMS - 1); 
        TS_ASSERT_EQUALS(target_vec.back(), ELEMS - 2); 
    } 

    void test_SmallVec_insert()
    { 
        SVec& target_vec = *(vec3);
        target_vec.append(ELEMS, 1);

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS * 2);
        unsigned x = 0;
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], x);
        }
        for(; x < (ELEMS * 2); ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1u);
        }
    } 

    void test_SmallVec_append()
    { 
        SVec& target_vec = *(vec3);
        target_vec.append(ELEMS, 1);

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS * 2);
        unsigned x = 0;
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], x);
        }
        for(; x < (ELEMS * 2); ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1u);
        }
    } 

    void test_SmallVec_append_iters()
    { 
        SVec append_vec(ELEMS, 1);
        SVec& target_vec = *(vec3);
        target_vec.append(append_vec.begin(), append_vec.end());

        TS_ASSERT(!target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS * 2);
        unsigned x = 0;
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], x);
        }
        for(; x < (ELEMS * 2); ++x)
        {
            TS_ASSERT_EQUALS(target_vec[x], 1u);
        }
    } 
 
    void test_SmallVec_eraseIf()
    { 
        //SVec& target_vec = *(vec3);
        SVec target_vec;
        SVec evens_vec;

        for(unsigned int x = 0; x < ELEMS; x++)
        {
            target_vec.push_back(x);
            evens_vec.push_back(x & 1);
        }

        target_vec.eraseIf(evens_vec.begin(), evens_vec.end());
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS / 2 + ELEMS % 2);
        for(unsigned int x = 0; x < target_vec.size(); x++)
        {
            TS_ASSERT_EQUALS(target_vec[x], (x * 2));
        }
    } 

    void test_SmallVec_eraseUnless()
    { 
        //SVec& target_vec = *(vec3);
        SVec target_vec;
        SVec evens_vec;

        for(unsigned int x = 0; x < ELEMS; x++)
        {
            target_vec.push_back(x);
            evens_vec.push_back(!(x & 1));
        }

        target_vec.eraseUnless(evens_vec.begin(), evens_vec.end());
        TS_ASSERT_EQUALS(target_vec.size(), ELEMS / 2 + ELEMS % 2);
        for(unsigned int x = 0; x < target_vec.size(); x++)
        {
            TS_ASSERT_EQUALS(target_vec[x], (x * 2));
        }
    } 
   
    // the the ability to swap to vectors
    void test_SmallVec_swap()
    { 
        SVec target_vec(mPool);
        SVec& source_vec = *(vec3);

        /* target vector should be empty */
        TS_ASSERT(target_vec.empty()); 
        TS_ASSERT_EQUALS(target_vec.size(), 0u); 

        source_vec.swap(target_vec);

        /* target vector should now have content */
        verify_vec(&target_vec);

        /* source vector should now be empty */
        TS_ASSERT(source_vec.empty()); 
        TS_ASSERT_EQUALS(source_vec.size(), 0u); 
    } 

    // test the ability to swap with a vector
    // in a different memory pool
    void test_SmallVec_swap2()
    {
        MempoolOwner<value_type> pool;
        SVec* vPtr = new SVec(pool);
        SVec& target_vec = *(vPtr);
        SVec& source_vec = *(vec3);

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

    // test the operators of SmallVec
    void test_SmallVec_operators()
    { 
        SVec smaller_vec(mPool);
        SVec& bigger_vec = *(vec3);

        TS_ASSERT(bigger_vec != smaller_vec);
        TS_ASSERT(!(bigger_vec == smaller_vec));
        TS_ASSERT(bigger_vec > smaller_vec);
        TS_ASSERT(bigger_vec >= smaller_vec);
        TS_ASSERT(smaller_vec < bigger_vec);
        TS_ASSERT(smaller_vec <= bigger_vec);
    } 
    
    // test iterator access
    void test_SmallVec_iterators()
    { 
        SVec& vec = *(vec3);
        SVec::iterator itr;
        SVec::reverse_iterator ritr;

        // test forward access via dereference
        unsigned int x = 0;
        for(itr = vec.begin(); itr != vec.end(); ++itr)
        {
            TS_ASSERT(*itr == x++);
        }

        // test forward access via index
        itr = vec.begin();
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT(itr[x] == x);
        }

        // test reverse access via dereference
        x = ELEMS - 1;
        for(ritr = vec.rbegin(); ritr != vec.rend(); ++ritr)
        {
            TS_ASSERT_EQUALS(*ritr,  x--);
        }

        // test reverse access via index
        ritr = vec.rbegin();
        for(x = 0; x < ELEMS; ++x)
        {
            TS_ASSERT_EQUALS(ritr[x],  ELEMS - x -1);
        }
    }

    void test_push_front()
    {
        SVec vec{1u,2u,3u};
        vec.push_front(4u);
        SVec expect{4u,1u,2u,3u};
        TS_ASSERT_EQUALS(vec,expect);
    }

    void test_eraseValue()
    {
        SVec vec{1u,2u,0u,2u,3u,1u,2u};
        EraseValue(vec,1u);
        SVec expect1{2u,0u,2u,3u,2u};
        TS_ASSERT_EQUALS(vec,expect1);
        EraseValue(vec,0u);
        SVec expect2{2u,2u,3u,2u};
        TS_ASSERT_EQUALS(vec,expect2);
        EraseValue(vec,3u);
        SVec expect3{2u,2u,2u};
        TS_ASSERT_EQUALS(vec,expect3);
        EraseValue(vec,2u);
        SVec expect4;
        TS_ASSERT_EQUALS(vec,expect4);
    }

private:
    MempoolOwner<value_type> mPool;
    SVec* vec3;
};

#endif // __SMALLVECTESTS_H
