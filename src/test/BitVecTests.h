#ifndef __BASEVEC_TESTS_H
#define __BASEVEC_TESTS_H

#include "feudal/BitVec.h"
#include "dna/Bases.h"
#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>

typedef BitVec::size_type BitVec_size_type;

class BitVecTests : public CxxTest::TestSuite
{
public:
    void setUp()
    { }

    void tearDown()
    { }

    BitVec build_vec(const String &s)
    {
        BitVec bv;
        for (String::const_iterator itr = s.begin(); itr != s.end(); ++itr)
        {
            bv.push_back((*itr == '1') ? 1 : 0);
        }
        return bv;
    }

    void verify_vec(const BitVec& bv, const String& s)
    {
        String::const_iterator str_itr = s.begin();

        TS_ASSERT_EQUALS(bv.size(), s.size());
        for(unsigned int i = 0; i < bv.size(); ++i)
        {
            TS_ASSERT_EQUALS(bv[i], ((*str_itr == '1') ? 1 : 0));
            ++str_itr;
        }
    }

    void test_BitVec_helpers()
    { 
        BitVec bv1 = build_vec("0101");
        verify_vec(bv1, "0101");
    }

    void test_BitVec_sum()
    { 
        BitVec bv1 = build_vec("0000101");
        TS_ASSERT_EQUALS(bv1.Sum(), 2u);
        bv1 = build_vec("100000101");
        bv1.resize(43,false);
        TS_ASSERT_EQUALS(bv1.Sum(), 3u);
        bv1.resize(57,true);
        TS_ASSERT_EQUALS(bv1.Sum(), 17u);
        bv1 = build_vec("00000000000");
        TS_ASSERT_EQUALS(bv1.Sum(), 0u);
    }

    void test_BitVec_invert()
    {
        BitVec bv1 = build_vec("101");
        bv1.invert();
        verify_vec(bv1, "010");
    }

    void test_BitVec_reverse()
    { 
        BitVec bv1 = build_vec("1010110010");
        bv1.ReverseMe();
        verify_vec(bv1, "0100110101");
    }

    void test_BitVec_nor()
    {
        BitVec bv1 = build_vec("1010110010");
        BitVec bv2 = build_vec("0110011001");
        BitVec bv3 = nor(bv1,bv2);
        verify_vec(bv3, "0001000100");
    }

    void test_BitVec_zero()
    { 
        BitVec bv1 = build_vec("111");
        bv1.Zero();
        verify_vec(bv1, "000");
    }

    void test_BitVec_WackySet()
    { 
        BitVec target = build_vec("111111");
        target.Set(2, 4, 0);
        verify_vec(target, "110011");
        target = build_vec("111111");
    }
};

#endif
