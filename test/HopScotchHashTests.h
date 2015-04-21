#ifndef __BASEVEC_TESTS_H
#define __BASEVEC_TESTS_H

#include "feudal/HashSet.h"
#include <cxxtest/TestSuite.h>
#include <algorithm>
#include <vector>

struct NopHash : public std::unary_function<size_t,size_t>
{ size_t operator()( size_t key ) const { return key; } };

typedef HopscotchHashSet<size_t, NopHash> HS;
typedef std::vector<size_t> Vec;

class HopScotchHashTests : public CxxTest::TestSuite
{
public:
    void test_basic_storage()
    {
        HCF<size_t,NopHash> hcf;
        HS hht(0x200,hcf);
        Vec clist;
        clist.reserve(0x100);
        for ( size_t elem = 0; elem < 0x100; ++elem )
        {
            size_t key = rand();
            while ( !hht.add(key) )
            {
                TS_ASSERT(std::find(clist.begin(),clist.end(),key) != clist.end());
                key = rand();
            }
            clist.push_back(key);
        }

        TS_ASSERT_EQUALS(hht.size(),clist.size());

        for ( Vec::iterator itr = clist.begin(); itr != clist.end(); ++itr )
        {
            TS_ASSERT(hht.contains(*itr));
            TS_ASSERT(hht.remove(*itr));
            TS_ASSERT(!hht.contains(*itr));
        }

        TS_ASSERT_EQUALS(hht.size(),0u);
    }
};

#endif
