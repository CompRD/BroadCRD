///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FlatIteratorTest.h
 * \author tsharpe
 * \date Jan 13, 2010
 *
 * \brief
 */
#ifndef FLATITERATORTEST_H_
#define FLATITERATORTEST_H_

#include <cxxtest/TestSuite.h>
#include "Basevector.h"

class FlatIteratorTest : public CxxTest::TestSuite
{
public:
    void test_FlatIterator()
    {
        vecbvec vbv;
        vbv.push_back(bvec("A"));
        vbv.push_back(bvec("CC"));
        vbv.push_back(bvec("GGG"));

        String bases;
        for ( FlatIterator itr(vbv), end(vbv,vbv.size()); itr != end; ++itr )
            bases.push_back(Base::val2Char(*itr));
        TS_ASSERT_EQUALS( bases, "ACCGGG" );

        vbv.push_back(bvec());
        vbv.push_back(bvec("TTTT"));

        bases.clear();
        for ( FlatIterator itr(vbv), end(vbv,vbv.size()); itr != end; ++itr )
            bases.push_back(Base::val2Char(*itr));
        TS_ASSERT_EQUALS( bases, "ACCGGGTTTT" );

        vbv.push_back(bvec());

        bases.clear();
        for ( FlatIterator itr(vbv), end(vbv,vbv.size()); itr != end; ++itr )
            bases.push_back(Base::val2Char(*itr));
        TS_ASSERT_EQUALS( bases, "ACCGGGTTTT" );
    }
};
#endif /* FLATITERATORTEST_H_ */
