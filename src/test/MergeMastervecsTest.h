////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file MergeMastervecsTest.h
 * \author tsharpe
 * \date Nov 4, 2009
 *
 * \brief
 */
#ifndef MERGEMASTERVECSTEST_H_
#define MERGEMASTERVECSTEST_H_

#include <cxxtest/TestSuite.h>
#include "Basevector.h"

class MergeMastervecsTest : public CxxTest::TestSuite
{
public:
    void test_merge()
    {
        bvec bv1("ACGTACGTACG");
        bvec bv2("TGCATGC");
        vecbvec v1;
        v1.push_back(bv1).push_back(bv2);
        vecbvec v2;
        v1.push_back(bv2).push_back(bv1).push_back(bv2);
        v1.WriteAll("/tmp/file1");
        v2.WriteAll("/tmp/file2");
        vecbvec v3;
        v3.append(v1.begin(),v1.end());
        v3.append(v2.begin(),v2.end());
        v3.append(v1.begin(),v1.end());
        std::vector<String> filenames;
        filenames.push_back("/tmp/file1");
        filenames.push_back("/tmp/file2");
        filenames.push_back("/tmp/file1");
        MergeMastervecs("/tmp/file3",filenames);
        vecbvec v4("/tmp/file3");
        remove("/tmp/file1");
        remove("/tmp/file2");
        remove("/tmp/file3");
        TS_ASSERT_EQUALS(v3,v4);
    }
};

#endif /* MERGEMASTERVECSTEST_H_ */
