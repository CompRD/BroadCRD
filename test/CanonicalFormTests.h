///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * CanonicalFormTests.h
 *
 *  Created on: Nov 20, 2013
 *      Author: tsharpe
 */

#ifndef CANONICALFORMTESTS_H_
#define CANONICALFORMTESTS_H_

#include <cxxtest/TestSuite.h>
#include "Basevector.h"
#include "dna/CanonicalForm.h"

class CanonicalFormTests : public CxxTest::TestSuite
{
public:
    void testIsRC()
    {
        bvec bv1("AACCTT");
        bvec bv2("AAGGTT");
        unsigned const K = 6u;
        bool result = CF<K>::isRC(bv1.begin(),bv1.begin());
        TS_ASSERT(!result);
        result = CF<K>::isRC(bv1.begin(),bv2.begin());
        TS_ASSERT(result);
    }
};

#endif /* CANONICALFORMTESTS_H_ */
