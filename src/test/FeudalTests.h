////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file FeudalTests.h
 * \author tsharpe
 * \date Sep 2, 2009
 *
 * \brief Tests feudal file reading and writing.
 */
#ifndef FEUDALTESTS_H_
#define FEUDALTESTS_H_

#include <cxxtest/TestSuite.h>
#include "Qualvector.h"
#include "feudal/VirtualMasterVec.h"
#include "feudal/FieldVec.h"

class FeudalTests : public CxxTest::TestSuite
{
public:
    void testWritingAndReadingQuals()
    {
        char const* filename = "/tmp/feudalTest.fastq";
        vecqvec vqv;
        vqv.reserve(OUTERVEC_SIZE);
        unsigned long nnn = OUTERVEC_SIZE;
        while ( nnn-- )
            vqv.push_back(qvec(static_cast<qvec::value_type>(nnn),nnn));
        vqv.WriteAll(filename);

        vecqvec vqv2(filename);
        TS_ASSERT_EQUALS(vqv,vqv2);

        remove(filename);
    }

    void testWritingAndReadingQualsVirtually()
    {
        char const* filename = "/tmp/feudalTest.fastq";
        vecqvec vqv;
        vqv.reserve(OUTERVEC_SIZE);
        unsigned long nnn = OUTERVEC_SIZE;
        while ( nnn-- )
            vqv.push_back(qvec(static_cast<qvec::value_type>(nnn),nnn));
        vqv.WriteAll(filename);

        VirtualMasterVec<qvec> vmv(filename);
        TS_ASSERT_EQUALS(vqv.size(),vmv.size());
        vecqvec::const_iterator end(vqv.cend());
        vecqvec::const_iterator itr(vqv.cbegin());
        VirtualMasterVec<qvec>::const_iterator itr2(vmv.begin());
        for ( ; itr != end; ++itr, ++itr2 )
            TS_ASSERT_EQUALS(*itr,*itr2);

        remove(filename);
    }

    typedef FieldVec<4,MempoolAllocator<unsigned char> > NybbleVec;
    typedef MasterVec<NybbleVec> VecNybbleVec;

    void testWritingAndReadingNybbles()
    {
        char const* filename = "/tmp/feudalTest.nybbles";
        VecNybbleVec vqv;
        vqv.reserve(OUTERVEC_SIZE);
        unsigned long nnn = OUTERVEC_SIZE;
        while ( nnn-- )
            vqv.push_back(NybbleVec(nnn,static_cast<NybbleVec::value_type>(nnn%16)));
        vqv.WriteAll(filename);

        VecNybbleVec vqv2(filename);
        TS_ASSERT_EQUALS(vqv,vqv2);

        remove(filename);
    }

    void testWritingAndReadingNybblesVirtually()
    {
        char const* filename = "/tmp/feudalTest.nybbles";
        VecNybbleVec vqv;
        vqv.reserve(OUTERVEC_SIZE);
        unsigned long nnn = OUTERVEC_SIZE;
        while ( nnn-- )
            vqv.push_back(NybbleVec(nnn,static_cast<NybbleVec::value_type>(nnn%16)));
        vqv.WriteAll(filename);

        VirtualMasterVec<NybbleVec> vmv(filename);
        TS_ASSERT_EQUALS(vqv.size(),vmv.size());
        VecNybbleVec::const_iterator end(vqv.cend());
        VecNybbleVec::const_iterator itr(vqv.cbegin());
        VirtualMasterVec<NybbleVec>::const_iterator itr2(vmv.begin());
        for ( ; itr != end; ++itr, ++itr2 )
            TS_ASSERT_EQUALS(*itr,*itr2);

        remove(filename);
    }

    static unsigned long const OUTERVEC_SIZE = 100;
};

#include "feudal/OuterVecDefs.h"
template class OuterVec<FieldVec<4,MempoolAllocator<unsigned char> > >;

#endif /* FEUDALTESTS_H_ */
