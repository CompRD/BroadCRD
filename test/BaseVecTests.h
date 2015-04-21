#ifndef __BASEVEC_TESTS_H
#define __BASEVEC_TESTS_H

#include <cxxtest/TestSuite.h>
#include "dna/Bases.h"
#include <iostream>
#include <fstream>

// #define TEST_BASEVECTOR 1

#ifndef TEST_BASEVECTOR
#include "feudal/BaseVec.h"
typedef BaseVec::size_type BaseVec_size_type;
#else
// #include "Basevector.h"
typedef basevector BaseVec;
typedef unsigned int BaseVec_size_type;
#endif // TEST_BASEVECTOR

class BaseVecTests : public CxxTest::TestSuite
{
public:
    void setUp()
    { }

    void tearDown()
    { }

    BaseVec build_vec(const String &s)
    {
        BaseVec bv;
        for (String::const_iterator itr = s.begin(); itr != s.end(); ++itr)
        {
            bv.AppendBase(as_char(*itr));
        }
        return bv;
    }

    void verify_vec(const BaseVec& bv, const String& s)
    {
        BaseVec::const_iterator bv_itr = bv.Begin();
        String::const_iterator str_itr = s.cbegin();

        TS_ASSERT_EQUALS(bv.size(), s.size());
        while(str_itr != s.end())
        {
            TS_ASSERT_EQUALS(*bv_itr, as_char(*str_itr));
            ++bv_itr;
            ++str_itr;
        }
    }

    void verify_cstr(const char* cstr1, const char* cstr2)
    {
        TS_ASSERT_EQUALS(strlen(cstr1), strlen(cstr2));
        while(*cstr1)
        {
            TS_ASSERT_EQUALS(*cstr1, *cstr2);
            ++cstr1;
            ++cstr2;
        }
    }


    void test_BaseVec_constructor_subset()
    { 
        BaseVec bv1 = build_vec("GATC");
        BaseVec bv2(bv1, 1, 3);
        verify_vec(bv2, "ATC");
    }

    void test_BaseVec_constructor_sized()
    { 
        BaseVec_size_type cap = 100;
        BaseVec bv1(cap);
        TS_ASSERT_LESS_THAN_EQUALS(cap, bv1.size());
    }

    void test_BaseVec_constructor_string()
    { 
        BaseVec bv("GATC");
        verify_vec(bv, "GATC");
    }

    void test_BaseVec_IsHomopolymer()
    { 
        BaseVec homovec("AAAA");
        BaseVec heterovec("AGTC");
        TS_ASSERT(homovec.IsHomopolymer());
        TS_ASSERT(!heterovec.IsHomopolymer());
    }

    void test_BaseVec_HomopolPercent()
    { 
        BaseVec homovec("AAAAAAAAAAGGGTTTCCCC");
        pair<float, unsigned char> p;
        p = homovec.HomopolPercent();
        TS_ASSERT(p.first == 50.0);
        TS_ASSERT(p.second == as_char('A'));
    }

    void test_BaseVec_Homopol()
    { 
        BaseVec homovec("AAAGGGGGAAAAA");
        int ret;
        ret = homovec.Homopol(6);
        TS_ASSERT_EQUALS(ret, 5);
    }

    void test_BaseVec_SetFromStringWithNs()
    { 
        String s = "NNNNNNNNNN";
        BaseVec bvec;
        bvec.SetFromStringWithNs(s);
        BaseVec::const_iterator bv_itr = bvec.Begin();
        
        TS_ASSERT_EQUALS(bvec.size(), s.size());
        while(bv_itr != bvec.End())
        {
            const Base& b = Base::fromVal(*bv_itr);
            TS_ASSERT(!b.isAmbiguous());
            bv_itr++;
        }
    }

    void test_BaseVec_Reverse()
    { 
        BaseVec bvec("TAGTATCGGATTTTATCCGA");
        bvec.Reverse();
        verify_vec(bvec, "AGCCTATTTTAGGCTATGAT");
        bvec.Reverse();
        verify_vec(bvec, "TAGTATCGGATTTTATCCGA");
    }

    void test_BaseVec_ReverseComplement()
    {
        BaseVec bvec("TTGGAACC");
        bvec.ReverseComplement();
        verify_vec(bvec, "GGTTCCAA");
        bvec.ReverseComplement();
        verify_vec(bvec, "TTGGAACC");
        BaseVec bvec1("GGTTCAA");
        bvec1.ReverseComplement();
        verify_vec(bvec1, "TTGAACC");
        bvec1.ReverseComplement();
        verify_vec(bvec1, "GGTTCAA");
    }

    void test_BaseVec_Canonicalize()
    { 
        BaseVec bvec("GTAC");
        TS_ASSERT_EQUALS(CanonicalForm::PALINDROME, bvec.Canonicalize());
        bvec = BaseVec("GGTTCCAT");
        TS_ASSERT_EQUALS(CanonicalForm::REV, bvec.Canonicalize());
        // should of RC'd
        verify_vec(bvec, "ATGGAACC");
        TS_ASSERT_EQUALS(CanonicalForm::FWD, bvec.Canonicalize());
    }

    void test_BaseVec_ToString()
    { 
        String s("GATC");
        BaseVec bvec(s);
        TS_ASSERT_EQUALS(bvec.ToString(), s);
    }

    void test_BaseVec_PrintBases()
    { 
        String s_src = "GGGGGG";
        String s_target = "CCCCC\nC\n";
        BaseVec bvec(s_src);
        stringstream os;
        bvec.PrintBases(os, 0, bvec.size(), 1, 5);
        verify_cstr(os.str().c_str(), s_target.c_str());
    }

    void test_BaseVec_Print()
    { 
        String s_src = "GGGGGG";
        stringstream os;
        BaseVec bvec(s_src);
        s_src += "\n";
        bvec.Print(os);
        verify_cstr(os.str().c_str(), s_src.c_str());
    }

    void test_BaseVec_Print_FastA()
    { 
        String s_src = "GGGGGG";
        stringstream os;
        BaseVec bvec(s_src);
        s_src = ">sequence_1\n" + s_src + "\n";
        bvec.Print(os, 1);
        verify_cstr(os.str().c_str(), s_src.c_str());
    }

    void test_BaseVec_PrintCol()
    { 
        String s_src = "GGGGGG";
        String s_target = "GGGGG\nG\n";
        BaseVec bvec(s_src);
        stringstream os;
        bvec.PrintCol(os, 5);
        verify_cstr(os.str().c_str(), s_target.c_str());
    }

    void test_BaseVec_PrintCol_FastA()
    { 
        String s_src = "GGGGGG";
        String s_target = ">testseq\nGGGGG\nG\n";
        BaseVec bvec(s_src);
        stringstream os;
        bvec.PrintCol(os, "testseq", 5);
        verify_cstr(os.str().c_str(), s_target.c_str());
    }

    void test_BaseVec_Find()
    { 
        BaseVec target("AAAGTCAA");
        BaseVec token("AGTC");
        TS_ASSERT_EQUALS(target.Find(token, 0, target.size()), 2u);
        token = BaseVec("GGAC");
        TS_ASSERT_EQUALS(target.Find(token, 0, target.size()), target.size());
    }

    void test_BaseVec_FindAll()
    { 
        typedef vec<unsigned int> VecUInt;
        BaseVec target("AAAGGGAAAGGGAAA");
        BaseVec token("AAA");
        VecUInt ret = target.FindAll(token);
        TS_ASSERT_EQUALS(ret.size(), 3u);
        for(VecUInt::iterator itr = ret.begin(); itr != ret.end(); itr++)
        {
            using std::distance;
            TS_ASSERT_EQUALS(*itr, distance(ret.begin(), itr) * 6);
        }
    }

    void test_BaseVec_Cap()
    { 
        BaseVec src("AAAGGGAAAGGGAAA");
        String s_target("AGAGA");
        src.Cap(1);
        TS_ASSERT_EQUALS(src.ToString(), s_target);
        src.Cap(0);
        TS_ASSERT_EQUALS(src.ToString(), "");
    }

    void test_BaseVec_CopyBases()
    { 
        BaseVec src("AAAA");
        BaseVec target("GGGG");
        CopyBases(src, 2, target, 2, 2, True);
        TS_ASSERT_EQUALS(target.ToString(), String("GGTT"));
    }

    void test_BaseVec_StringReverseComplement()
    { 
        String input = "TTGGAACC";
        String output;
        String expected = "GGTTCCAA";
        StringReverseComplement(input, output);
        TS_ASSERT_EQUALS(output, expected);
    }

    void test_BaseVec_GcBases()
    { 
        BaseVec target("GCATGCAT");
        TS_ASSERT_EQUALS(target.GcBases(0, 4), 2u);
        TS_ASSERT_EQUALS(target.GcBases(4, 8), 2u);
        TS_ASSERT_EQUALS(target.GcBases(0, 8), 4u);
        TS_ASSERT_EQUALS(GcBases(target, 0, 8), 4u);
    }

    void test_BaseVec_SetToSubOf()
    { 
        BaseVec target("");
        BaseVec    src("AGTCCTGA");
        target.SetToSubOf(src, 0, -1);
        verify_vec(target, "AGTCCTGA");

        target.SetToSubOf(src, 0, 8);
        verify_vec(target, "AGTCCTGA");
        target.SetToSubOf(src, 0, 7);
        verify_vec(target, "AGTCCTG");
        target.SetToSubOf(src, 0, 6);
        verify_vec(target, "AGTCCT");
        target.SetToSubOf(src, 0, 5);
        verify_vec(target, "AGTCC");
        target.SetToSubOf(src, 0, 4);
        verify_vec(target, "AGTC");
        target.SetToSubOf(src, 0, 3);
        verify_vec(target, "AGT");
        target.SetToSubOf(src, 0, 2);
        verify_vec(target, "AG");
        target.SetToSubOf(src, 0, 1);
        verify_vec(target, "A");
        target.SetToSubOf(src, 0, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 1, 7);
        verify_vec(target, "GTCCTGA");
        target.SetToSubOf(src, 1, 6);
        verify_vec(target, "GTCCTG");
        target.SetToSubOf(src, 1, 5);
        verify_vec(target, "GTCCT");
        target.SetToSubOf(src, 1, 4);
        verify_vec(target, "GTCC");
        target.SetToSubOf(src, 1, 3);
        verify_vec(target, "GTC");
        target.SetToSubOf(src, 1, 2);
        verify_vec(target, "GT");
        target.SetToSubOf(src, 1, 1);
        verify_vec(target, "G");
        target.SetToSubOf(src, 1, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 2, 6);
        verify_vec(target, "TCCTGA");
        target.SetToSubOf(src, 2, 5);
        verify_vec(target, "TCCTG");
        target.SetToSubOf(src, 2, 4);
        verify_vec(target, "TCCT");
        target.SetToSubOf(src, 2, 3);
        verify_vec(target, "TCC");
        target.SetToSubOf(src, 2, 2);
        verify_vec(target, "TC");
        target.SetToSubOf(src, 2, 1);
        verify_vec(target, "T");
        target.SetToSubOf(src, 2, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 3, 5);
        verify_vec(target, "CCTGA");
        target.SetToSubOf(src, 3, 4);
        verify_vec(target, "CCTG");
        target.SetToSubOf(src, 3, 3);
        verify_vec(target, "CCT");
        target.SetToSubOf(src, 3, 2);
        verify_vec(target, "CC");
        target.SetToSubOf(src, 3, 1);
        verify_vec(target, "C");
        target.SetToSubOf(src, 3, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 4, 4);
        verify_vec(target, "CTGA");
        target.SetToSubOf(src, 4, 3);
        verify_vec(target, "CTG");
        target.SetToSubOf(src, 4, 2);
        verify_vec(target, "CT");
        target.SetToSubOf(src, 4, 1);
        verify_vec(target, "C");
        target.SetToSubOf(src, 4, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 5, 3);
        verify_vec(target, "TGA");
        target.SetToSubOf(src, 5, 2);
        verify_vec(target, "TG");
        target.SetToSubOf(src, 5, 1);
        verify_vec(target, "T");
        target.SetToSubOf(src, 5, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 6, 2);
        verify_vec(target, "GA");
        target.SetToSubOf(src, 6, 1);
        verify_vec(target, "G");
        target.SetToSubOf(src, 6, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 7, 1);
        verify_vec(target, "A");
        target.SetToSubOf(src, 7, 0);
        verify_vec(target, "");

        target.SetToSubOf(src, 8, 0);
        verify_vec(target, "");
    }

    void test_BaseVec_GcParent()
    { 
        BaseVec target("AGCT");
        TS_ASSERT_EQUALS(GcPercent("AGCT"), 50.0)
#ifndef TEST_BASEVECTOR
        TS_ASSERT_EQUALS(target.GcPercent(0, 4), 50.0)
#endif
    }

    void test_BaseVec_Overlap()
    { 
        BaseVec target("AAAATTTT");
        BaseVec otest("TTTTT");
        TS_ASSERT(Overlap(target, otest, 4));
        TS_ASSERT(!Overlap(target, otest, 5));
#ifndef TEST_BASEVECTOR
        TS_ASSERT(target.Overlap(otest, 4));
#endif
    }

    void test_BaseVec_LargestOverlap()
    { 
        BaseVec target("AAAATTTT");
        BaseVec otest("TTTTT");
        TS_ASSERT(LargestOverlap(target, otest));
    }

    void test_BaseVec_Cat()
    { 
        BaseVec b1("AAAA");
        BaseVec b2("GGGG");
        BaseVec b3("TTTT");
        BaseVec target;
        target = Cat(b1, b2);
        verify_vec(target, "AAAAGGGG");
        target = Cat(b1, b2, b3);
        verify_vec(target, "AAAAGGGGTTTT");
    }

    void test_Kmers()
    {
        BaseVec bv("ACCGTAACGCCCTATGCCTAACGT");
        unsigned int kmer = 0x0190E;
        TS_ASSERT_EQUALS(bv.extractKmer(3,7),kmer);
        bv.assignBaseBits(7,reinterpret_cast<unsigned char const*>(&kmer));
        TS_ASSERT_EQUALS(bv,BaseVec("GTAACGC"));
    }
};

#endif
