/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
* \file FeudalStringTests.h
* \author ghall
* \date Oct 1, 2009
*
* \brief Tests for FeudalString (as CharString)
*/
#ifndef __FEUDALSTRING_TESTS_H
#define __FEUDALSTRING_TESTS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cxxtest/TestSuite.h>
#include "dna/Bases.h"
#include "feudal/CharString.h"
#include "system/ParsedArgs.h"

class StringTests : public CxxTest::TestSuite
{
public:
    static const String globalTest;

    void setUp()
    { }

    void tearDown()
    { }

    void verify_str(const String& tstr, const char* target)
    {
        const char *actual = tstr.c_str();
        TS_ASSERT_EQUALS(strlen(actual), strlen(target));
        while(*actual && *target)
        {
            TS_ASSERT_EQUALS(*actual++, *target++);
        }
    }

    void test_String_Ctor()
    { 
        const char hworld[] = "Hello World!";
        String tstr1(hworld);
        /* FeudalString() */
        { String tstr1; verify_str(tstr1, ""); }
        /* FeudalString(const char* cstr) */
        { String tstr2(hworld); verify_str(tstr2, hworld); }
        /* FeudalString(Itr first, Itr const& last) */
        { String tstr2(tstr1.begin(), tstr1.end()); verify_str(tstr2, "Hello World!"); }
        /* FeudalString(const FeudalString& str, size_type pos) */
        { String tstr2(tstr1, 6); verify_str(tstr2, "World!"); }
        /* FeudalString(const FeudalString& str, size_type pos, size_type n) */
        { String tstr2(tstr1, 6, 5); verify_str(tstr2, "World"); }
        /* FeudalString(const char* cstr, size_type n) */
        { String tstr2(hworld, 5); verify_str(tstr2, "Hello"); }
        /* FeudalString(const char* cstr, size_type n) */
        { String tstr2(hworld, 6); verify_str(tstr2, "Hello "); }
        /* FeudalString(const char* cstr, size_type n) */
        { String tstr2(hworld + 6, 5); verify_str(tstr2, "World"); }
        /* FeudalString(value_type c, const Allocator& a = Allocator()) */
        { String tstr2(4, 'H'); verify_str(tstr2, "HHHH"); }
        /* FeudalString(size_type n, const Allocator& a = Allocator()) */
        { String tstr2(4); verify_str(tstr2, ""); }
    }

    /*
     * test STL methods
     */

    void test_String_operator_equals()
    { 
        String tstr1;
        { tstr1 = "Hello!"; verify_str(tstr1, "Hello!"); }
        { tstr1 = "Goodbye!"; verify_str(tstr1, "Goodbye!"); }
        { tstr1 = String("Allo?"); verify_str(tstr1, "Allo?"); }
        { tstr1 = 'c'; verify_str(tstr1, "c"); }
    }

    void test_String_operator_IN()
    { 
        String tstr1 = "  Hello World!  ";
        std::istringstream ss(tstr1.c_str());
        String target;
        ss >> target;
        verify_str(target, "Hello");
    }

    void test_String_operator_OUT()
    { 
        String tstr1 = "Hello World!";
        std::ostringstream ss;
        ss << tstr1;
        String target(ss.str().c_str());
        verify_str(target, "Hello World!");
    }

    void test_String_c_str()
    { 
        String tstr("");
        verify_str(tstr, "");
        tstr = "Hello World!";
        verify_str(tstr, "Hello World!");
        tstr.push_back('!');
        verify_str(tstr, "Hello World!!");
    }

    void test_String_concat()
    { 
        {   String tstr1 = "Hello";
            String tstr2 = "World!";
            tstr1 += string(" ");
            tstr1 += tstr2;
            verify_str(tstr1, "Hello World!"); }
        {   String tstr1 = "Hello";
            String tstr2 = "World";
            String tstr3 = tstr1 + ' ' + tstr2 + string("!");
            verify_str(tstr3, "Hello World!"); }
        {   string FILENAME = "DIARY.TXT";
            String str = String("gzip -1 > ") + FILENAME + ".gz";
            verify_str(str, "gzip -1 > DIARY.TXT.gz"); }
    }

    void test_String_append()
    { 
        // append(const FeudalString& str)
        {   String tstr1 = "Hello";
            String tstr2 = " World!";
            tstr1.append(tstr2);
            verify_str(tstr1, "Hello World!"); }
        // append(const FeudalString& str, size_type pos, size_type n)
        {   String tstr1 = "Hello";
            tstr1.append(" World!", 1, 5);
            verify_str(tstr1, "HelloWorld"); }
        // append(const FeudalString& str, size_type pos, size_type n)
        {   String tstr1 = "Hello";
            String tstr2 = " World!";
            tstr1.append(tstr2, 1, 5);
            verify_str(tstr1, "HelloWorld"); }
        // append(const_pointer cstr, size_type n)
        {   String tstr1 = "Hello";
            tstr1.append(" World!", 6);
            verify_str(tstr1, "Hello World"); }
        // append(const_pointer cstr)
        {   String tstr1 = "Hello";
            tstr1.append(" World!");
            verify_str(tstr1, "Hello World!"); }
        // append(size_type sz, const value_type c)
        {   String tstr1 = "Hello";
            tstr1.append(4, 'H');
            verify_str(tstr1, "HelloHHHH"); }
        // assign(InputIterator first, InputIterator last)
        {   String tstr1 = "Hello";
            String tstr2 = " World!";
            tstr1.append(tstr2.begin(), tstr2.end());
            verify_str(tstr1, "Hello World!"); }
    }

    void test_String_assign()
    { 
        // assign(const FeudalString& str)
        {   String tstr1 = "not me, not me, not me, not me";
            String tstr2 = "Hello World!";
            tstr1.assign(tstr2);
            verify_str(tstr1, "Hello World!"); }
        // assign(const FeudalString& str, size_type pos, size_type n)
        {   String tstr1 = "not me, not me, not me, not me";
            String tstr2 = "Hello World!";
            tstr1.assign(tstr2, 6, 5);
            verify_str(tstr1, "World"); }
        // assign(const_pointer cstr, size_type n)
        {   String tstr1 = "not me, not me, not me, not me";
            String tstr2 = "Hello World!";
            tstr1.assign("Hello World!", 5);
            verify_str(tstr1, "Hello"); }
        // assign(const_pointer cstr)
        {   String tstr1 = "not me, not me, not me, not me";
            tstr1.assign("Hello World!");
            verify_str(tstr1, "Hello World!"); }
        // assign(size_type n, value_type c)
        {   String tstr1 = "not me, not me, not me, not me";
            tstr1.assign(4, 'H');
            verify_str(tstr1, "HHHH"); }
        // assign(InputIterator first, InputIterator last)
        {   String tstr1 = "not me, not me, not me, not me";
            String tstr2 = "Hello World!";
            tstr1.assign(tstr2.begin(), tstr2.end());
            verify_str(tstr1, "Hello World!"); }
    }

    void test_String_insert()
    { 
        // insert(size_type pos1, const FeudalString& str)
        {   String tstr1 = "Hello !";
            String tstr2 = "World";
            tstr1.insert(6, tstr2);
            verify_str(tstr1, "Hello World!"); }
        // insert(size_type pos1, const FeudalString& str, 
        //                      size_type pos2, size_type n)
        {   String tstr1 = "Hello !";
            String tstr2 = "!World!";
            tstr1.insert(6, tstr2, 1, 5);
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& insert(size_type pos1, const charT* s, size_type n)
        {   String tstr1 = "Hello !";
            const char *cstr = "World!!";
            tstr1.insert(6, cstr, 5);
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& insert(size_type pos1, const T* cstr)
        {   String tstr1 = "Hello !";
            tstr1.insert(6, 4, 'H');
            verify_str(tstr1, "Hello HHHH!"); }
        // FeudalString& insert(size_type pos1, size_type n, value_type c)
        {   String tstr1 = "Hello !";
            const char *cstr = "World";
            tstr1.insert(6, cstr);
            verify_str(tstr1, "Hello World!"); }
        // iterator insert(iterator p, value_type c)
        {   String tstr1 = "Hello !";
            tstr1.insert(tstr1.begin(6), 'W');
            verify_str(tstr1, "Hello W!"); }
        // void insert(iterator p, size_type n, value_type c)
        {   String tstr1 = "Hello !";
            tstr1.insert(tstr1.begin(6), 4, 'H');
            verify_str(tstr1, "Hello HHHH!"); }
        // void insert(iterator p, InputIterator first, InputIterator last)
        {   String tstr1 = "Hello !";
            String tstr2 = "World";
            tstr1.insert(tstr1.begin(6), tstr2.begin(), tstr2.end());
            verify_str(tstr1, "Hello World!"); }
    }

    void test_String_erase()
    { 
        // FeudalString& erase(size_type pos = 0, size_type n = npos)
        {   String tstr1 = "## Hello World! ##";
            tstr1.erase(0, 3);
            tstr1.erase(tstr1.size() - 3);
            verify_str(tstr1, "Hello World!"); }
        // iterator erase(iterator p)
        {   String tstr1 = "#Hello World!#";
            tstr1.erase(tstr1.begin());
            tstr1.erase(tstr1.begin(tstr1.size()-1));
            verify_str(tstr1, "Hello World!"); }
        // iterator erase(iterator first, iterator last)
        {   String tstr1 = "Hello ###World!";
            tstr1.erase(tstr1.begin(6), tstr1.begin(9));
            verify_str(tstr1, "Hello World!"); }
    }

    void test_String_replace()
    { 
        // FeudalString& replace(size_type pos1, size_type n1, 
        //                                  const FeudalString& str)
        {   String tstr1 = "Hello Moon!";
            String tstr2 = "World";
            tstr1.replace(6, 4, tstr2);
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& replace(iterator i1, iterator i2, 
        //                                  const FeudalString& str)
        {   String tstr1 = "Hello Moon!";
            String tstr2 = "World";
            tstr1.replace(tstr1.begin(6), tstr1.begin(10), tstr2);
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& replace(size_type pos1, size_type n1, 
        //                          const FeudalString& str,
        //                          size_type pos2, size_type n2)
        {   String tstr1 = "Hello Moon!";
            String tstr2 = "#World#";
            tstr1.replace(6, 4, tstr2, 1, 5);
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& replace(size_type pos1, size_type n1, 
        //                          const T* cstr, size_type n2)
        {   String tstr1 = "Hello Moon!";
            tstr1.replace(6, 4, "World$", 5);
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& replace(iterator i1, iterator i2, 
        //                          const T* s, size_type n2)
        {   String tstr1 = "Hello Moon!";
            tstr1.replace(tstr1.begin(6), tstr1.begin(10), "World$", 5);
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& replace(size_type pos1, size_type n1, const T* cstr)
        {   String tstr1 = "Hello Moon!";
            tstr1.replace(6, 4, "World");
            verify_str(tstr1, "Hello World!"); }
        // FeudalString& replace(iterator i1, iterator i2, const T* cstr)
        {   String tstr1 = "Hello Moon!";
            tstr1.replace(tstr1.begin(6), tstr1.begin(10), "World");
            verify_str(tstr1, "Hello World!"); }
        //  FeudalString& replace(size_type pos1, size_type n1, 
        //                          size_type n2, value_type c)
        {   String tstr1 = "Hello Moon!";
            tstr1.replace(6, 4, 4, 'H');
            verify_str(tstr1, "Hello HHHH!"); }
        // FeudalString& replace(iterator i1, iterator i2, 
        //                          size_type n2, value_type c)
        {   String tstr1 = "Hello Moon!";
            tstr1.replace(tstr1.begin(6), tstr1.begin(10), 4, 'H');
            verify_str(tstr1, "Hello HHHH!"); }
        // FeudalString& replace(iterator i1, iterator i2, 
        //              InputIterator first, InputIterator last)
        {   String tstr1 = "Hello Moon!";
            String tstr2 = "World";
            tstr1.replace(tstr1.begin(6), tstr1.begin(10), tstr2.begin(), tstr2.end());
            verify_str(tstr1, "Hello World!"); }
    }

    void test_String_find()
    { 
        {   String line = "zzzizzzz";
            TS_ASSERT_EQUALS(line.find("zzz", 5u), 5u); }
        {   String line = "> love";
            TS_ASSERT(line.Contains(">", 0)); }
        // size_type find(const string& str, size_type pos = 0) const
        {   String tstr1 = "Hello World!";
            String needle = "World";
            TS_ASSERT_EQUALS(tstr1.find(needle), 6u);
            needle = "World!";
            TS_ASSERT_EQUALS(tstr1.find(needle), 6u);
            needle = "Hello";
            TS_ASSERT_EQUALS(tstr1.find(needle), 0u);
            needle = "Hello World!";
            TS_ASSERT_EQUALS(tstr1.find(needle), 0u);
            needle = "";
            TS_ASSERT_EQUALS(tstr1.find(needle), 0u);
            needle = "NONSENSE";
            TS_ASSERT_EQUALS(tstr1.find(needle), String::npos); }
        // size_type find(const char* cstr, size_type pos, size_type n) const
        {   String tstr1 = "Hello World!";
            TS_ASSERT_EQUALS(tstr1.find("World#", 0, 5), 6u); }
        // size_type find(const char* cstr, size_type pos = 0) const
        {   String tstr1 = "Hello World!";
            TS_ASSERT_EQUALS(tstr1.find("World"), 6u); }
        // size_type find(value_type c, size_type pos = 0) const
        {   String tstr1 = "Hello World!";
            TS_ASSERT_EQUALS(tstr1.find('W'), 6u); }
        // npos
        {   String tstr1 = "Hello World!";
            TS_ASSERT_EQUALS(tstr1.find("Moon"), String::npos); }
    }

    void test_String_rfind()
    { 
        // size_type rfind(const FeudalString& str, size_type pos = npos) const
        {   String tstr1 = "Welcome, Welcome...";
            String needle = "Welcome";
            TS_ASSERT_EQUALS(tstr1.rfind(needle), 9u);
            needle = "Welcome...";
            TS_ASSERT_EQUALS(tstr1.rfind(needle), 9u);
            needle = "Welcome,";
            TS_ASSERT_EQUALS(tstr1.rfind(needle), 0u);
            needle = "";
            TS_ASSERT_EQUALS(tstr1.rfind(needle), tstr1.size());
            needle = "NONSENSE";
            TS_ASSERT_EQUALS(tstr1.rfind(needle), String::npos); }
        // size_type rfind(const char* cstr, size_type pos, size_type n) const
        {   String tstr1 = "Welcome, Welcome...";
            TS_ASSERT_EQUALS(tstr1.rfind("Welcome!", String::npos, 7), 9u); }
        // size_type rfind(const T* cstr, size_type pos = npos) const
        {   String tstr1 = "Welcome, Welcome...";
            TS_ASSERT_EQUALS(tstr1.rfind("Welcome"), 9u); }
        // size_type rfind(const T* cstr, size_type pos = npos) const
        {   String tstr1 = "Welcome, Welcome...";
            TS_ASSERT_EQUALS(tstr1.rfind('W'), 9u); }
        // npos
        {   String tstr1 = "Welcome, Welcome...";
            String needle = "Goodbye";
            TS_ASSERT_EQUALS(tstr1.rfind(needle), String::npos); }
    }

    void test_String_find_first_of()
    { 
        // size_type rfind(const FeudalString& str, size_type pos = npos) const
        {   String tstr1 = "zzzizzz";
            String needle = "aeiou";
            TS_ASSERT_EQUALS(tstr1.find_first_of(needle), 3u);
        }
        // size_type find_first_of(const charT* cstr, size_type pos, size_type n) const
        {   String tstr1 = "zzzizzz";
            TS_ASSERT_EQUALS(tstr1.find_first_of("aeiouz", 0, 5u), 3u);
        }
        // size_type find_first_of(const charT* cstr, size_type pos = 0) const
        {   String tstr1 = "zzzizzz";
            TS_ASSERT_EQUALS(tstr1.find_first_of("aeiou"), 3u);
        }
        // size_type find_first_of(char c, size_type pos = 0) const
        {   String tstr1 = "zzzizzz";
            TS_ASSERT_EQUALS(tstr1.find_first_of('i'), 3u);
        }
        // npos
        {   String tstr1 = "zzzizzz";
            String needle = "nope";
            TS_ASSERT_EQUALS(tstr1.find_first_of(needle), String::npos);
        }
    }

    void test_String_find_last_of()
    { 
        // size_type rfind(const FeudalString& str, size_type pos = npos) const
        {   String tstr1 = "azzizzz";
            String needle = "aeiou";
            TS_ASSERT_EQUALS(tstr1.find_last_of(needle), 3u);
        }
        // size_type find_last_of(const charT* cstr, size_type pos, size_type n) const
        {   String tstr1 = "azzizzz";
            TS_ASSERT_EQUALS(tstr1.find_last_of("aeiouz", String::npos, 5u), 3u);
        }
        // size_type find_last_of(const charT* cstr, size_type pos = 0) const
        {   String tstr1 = "azzizzz";
            TS_ASSERT_EQUALS(tstr1.find_last_of("aeiou"), 3u);
        }
        // size_type find_last_of(char c, size_type pos = 0) const
        {   String tstr1 = "azzizzz";
            TS_ASSERT_EQUALS(tstr1.find_last_of('i'), 3u);
        }
        // npos
        {   String tstr1 = "azzizzz";
            String needle = "nope";
            TS_ASSERT_EQUALS(tstr1.find_last_of(needle), String::npos);
        }
    }

    void test_String_find_first_not_of()
    { 
        // size_type find_first_not_of(const FeudalString& str, size_type pos = npos) const
        {   String tstr1 = "zzzizzz";
            String needle = "abz";
            TS_ASSERT_EQUALS(tstr1.find_first_not_of(needle), 3u);
        }
        // size_type find_first_not_of(const charT* cstr, size_type pos, size_type n) const
        {   String tstr1 = "zzzizzz";
            TS_ASSERT_EQUALS(tstr1.find_first_not_of("abzi", 0, 3u), 3u);
        }
        // size_type find_first_not_of(const charT* cstr, size_type pos = 0) const
        {   String tstr1 = "zzzizzz";
            TS_ASSERT_EQUALS(tstr1.find_first_not_of("abz"), 3u);
        }
        // size_type find_first_not_of(char c, size_type pos = 0) const
        {   String tstr1 = "zzzizzz";
            TS_ASSERT_EQUALS(tstr1.find_first_not_of('z'), 3u);
        }
        // npos
        {   String tstr1 = "zzzizzz";
            String needle = "zi";
            TS_ASSERT_EQUALS(tstr1.find_first_not_of(needle), String::npos);
        }
    }

    void test_String_find_last_not_of()
    { 
        // size_type rfind(const FeudalString& str, size_type pos = npos) const
        {   String tstr1 = "azzizzz";
            String needle = "abz";
            TS_ASSERT_EQUALS(tstr1.find_last_not_of(needle), 3u);
        }
        // size_type find_last_not_of(const charT* cstr, size_type pos, size_type n) const
        {   String tstr1 = "azzizzz";
            TS_ASSERT_EQUALS(tstr1.find_last_not_of("abzi", String::npos, 3u), 3u);
        }
        // size_type find_last_not_of(const charT* cstr, size_type pos = 0) const
        {   String tstr1 = "azzizzz";
            TS_ASSERT_EQUALS(tstr1.find_last_not_of("abz"), 3u);
        }
        // size_type find_last_not_of(char c, size_type pos = 0) const
        {   String tstr1 = "azzizzz";
            TS_ASSERT_EQUALS(tstr1.find_last_not_of('z'), 3u);
        }
        // npos
        {   String tstr1 = "azzizzz";
            String needle = "aziz";
            TS_ASSERT_EQUALS(tstr1.find_last_not_of(needle), String::npos);
        }
    }

    void test_String_compare()
    {
        {   String down = "abc";
            String up = "ABC";
            String up2 = "ABC1";
            TS_ASSERT(down.compare(up) > 0);
            TS_ASSERT(up.compare(down) < 0);
            TS_ASSERT_EQUALS(up.compare(up), 0);
            TS_ASSERT(up2.compare(up)>0);
            TS_ASSERT(up.compare(up2)<0);
        }
        {   String down = "abc";
            String up = "ABC";
            TS_ASSERT(down.compare("ABC") > 0);
            TS_ASSERT(up.compare("abc") < 0);
            TS_ASSERT_EQUALS(up.compare("ABC"), 0);
        }
        {   String down = "abc";
            String up = "ABC";
            TS_ASSERT(down.compare("ABC") > 0);
            TS_ASSERT(up.compare("abc") < 0);
            TS_ASSERT_EQUALS(up.compare("ABC"), 0);
        }
        {   String down = "abcde";
            String up = "abcde123";
            TS_ASSERT(up.compare(0, 4, down) < 0);
            TS_ASSERT_EQUALS(up.compare(0, 5, down), 0);
            TS_ASSERT(up.compare(0, 6, down) > 0);
        }
        {   String down = "abcde987";
            String up = "abcde123";
            TS_ASSERT_EQUALS(down.compare(0, 4, up, 0, 4), 0);
            TS_ASSERT(down.compare(0, 5, up, 0, 4) > 0);
            TS_ASSERT(down.compare(0, 4, up, 0, 5) < 0);
        }
        {   String down = "abcde987";
            TS_ASSERT_EQUALS(down.compare(0, 4, "abcde123", 4), 0);
            TS_ASSERT(down.compare(0, 5, "abcde123", 4) > 0);
            TS_ASSERT(down.compare(0, 4, "abcde123", 5) < 0);
        }
    }


    void test_String_copy()
    { 
        String tstr1 = "!Hello World!";
        char target[256];
        tstr1.copy(target, 12, 1);
        target[12] = 0;
        String tstr2(target);
        verify_str(tstr2, "Hello World!");
    }

    void test_String_swap()
    { 
        String tstr1 = "Hello World!";
        String tstr2 = "Goodbye World!";
        tstr1.swap(tstr2);
        verify_str(tstr2, "Hello World!");
    }

    void test_String_operator_getline()
    { 
        String target;
        String tstr1 = "1\r\n2\n3\r4";
        std::istringstream ss(tstr1.c_str());
        TS_ASSERT(getline(ss, target));
        verify_str(target, "1");
        TS_ASSERT(getline(ss, target));
        verify_str(target, "2");
        TS_ASSERT(getline(ss, target));
        verify_str(target, "3");
        TS_ASSERT(!getline(ss,target));
        verify_str(target, "4");
    }

    /*
     * test add on Broad methods
     */

    void test_String_operator_string()
    { 
        String tstr1 = "Hello World!";
        std::string s = tstr1;
        String tstr2 = s.c_str();
        verify_str(tstr2, "Hello World!");
    }

    void test_String_StartsWith()
    { 
        {   String tstr1 = "Hello World!";
            String needle = "Hello ";
            TS_ASSERT(tstr1.StartsWith(needle));
            needle.push_back('!');
            TS_ASSERT(!tstr1.StartsWith(needle));
        }
    }

    void test_String_EndsWith()
    { 
        {   String tstr1 = "Hello World!";
            String needle = "World!";
            TS_ASSERT(tstr1.EndsWith(needle));
            needle.push_back('!');
            TS_ASSERT(!tstr1.EndsWith(needle));
        }
    }

    void test_String_Contains()
    { 
        String tstr1 = "Hello World!";
        String needle = "ello ";
        String space = " ";
        TS_ASSERT(tstr1.Contains(needle));
        TS_ASSERT(tstr1.Contains(space)); 
        TS_ASSERT(tstr1.Contains(needle, 1));
        TS_ASSERT(!tstr1.Contains(needle, 2));
        TS_ASSERT(!tstr1.Contains(space, 0)); // Contains is position specific
        TS_ASSERT(tstr1.Contains(space, 5));  // Contains is position specific
        needle.push_back('!');
        TS_ASSERT(!tstr1.Contains(needle));
        needle = "World!";
        TS_ASSERT(tstr1.Contains(needle, -1));
        needle.push_back('!');
        TS_ASSERT(!tstr1.Contains(needle, -1));
        TS_ASSERT(!String("foo").Contains("xyzzy",0));
    }

    void test_String_Position()
    { 
        String tstr1 = "Hello World!";
        TS_ASSERT_EQUALS(tstr1.Position(' '), 5);
        TS_ASSERT_EQUALS(tstr1.Position(" ", 4), -1);
        TS_ASSERT_EQUALS(tstr1.Position(" ", 6), 5);
        TS_ASSERT_EQUALS(tstr1.Position(" ", 5), -1);
        TS_ASSERT_EQUALS(tstr1.Position("Hello World!"), 0);
        TS_ASSERT_EQUALS(tstr1.Position("ello World!", 11), -1);
        TS_ASSERT_EQUALS(tstr1.Position("ello World!", 12), 1);
        TS_ASSERT_EQUALS(tstr1.Position('#'), static_cast<int>(String::npos));
        TS_ASSERT_EQUALS(tstr1.Position('#'), -1);
        TS_ASSERT(tstr1.Position('#') < 0);
    }

    void test_String_Before()
    { 
        {
            String tstr1 = "Hello World! THIS IS GARBAGE";
            String tstr2 = tstr1.Before(" THIS IS GARBAGE");
            verify_str(tstr2, "Hello World!");
        }
        {
            String tstr1 = "genome.fastb";
            String tstr2 = tstr1.SafeBefore(".fastb");
            verify_str(tstr2, "genome");
        }

    }

    void test_String_After()
    { 
        {
            String tstr1 = "One Two Three";
            String tstr2 = tstr1.After("One ");
            verify_str(tstr2, "Two Three");
        }
        {
            String tstr1 = "genome.fastb";
            String tstr2 = tstr1.SafeAfter("genome.");
            verify_str(tstr2, "fastb");
        }
    }

    void test_String_DeleteLeading()
    { 
        String tstr1 = "  One Two Three";
        DeleteLeadingWhiteSpace(tstr1);
        verify_str(tstr1, "One Two Three");
        tstr1 = "   ";
        DeleteLeadingWhiteSpace(tstr1);
        verify_str(tstr1, "");
    }

    void test_String_DeleteTrailing()
    { 
        String tstr1 = "One Two Three  ";
        DeleteTrailingWhiteSpace(tstr1);
        verify_str(tstr1, "One Two Three");
        tstr1 = "   ";
        DeleteTrailingWhiteSpace(tstr1);
        verify_str(tstr1, "");
    }

    void test_String_IsInt()
    { 
        { String tstr1 = ""; TS_ASSERT(tstr1.IsInt()); }
        { String tstr1 = "10"; TS_ASSERT(tstr1.IsInt()); }
        { String tstr1 = "-10K"; TS_ASSERT(tstr1.IsInt()); }
        { String tstr1 = "+1000M"; TS_ASSERT(tstr1.IsInt()); }
        { String tstr1 = "X"; TS_ASSERT(!tstr1.IsInt()); }
        { String tstr1 = "++Hi MoM"; TS_ASSERT(!tstr1.IsInt()); }
    }

    void test_String_Int()
    { 
        { String tstr1 = "123"; TS_ASSERT_EQUALS(tstr1.Int(),123); }
        { String tstr1 = "-10K"; TS_ASSERT_EQUALS(tstr1.Int(), -10000); }
        { String tstr1 = "+10G"; TS_ASSERT_EQUALS(tstr1.Int(), 10e9); }
        { String tstr1 = "+1000M"; TS_ASSERT_EQUALS(tstr1.Int(), 1000e6); }
    }

    void test_String_IsDouble()
    { 
        { String tstr1 = ""; TS_ASSERT(tstr1.IsDouble()); }
        { String tstr1 = "0.123"; TS_ASSERT(tstr1.IsDouble()); }
        { String tstr1 = "-10.1K"; TS_ASSERT(tstr1.IsDouble()); }
        { String tstr1 = "+1000.1M"; TS_ASSERT(tstr1.IsDouble()); }
        { String tstr1 = "+1000M"; TS_ASSERT(tstr1.IsDouble()); }
        { String tstr1 = "X"; TS_ASSERT(!tstr1.IsDouble()); }
        { String tstr1 = "++Hi MoM"; TS_ASSERT(!tstr1.IsDouble()); }
    }

    void test_String_Double()
    { 
        { String tstr1 = "0."; TS_ASSERT_EQUALS(tstr1.Double(), 0.); }
        { String tstr1 = "-10.1K"; 
          TS_ASSERT_EQUALS(tstr1.Double(), 1e3 * atof("-10.1")); }
        { String tstr1 = "+10.1G"; 
          TS_ASSERT_EQUALS(tstr1.Double(), 1e9 * atof("10.1")); }
        { String tstr1 = "+1000.1M"; 
          TS_ASSERT_EQUALS(tstr1.Double(), 1e6 * atof("1000.1")); }
        { String tstr1 = "1500.1%"; 
          double res = .01 * atof("1500.1");
          TS_ASSERT_EQUALS((long) tstr1.Double(), (long) res); }
    }

    void test_String_IsBool()
    { 
        { String tstr1 = "False"; TS_ASSERT(tstr1.IsBool()); }
        { String tstr1 = "True"; TS_ASSERT(tstr1.IsBool()); }
        { String tstr1 = "bla"; TS_ASSERT(!tstr1.IsBool()); }
    }

    void test_String_ToBool()
    { 
        { String tstr1 = "False"; TS_ASSERT(!tstr1.ToBool()); }
        { String tstr1 = "True"; TS_ASSERT(tstr1.IsBool()); }
    }

    void test_String_ToLower()
    { 
        { String tstr1 = "FALSE"; tstr1.ToLower(); verify_str(tstr1, "false"); }
        { String tstr1 = "false"; tstr1.ToLower(); verify_str(tstr1, "false"); }
    }

    void test_String_ToUpper()
    { 
        { String tstr1 = "false"; tstr1.ToUpper(); verify_str(tstr1, "FALSE"); }
        { String tstr1 = "FALSE"; tstr1.ToUpper(); verify_str(tstr1, "FALSE"); }
    }

    void test_String_Freq()
    { 
        {   String tstr1 = "zzzizzzz";
            TS_ASSERT_EQUALS(tstr1.Freq("zzz"), 3u);
            TS_ASSERT_EQUALS(tstr1.Freq("zzzz"), 1u);
            TS_ASSERT_EQUALS(tstr1.Freq("a"), 0u);
        }
    }


    void test_String_evaluate()
    {
        {   String expr = "8+20/5";
            TS_ASSERT_EQUALS(Evaluate(expr), "12");
            expr = "-20*8/4*(30-4+2)";
            TS_ASSERT_EQUALS(Evaluate(expr), "-1120");
        }
    }

    void test_String_GlobalReplaceBy()
    {
        {   String replaceMe = "I am a bad, bad, string.";
            replaceMe.GlobalReplaceBy("bad", "sad");
            verify_str(replaceMe, "I am a sad, sad, string.");
        }
    }

    void test_String_ReplaceBy()
    {
        {   String replaceMe = "I am a bad, bad, string.";
            replaceMe.ReplaceBy("bad", "sad");
            verify_str(replaceMe, "I am a sad, bad, string.");
        }
    }

    void test_String_ToStringAbbrev()
    {
        long num = 1234567890;
        String val = ToStringAbbrev(num);
        verify_str(val, "1.23G");
        num = 123;
        val = ToStringAbbrev(num);
        verify_str(val, "123");
        num = 1234;
        val = ToStringAbbrev(num);
        verify_str(val, "1K");
    }

    void test_String_ToStringAddCommas()
    { 
        long num = 1234567890;
        String val = ToStringAddCommas(num);
        verify_str(val, "1,234,567,890");
        num = 123;
        val = ToStringAddCommas(num);
        verify_str(val, "123");
        num = 1234;
        val = ToStringAddCommas(num);
        verify_str(val, "1,234");
    }

    void test_String_Trim()
    { 
        String trimMe = "  Hello World!  ";
        verify_str(trimMe.Trim(" "), "Hello World!");
        trimMe = "Hello World!  ";
        verify_str(trimMe.Trim(" "), "Hello World!");
        trimMe = "  Hello World!";
        verify_str(trimMe.Trim(" "), "Hello World!");
        trimMe = "";
        verify_str(trimMe.Trim(" "), "");
        trimMe = " ";
        verify_str(trimMe.Trim("!"), " ");
        trimMe = "           ";
        verify_str(trimMe.Trim(" "), "");
    }

    void test_ReplaceExtension()
    {
        String fn = "DIARY.TXT";
        verify_str(fn.ReplaceExtension("TXT", "EXE"), "DIARY.EXE");
    }
};

const String StringTests::globalTest = "This is a test";

#endif // __FEUDALSTRING_TESTS_H
