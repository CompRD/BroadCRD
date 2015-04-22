/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file MD5Test.cc
 * \author tsharpe
 * \date Dec 9, 2008
 *
 * \brief Tests the MD5 implementation.
 */
#include <iostream>
#include "util/MD5.h"
#include "system/Assert.h"
#include "system/Types.h"

int main( int argc, char** argv )
{
    AssertEq(sizeof(MD5::uint),4u); // 4-byte ints
    AssertEq(sizeof(MD5::ulong),8u); // 8-byte longs
    AssertEq(*reinterpret_cast<unsigned int const*>("0123"),0x33323130u); // little endian

    char const*const testStrings[] =
    {
        "",
        "a",
        "abc",
        "message digest",
        "abcdefghijklmnopqrstuvwxyz",
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
        "12345678901234567890123456789012345678901234567890123456789012345678901234567890"
    };

    char const*const resultStrings[] =
    {
        "d41d8cd98f00b204e9800998ecf8427e",
        "0cc175b9c0f1b6a831c399e269772661",
        "900150983cd24fb0d6963f7d28e17f72",
        "f96b697d7cb7938d525a2f31aaf161d0",
        "c3fcd3d76192e4007dfb496cca67e13b",
        "d174ab98d277d9f5a5611c2c9f419d9f",
        "57edf4a22be3c955ac49da2e2107b67a"
    };

    uint const N_TESTS = sizeof(testStrings)/sizeof(char*);

    bool aOK = true;
    MD5 md5;
    for ( unsigned int iii = 0; iii < N_TESTS; ++iii )
    {
        char const* test = testStrings[iii];
        char const* result = resultStrings[iii];
        md5.reInit();
        md5.update(test,test+strlen(test));
        char const* digest = md5.getHexDigest();
        bool ok = !memcmp(result,digest,32);
        aOK = aOK && ok;
        std::cout << test << "\t" << result << "\t" << digest << "\t" << (ok?"OK":"NFG") << std::endl;
        md5.reInit();
        md5.update(test,strlen(test));
        digest = md5.getHexDigest();
        ok = !memcmp(result,digest,32);
        aOK = aOK && ok;
        std::cout << test << "\t" << result << "\t" << digest << "\t" << (ok?"OK":"NFG") << std::endl;
    }
    md5.reInit();
    char a = 'a';
    for ( int iii = 0; iii < 1000000; ++iii )
    {
        md5.update(a);
    }
    char const* digest = md5.getHexDigest();
    char const* result = "7707d6ae4e027c70eea2a935c2296f21";
    bool ok = !memcmp(result,digest,32);
    aOK = aOK && ok;
    cout << "a million a's\t" << result << "\t" << digest << "\t" << (ok?"OK":"NFG") << endl;

    return !aOK;
}
