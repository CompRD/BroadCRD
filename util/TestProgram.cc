///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library GMPXX

// require <cstddef> to use gmp in GCC-4.9
#include <cstddef>
#include <gmpxx.h>

#include <iostream>

int main( int argc, char**argv )
{
    mpz_t anInteger;
    mpz_init(anInteger);

    std::cout << "Test example compiled and built successfully!\n";

    return 0;
}
