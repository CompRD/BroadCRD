///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * TestMemoryMonitor.cc
 *
 *  Created on: Mar 25, 2013
 *      Author: tsharpe
 */

#include "system/MemoryMonitor.h"
#include <cstring>
#include <iostream>
#include <unistd.h>

int main( int argc, char** argv )
{
    char const* log = "log";
    if ( argc == 1 )
        log = nullptr;
    MemoryMonitor mm(1u,log);
    mm.trackUse();
    size_t const LEN = 1ul << 30;
    char* a1 = new char[LEN];
    memset(a1,255,LEN);
    sleep(10);
    std::cout << mm.getUse() << std::endl;
    char* a2 = new char[LEN];
    memset(a2,255,LEN);
    sleep(10);
    std::cout << mm.getUse() << std::endl;
    delete [] a2;
    delete [] a1;
}
