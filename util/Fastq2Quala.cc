///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Fastq2Quala.cc
 * \author tsharpe
 * \date Aug 12, 2009
 *
 * \brief Converts a Fastq file to FASTA file of quals.
 */
#include "system/Assert.h"
#include "system/RunTime.h"
#include <iostream>
#include <string>

using std::cin;
using std::cout;
using std::endl;

#define OFFSET 64
// if you have Sanger (or NCBI) format data, the OFFSET should be 33 instead
// just try it and if you get negative quality scores, change 64 to 33

int main( int argc, char **argv )
{
    RunTime();

    char buf[1024];
    while ( cin.getline(buf, sizeof(buf)) )
    {
        AssertEq(buf[0],'@');
        buf[0] = '>';
        cout << buf << '\n';
        cin.ignore(1024, '\n');
        cin.ignore(1024, '\n');
        cin.getline(buf, sizeof(buf));
        for ( char* ppp = buf; *ppp; ++ppp )
        {
            cout << (*ppp - OFFSET);
            cout << (ppp[1] ? ' ' : '\n');
        }
    }
}
