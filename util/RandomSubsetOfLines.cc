/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// RandomSubsetOfLines: read lines from stdin and return a random subset thereof.

#include <cstdlib>
#include <fstream>

#include "Basevector.h"
#include "MainTools.h"

int main( int argc, char *argv[] )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_Int(N);
    CommandArgument_String_OrDefault(IN, "/dev/stdin");
    CommandArgument_String_OrDefault(OUT, "/dev/stdout");
    CommandArgument_Int_OrDefault(SEED, 0);
    EndCommandArguments;

    srand(SEED);

    std::ifstream in(IN.c_str());
    ForceAssert(in);

    std::ofstream out(OUT.c_str());
    ForceAssert(out);

    String line;
    vec<String> lines;
    lines.reserve(N);
    int line_number = 0;
    while ( getline(in, line) )
    {
        if ( lines.isize() >= N )
        {
            if ( rand() % line_number < N )
            {
                int i = rand() % lines.size();
                lines[i] = line;
            }
        }
        else
        {
            lines.push_back(line);
        }
        line_number += 1;
    }
    ForceAssert(in.eof());

    printf("%d lines\n", lines.isize());
    fflush( stdout);
    ForceAssertEq(lines.isize(),N);

    for ( unsigned int i = 0; i < lines.size(); i++ )
    {
        out << lines[i] << '\n';
    }

    out.close();
    ForceAssert(out);
}

