///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/// RaiseFasta: upper-case bases in a fasta file
#include <cctype>
#include <iostream>
#include <string>

// Please Note:  This is safe to use as a filter because it doesn't use any
// CRD code, which may, by design, unpredictably spew information to cout.
// Please don't mix any of our junk in.
// Also, please also consider using sed or awk instead of this: Even awk is
// probably faster than this, but I haven't checked.

int main()
{
    std::string line;
    while ( getline(std::cin,line) )
    {
        if ( line.size() && line[0] != '>' )
        {
            typedef std::string::iterator Itr;
            for ( Itr itr(line.begin()), end(line.end()); itr != end; ++itr )
                *itr = toupper(*itr);
        }
        std::cout << line << '\n';
    }
    return !std::cin.eof();
}
