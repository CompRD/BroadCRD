///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file CleanFasta.cc
 * \author tsharpe
 * \date Feb 16, 2012
 *
 * \brief
 */

#include "util/CleanFasta.h"
#include "system/Assert.h"
#include <fstream>
#include <string>

void CleanFasta( String const& fastaIn, String const& fastaOut )
{
    std::ifstream in(fastaIn.c_str());
    std::ofstream out(fastaOut.c_str(),std::ios_base::out|std::ios_base::trunc);
    std::string line;
    while ( std::getline(in,line) )
    {
        if ( line.size() && line[0] == '>' )
            out << line << '\n';
        else
        {
            typedef std::string::iterator Itr;
            for ( Itr itr(line.begin()), end(line.end()); itr != end; ++itr )
            {
                char ccc = *itr;
                if ( ccc != 'n' && ccc != 'N' )
                    out << ccc;
            }
            out << '\n';
        }
    }
    out.close();
    ForceAssert(in.eof());
    ForceAssert(out);
}
