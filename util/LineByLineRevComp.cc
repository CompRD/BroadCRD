/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file LineByLineRevComp.cc
 * \author tsharpe
 * \date Jan 15, 2009
 *
 * \brief
 *
 *
 */
#include <cstring>
#include <iostream>
#include <string>

using namespace std;

char gRev[256];

int main( int argc, char** argv )
{
    memset(gRev,'?',sizeof(gRev));
    gRev['A'] = 'T';
    gRev['C'] = 'G';
    gRev['G'] = 'C';
    gRev['T'] = 'A';
    gRev['B'] = 'V';
    gRev['D'] = 'H';
    gRev['H'] = 'D';
    gRev['V'] = 'B';
    gRev['R'] = 'Y';
    gRev['Y'] = 'R';
    gRev['K'] = 'M';
    gRev['M'] = 'K';
    gRev['S'] = 'S';
    gRev['W'] = 'W';
    gRev['N'] = 'N';
    gRev['a'] = 't';
    gRev['c'] = 'g';
    gRev['g'] = 'c';
    gRev['t'] = 'a';
    gRev['b'] = 'v';
    gRev['d'] = 'h';
    gRev['h'] = 'd';
    gRev['v'] = 'b';
    gRev['r'] = 'y';
    gRev['y'] = 'r';
    gRev['k'] = 'm';
    gRev['m'] = 'k';
    gRev['s'] = 's';
    gRev['w'] = 'w';
    gRev['n'] = 'n';

    char buf[8192];
    while ( cin.good() )
    {
        streamsize len = cin.getline(buf,sizeof(buf)).gcount();
        if ( len )
        {
            string x(buf);
            string::reverse_iterator end = x.rend();
            for ( string::reverse_iterator itr = x.rbegin(); itr != end; ++itr )
            {
                cout << gRev[static_cast<unsigned char>(*itr)];
            }
            cout << '\n';
        }
    }
}
