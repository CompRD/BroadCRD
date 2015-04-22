// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// ReverseFasta: reverse sequences in a fasta file.  Inefficient.

#include <ctype.h>

#include "MainTools.h"
#include "dna/Bases.h"

void dumpBases( vec<char>& bases )
{
    int count = 0;
    vec<char>::reverse_iterator end(bases.rend());
    for ( vec<char>::reverse_iterator itr(bases.rbegin()); itr != end; ++itr )
    {
        cout << GeneralizedBase::complementChar(*itr);
        if ( ++count % 80 == 0 )
            cout << "\n";
    }
    if ( count % 80 != 0 )
        cout << "\n";
}

int main( int argc, char *argv[] )
{
    BeginCommandArguments;
    CommandArgument_String(PRE);
    CommandArgument_String(INPUT);
    EndCommandArguments;

    Ifstream( in, PRE + "/" + INPUT );

    vec<char> bases;
    bases.reserve(1000000);

    char c;
    while ( in >> c )
    {
        if ( c == '>' )
        {
            dumpBases(bases);
            bases.clear();
            String s;
            getline(in, s);
            cout << ">rc_" << s << "\n";
        }
        if ( isalpha(c) )
            bases.push_back(c);
    }
    dumpBases(bases);
}
