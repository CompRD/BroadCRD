// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// SplitFastaByAmbiguity: split a fasta file into two parts: those records
// which contain ambiguous bases, and those which do not.

#include "MainTools.h"
#include "FastIfstream.h"
#include "dna/Bases.h"

int main( int argc, char *argv[] )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String(IN);
    CommandArgument_String(OUT_CLEAN);
    CommandArgument_String(OUT_DIRTY);
    EndCommandArguments;

    vec<String> record(100000);
    fast_ifstream in(IN);
    Ofstream( clean, OUT_CLEAN );
    Ofstream( dirty, OUT_DIRTY );
    int count = 0;
    Bool amb = False;
    String line;
    while ( 1 )
    {
        getline(in, line);
        if ( in.fail() || line.Contains(">", 0) )
        {
            for ( int i = 0; i < count; i++ )
            {
                if ( !amb )
                    clean << record[i] << "\n";
                else
                    dirty << record[i] << "\n";
            }
            if ( in.fail() )
                break;
            record[0] = line;
            count = 1;
            amb = False;
        }
        else
        {
            if ( !amb )
                amb = !Base::areBases(line.begin(),line.end());
            if ( count == 100000 )
                FatalErr( "Fasta record too large." );
            record[count++] = line;
        }
    }
}
