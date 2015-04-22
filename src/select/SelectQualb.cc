// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// SelectQualb: select some of the entries from a qualb file, yielding as output
// a new qualb file.

#include "system/ParsedArgs.h"
#include "Qualvector.h"
#include "system/RunTime.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"

int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(INPUT);
     CommandArgument_String(OUTPUT);
     CommandArgument_String_OrDefault(IDS, "");
     EndCommandArguments;

     String input = PRE + "/" + INPUT;

     vec<int> ids;
     Ifstream( idin, IDS );
     while(1)
     {    int n;
          idin >> n;
          if ( !idin ) break;
          ids.push_back(n);    }
     sort( ids.begin( ), ids.end( ) );

     vecqualvector Q;
     Q.Read( input, ids, 0 );
     Q.WriteAll( OUTPUT );    }
