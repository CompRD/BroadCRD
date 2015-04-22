// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// Substitute: scan each input line for "FROM".  If it contains it, replace
// it by "TO".

#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FROM);
     CommandArgument_String(TO);
     CommandArgument_Bool_OrDefault(ALL, False);
     EndCommandArguments;

     TO.GlobalReplaceBy( "\\n", "\n" );
     String line;
     while(1)
     {    getline( cin, line );
          if ( !cin ) break;
          if ( line.Contains(FROM) ) 
          {    if ( !ALL ) line.ReplaceBy( FROM, TO );
               else line.GlobalReplaceBy( FROM, TO );    }
          cout << line << "\n";    }    }
     
