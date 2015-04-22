// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// GrepSmallFiles: check all files in directory D smaller than NBYTES to see 
// if they contain a particular pattern X.

#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(D);
     CommandArgument_UnsignedInt(NBYTES);
     CommandArgument_String(X);
     EndCommandArguments;

     vec<String> all = AllFiles(D);
     String quote = "\"";
     String com = "cd " + D + "; grep " + quote + X + quote;
     for ( int i = 0; i < all.isize( ); i++ )
     {    String fn = D + "/" + all[i];
          if ( IsRegularFile(fn) && FileSize(fn) < NBYTES ) 
               com += " " + quote + all[i] + quote;    }
     System(com);    }
