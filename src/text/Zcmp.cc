// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// Zcmp: check if two compressed files have equal uncompressed contents.
// This is like the standard utility "zcmp", but does the computation on the fly,
// so that only trivial memory and disk space are needed.
//
// If DOT_PER specified, print a dot every DOT_PER bytes (to show progress).

#include "MainTools.h"
#include "system/ProcBuf.h"
#include <cstddef>
#include <iostream>

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FILE1);
     CommandArgument_String(FILE2);
     CommandArgument_UnsignedInt_OrDefault(DOT_PER, 0);
     EndCommandArguments;

     String command1 = "gzip -dc " + FILE1;
     procbuf pbuf1(command1.c_str(),std::ios_base::in);
     istream is1(&pbuf1);

     String command2 = "gzip -dc " + FILE2;
     procbuf pbuf2(command2.c_str(),std::ios_base::in);
     istream is2(&pbuf2);

     size_t count = 0;
     while ( is1 && is2 )
     {
         if ( is1.get() != is2.get() )
         {
             std::cout << "Not equal." << std::endl;
             return 1;
         }
         if ( DOT_PER && !(++count % DOT_PER) )
         {
             std::cout << '.';
             if ( !(count % (DOT_PER*80) ) )
                 std::cout << '\n';
             std::cout.flush();
         }
     }
     if ( !is1.eof() || !is2.eof() )
     {
         std::cout << "Not equal." << std::endl;
         return 1;
     }
     std::cout << "Equal." << std::endl;
     return 0;
}
