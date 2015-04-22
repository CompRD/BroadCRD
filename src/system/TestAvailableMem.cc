///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Test memory to see if it's really available.  It seems that this has to be
// a separate executable: you can't just copy the code to inside another program
// without messing up its memory manager.

#include <malloc.h>
#include "MainTools.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArgumentsAcceptEmptyArgListNoHeader;
     CommandArgument_Double_OrDefault_Doc( MAX_MEM_GB, 0,
          "if specified, maximum allowed RAM use in GB; in some cases may be "
          "exceeded by our code");
     EndCommandArguments;

     if ( MAX_MEM_GB > 0 )
     {    int64_t max_bytes
               = int64_t( round( MAX_MEM_GB * 1024.0 * 1024.0 * 1024.0 ) );
          SetMaxMemory(max_bytes);    }
     int64_t P = GetMaxMemory( );
     int64_t d;
     for ( d = 1; d <= 100; d++ )
     {    int64_t M = (d*P)/100;
          void* ptr = malloc(M);
          if ( ptr == NULL )
          {    d--;
               cout << "\nAble to allocate " << d << "% of "
                    << "supposedly available memory (and not more)." << endl;
               cout << "Can access " << M / (1024*1024*1024) << " GB." << endl;
               if ( d < 90 )
               {   cout << "WARNING: you may have less available memory than "
                        << "you think!" << endl;    }
               else 
               {    cout << "But you're pretty close, so things should be OK."
                         << endl;    }
               M = (d*P)/100;
               free(ptr);
               break;    }
          free(ptr);    }
     if ( d > 100 )
     {    cout << "\nAble to allocate 100% of supposedly available memory."
               << endl;
          cout << "Can access at least " << P / (1024*1024*1024) << " GB." 
               << endl;    }
     cout << endl;    }
