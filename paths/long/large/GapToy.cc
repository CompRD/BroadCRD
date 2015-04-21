///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library JEMALLOC

#include "MainTools.h"
#include "paths/long/large/GapToyCore.h"

int main(int argc, char *argv[])
{
     RunTime( );
     GapToyCore( argc, argv );
     Scram(0);
}

#if 0
namespace std
{
    new_handler gCurHandler;
    new_handler set_new_handler( new_handler handler )
    { new_handler result = gCurHandler;
      gCurHandler = handler;
      return result; }
}

void* operator new( size_t siz )
{
    if ( !siz ) siz = 1ul;
    void* mem;
    while ( !(mem = malloc(siz)) )
    {
        if ( !std::gCurHandler )
            throw std::bad_alloc();
        (*std::gCurHandler)();
    }
    return mem;
}

void operator delete( void* mem )
{
    if ( mem ) free(mem);
}
#endif
