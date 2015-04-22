///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Remove duplicates from magic.

#include "Basevector.h"
#include "MainTools.h"

int main( )
{    vecbasevector bases( "/wga/scr4/macro/sample.ecoli.Feb6.fastb" );
     vec<Bool> to_delete( bases.size( ), False );
     for ( int i = 1; i < (int) bases.size( ); i++ )
          if ( bases[i] == bases[i-1] ) to_delete[i] = True;
     bases.EraseIf(to_delete);
     bases.WriteAll( "/wga/scr4/macro/sample.ecoli.Feb6.dedup.fastb" );    }
