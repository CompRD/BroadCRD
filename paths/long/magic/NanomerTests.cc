///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/long/magic/Nanomer.h"

int main( int argc, char* argv[] )
{
     RunTime( );

//     BeginCommandArguments;
//     CommandArgument_String(OUT_DIR);
//     EndCommandArguments;

     for (size_t i = 0; i < 1024; ++i )
         cout << "index " << i << ": " << Nanomer<5>(i).ToString() << endl;

     return 0;
}
