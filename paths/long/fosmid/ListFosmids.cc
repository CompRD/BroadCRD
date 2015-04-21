///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ListFosmids.  List Fosmids meeting certain criteria.

#include "MainTools.h"
#include "ParseSet.h"
#include "paths/long/fosmid/Fosmids.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(CLASS, "goods or oks or goods,oks");
     EndCommandArguments;

     vec<int> answer = expand_fosmids(CLASS);
     for ( int i = 0; i < answer.isize( ); i++ )
     {    if ( i > 0 ) cout << " ";
          cout << answer[i];    }
     cout << "\n";    }
     
