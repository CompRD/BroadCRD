///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Convert PerfStat output into standard form for dexter to read.
// Goes from standard input to standard output.

#include "MainTools.h"

int main( )
{
     RunTime( );

     String line, long_name, short_name, value;
     while(1)
     {    getline( cin, line );
          if ( cin.fail( ) ) break;
          if ( !line.Contains( "PERFSTAT: ", 0 ) ) continue;
          long_name = line.Between( "PERFSTAT: ", "[" );
          short_name = line.Between( "[", "]" );
          value = line.After( "] = " );
          cout << short_name << '\t' << value << '\t' 
               << long_name << '\n';    }    }
     
