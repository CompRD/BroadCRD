///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FastbSize: report the size of a fastb record.

#include "Basevector.h"
#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FILE);
     CommandArgument_Int(ID);
     EndCommandArguments;

     if ( !FILE.Contains(".fastb") )
          cout << "Warning: " << FILE << " does not end in fastb." << endl;
     if ( !IsRegularFile(FILE) ) FatalErr( "File " << FILE << " doesn't exist" );

     vecbasevector b;
     b.ReadOne( FILE, ID );
     cout << b[0].size( ) << endl;    }
