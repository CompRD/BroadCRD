///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PrintAmb: print out the ambiguities in an efasta file F.

#include "MainTools.h"
#include "efasta/EfastaTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(F);
     EndCommandArguments;

     VecEFasta x;
     LoadEfastaIntoStrings( F, x );
     for ( size_t i = 0; i != x.size( ); i++ )
     {    cout << "\nrecord " << i << endl;
          x[i].PrintAmbiguities( cout, ToString(i) );    }    }
