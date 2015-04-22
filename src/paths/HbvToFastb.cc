///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Read a HyperBasevector and write a fastb file of its edges.

#include "MainTools.h"
#include "paths/HyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     EndCommandArguments;

     HyperBasevector hb;
     BinaryReader::readFile( IN, &hb );
     vecbasevector edges( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          edges[e] = hb.EdgeObject(e);
     edges.WriteAll(OUT);    }

