///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Measure the locality of a HyperBasevector by computing the median difference
// between the ids of two adjacent vertices, and between the ids of two adjacent 
// edges.

#include "MainTools.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(HBV, "HyperBasevector assembly");
     EndCommandArguments;

     HyperBasevector hb;
     BinaryReader::readFile( HBV, &hb );

     vec<int> vdist, edist;
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int w = hb.From(v)[j];
               int e = hb.EdgeObjectIndexByIndexFrom( v, j );
               vdist.push_back( Abs( v - w ) );
               for ( int k = 0; k < hb.From(w).isize( ); k++ )
               {    int f = hb.EdgeObjectIndexByIndexFrom( w, k );
                    edist.push_back( Abs( e - f ) );    }    }    }
     Sort(vdist), Sort(edist);
     cout << "median vertex id difference = " << Median(vdist) << endl;
     cout << "median edge id difference = " << Median(edist) << endl;    }
