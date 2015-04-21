///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Count paths e1..e2 with length at most L in the middle.

#include "MainTools.h"
#include "paths/HyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(e1);
     CommandArgument_Int(e2);
     CommandArgument_Int(L);
     EndCommandArguments;

     // Define directories.

     String dir = "/wga/scr4/jaffe/GapToy/51400.newchem/a.final";

     cout << "\n" << Date( ) << ": loading assembly" << endl;
     HyperBasevector hb;
     BinaryReader::readFile( dir + "/a.hbv", &hb );

     vec<int> len;
     for ( int i = 0; i < hb.E( ); i++ )
          len.push_back( hb.Kmers(i) );
     digraphE<int> G( hb, len );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     int v = to_right[e1], w = to_left[e2];

     vec<vec<int>> paths;
     double clock = WallClockTime( );
     cout << Date( ) << ": start" << endl;
     G.AllPathsLengthRange( v, w, 0, L, to_right, paths );
     cout << Date( ) << ": stop" << endl;
     cout << "found " << paths.size( ) << " paths, time used = "
          << TimeSince(clock) << endl;
     Scram(0);    }
