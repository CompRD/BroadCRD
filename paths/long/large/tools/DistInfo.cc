///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DistInfo.  Compute the distance in edges between two edges in a GapToy assembly.

#include "MainTools.h"
#include "paths/HyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", "looks for DIR/a.hbx");
     CommandArgument_String_Doc(E, "e1..e2");
     CommandArgument_Int_OrDefault_Doc(MAX_INT, 500, "maximum intermediate kmers");
     CommandArgument_Int_OrDefault_Doc(MAX_DIST, 30, "maximum distance to look");
     EndCommandArguments;

     int e1 = E.Before( ".." ).Int( ), e2 = E.After( ".." ).Int( );

     HyperBasevectorX hb;
     BinaryReader::readFile( DIR + "/a.hbx", &hb );

     vec<int> x = {e1};
     vec<int> d = {-1};
     vec<int> k = {0};

     double clock = WallClockTime( );
     for ( int j = 0; j < x.isize( ); j++ )
     {    int e = x[j];
          if ( e == e2 )
          {    cout << d[j] << endl;
               cout << TimeSince(clock) << " used" << endl;
               Scram(0);    }
          if ( k[j] > MAX_INT ) continue;
          if ( d[j] == MAX_DIST ) continue;
          int v = hb.ToRight(e);
          for ( int l = 0; l < (int) hb.From(v).size( ); l++ )
          {    int e = hb.IFrom( v, l );
               x.push_back(e);
               d.push_back( d[j] + 1 );
               k.push_back( k[j] + hb.Kmers(e) );    }    }
     cout << "Did not find e2." << endl;
     cout << TimeSince(clock) << " used" << endl;
     Scram(0);    }
