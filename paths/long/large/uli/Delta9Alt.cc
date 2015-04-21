///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Uli mice - find the mutations.

// order  assembly  true   state   sample
// 0      mouse1    #1     0       Balb/c, male
// 1      mouse4    #4b    1       son, responder
// 2      mouse5    #5     1       son, responder
// 3      mouse6    #6     0       son, nonresponder
// 4      mousea    #4a    1       son, responder
// 5      mouseb    #7     0       son, nonresponder
// 6      mousec    #3     1       ~B6 mother, responder
// 7      moused    #2     0       B6 father
// 8      mousee    DUP    1       reads reordered from mousea

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"

int main( )
{    RunTime( );
     double clock = WallClockTime( );
     cout << Date( ) << ": loading" << endl;
     String dir = "/wga/scr4/jaffe/GapToy/mouse.1456abcde";
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.final/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.final/a.inv", &inv );
     vec<vec<int>> count;
     BinaryReader::readFile( dir + "/a.final/a.countsb", &count );

     // Add in duplicated reads.

     for ( int e = 0; e < hb.E( ); e++ )
          count[1][e] += count[4][e];

     // Heuristics.

     const int max_look = 10;
     const double pceil = 0.02;
     const int min_high = 7;
     const double del = 0.04;
     const int passes = 15;
     const int min_peers = 1000;
     const int max_poly = 10;

     vec<int> resp = {1,2,6};
     vec<int> nonresp = {0,3,5,7};

     // Find edges having out of whack stats.

     cout << Date( ) << ": searching edges" << endl;
     cout << "\nThere are " << ToStringAddCommas( hb.E( ) ) << " edges." << endl;
     vec< triple<double,int,int> > calls;
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    int n = hb.Bases(e);
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( inv[e] < e ) continue;

          int m = 1000000000;
          for ( int j = 0; j < resp.isize( ); j++ )
               m = Min( m, count[ resp[j] ][e] );
          if ( m < min_high ) continue;

          int rsum = 0;
          for ( int j = 0; j < resp.isize( ); j++ )
               rsum += count[ resp[j] ][e];

          int nsum = 0;
          for ( int j = 0; j < nonresp.isize( ); j++ )
               nsum += count[ nonresp[j] ][e];

          if ( rsum < 20 * nsum ) continue;

          // if ( nsum > 2 ) continue;

          #pragma omp critical
          {    calls.push( m, nsum, e );    }    }

     // Print summary stats.

     Sort(calls);
     for ( int i = 0; i < calls.isize( ); i++ )
     {    int m = calls[i].first;
          int sum = calls[i].second;
          int e = calls[i].third;    
          cout << "[" << i+1 << "] ";
          cout << "e = " << e << ", inv[e] = " << inv[e] << ", m = " << m
               << ", sum = " << sum << ", resp =";
          for ( int j = 0; j < resp.isize( ); j++ )
               cout << " " << count[ resp[j] ][e];
          cout << ", nonresp =";
          for ( int j = 0; j < nonresp.isize( ); j++ )
               cout << " " << count[ nonresp[j] ][e];
          cout << endl;    }    }
