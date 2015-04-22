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

     // Build statistics table.

     cout << Date( ) << ": building statistics" << endl;
     int ns = count.size( ); // number of samples
     vec< triple< int, int, vec<int> > > C;
     for ( int e = 0; e < hb.E( ); e++ )
     {    int nk = hb.Kmers(e), sum = 0;
          for ( int s = 0; s < ns; s++ )
          {    if ( s == 4 || s == 8 ) continue;
               sum += count[s][e];    }
          vec<int> x(ns);
          for ( int s = 0; s < ns; s++ )
          {    if ( s == 4 || s == 8 ) continue;
               x[s] = count[s][e];    }
          C.push( nk, sum, x );    }
     Sort(C);

     // Heuristics.

     const int max_look = 10;
     const double pceil = 0.02;
     const int min_high = 5;
     const double del = 0.04;
     const int passes = 15;
     const int min_peers = 1000;
     const int max_poly = 10;

     vec<int> resp = {1,2,6};
     vec<int> nonresp = {0,3,5,7};

     // swap( resp, nonresp );

     // Find edges having out of whack stats.

     cout << Date( ) << ": searching edges" << endl;
     cout << "\nThere are " << ToStringAddCommas( hb.E( ) ) << " edges." << endl;
     vec< triple<double,int,int> > calls;
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    int n = hb.Bases(e);
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( inv[e] < e ) continue;

          // Do preliminary tests of counts.

          Bool OK = True;
          for ( int j = 0; j < resp.isize( ); j++ )
               if ( count[ resp[j] ][e] < min_high ) OK = False;
          for ( int j = 0; j < nonresp.isize( ); j++ )
               if ( count[ nonresp[j] ][e] > max_look ) OK = False;
          if ( !OK ) continue;

          int minc = 1000000000, maxc = 0;
          for ( int s = 0; s < ns; s++ )
          {    if ( s == 4 || s == 8 ) continue;
               minc = Min( minc, count[s][e] );
               maxc = Max( maxc, count[s][e] );    }
          if ( minc > max_look || maxc < min_high ) continue;

          // Exclude edges starting or ending with long homopolymers.
          // Commented out.

          const basevector& E = hb.O(e);
          int mp = 0;
          for ( int j = 1; j < E.isize( ); j++ )
          {    if ( E[j] != E[0] ) break;
               mp = Max( mp, j + 1 );    }
          for ( int j = E.isize( ) - 2; j >= 0; j-- )
          {    if ( E[j] != E[ E.isize( ) - 1 ] ) break;
               mp = Max( mp, E.isize( ) - j );    }
          // if ( mp > max_poly ) continue;

          // Proceed.

          int nk = hb.Kmers(e), sum = 0;
          for ( int s = 0; s < ns; s++ )
          {    if ( s == 4 || s == 8 ) continue;
               sum += count[s][e];    }

          // Find peer group for (nk,sum).

          int64_t low, high;
          int goods, sum1, sum2;
          for ( int pass = 1; pass <= passes; pass++ )
          {    double lm = 1.0 - ( del * pass ), hm = 1.0 + ( del * pass );
               int nk1 = int( floor( nk * lm ) ), nk2 = int( ceil( nk * hm ) );
               sum1 = int( floor( sum * lm ) ), sum2 = int( ceil ( sum * hm ) );
               low = LowerBound1( C, nk1 ), high = UpperBound1( C, nk2 );
               goods = 0;
               for ( int64_t i = low; i < high; i++ )
                    if ( C[i].second >= sum1 && C[i].second <= sum2 ) goods++; 
               if ( goods >= min_peers ) break;    }

          // Test each sample.

          Bool pass = True;
          // ostringstream out;
          double max_p = 0;
          int min_pos = 1000000000;
          for ( int s = 0; s < ns; s++ )
          {    if ( s == 4 || s == 8 ) continue;
               if ( !Member( nonresp, s ) ) 
                    min_pos = Min( min_pos, count[s][e] );    }
          for ( int s = 0; s < ns; s++ )
          {    if ( !Member( nonresp, s ) ) continue;
               int c = count[s][e];

               // Get mean for peer group.

               double mean = 0;
               for ( int64_t i = low; i < high; i++ )
               {    if ( C[i].second >= sum1 && C[i].second <= sum2 )
                         mean += C[i].third[s];    }
               mean /= goods;

               // Given a Poisson random variable with lambda = mean, what's
               // the probability of seeing <= c events?

               double p = PoissonCdfLong( mean, c );

               // Test for low probability event.

               if ( p <= pceil )
               {    max_p = Max( max_p, p );
                    // out << "\n";
                    // PRINT8_TO( out, e, nk, sum, goods, s, c, mean, p );    
                         }
               else 
               {    pass = False;
                    break;    }    }
          if (pass)
          {
               #pragma omp critical
               {    calls.push( max_p, min_pos, e );    }    }    }

     // Print summary stats.

     Sort(calls);
     for ( int i = 0; i < calls.isize( ); i++ )
     {    double max_p = calls[i].first;
          int min_pos = calls[i].second;
          int e = calls[i].third;    
          cout << "[" << i+1 << "] ";
          cout << "e = " << e << ", inv[e] = " << inv[e] << ", ";
          PRINT2( max_p, min_pos );    }    }

