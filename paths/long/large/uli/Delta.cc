
// Uli mice -- two sample control.

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
     String dir = "/wga/scr4/jaffe/GapToy/mouse.a4";
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.final/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.final/a.inv", &inv );
     vec<vec<int>> count;
     BinaryReader::readFile( dir + "/a.final/a.countsb", &count );

     // Build statistics table.

     cout << Date( ) << ": building statistics" << endl;
     int ns = count.size( ); // number of samples
     vec< triple< int, int, vec<int> > > C;
     for ( int e = 0; e < hb.E( ); e++ )
     {    int nk = hb.Kmers(e), sum = 0;
          for ( int s = 0; s < ns; s++ )
               sum += count[s][e];
          vec<int> x(ns);
          for ( int s = 0; s < ns; s++ )
               x[s] = count[s][e];
          C.push( nk, sum, x );    }
     Sort(C);

     // Heuristics.

     const int max_look = 10;
     const double pceil = 0.001;
     const int min_high = 10;
     const double del = 0.04;
     const int passes = 15;
     const int min_peers = 1000;
     const int max_poly = 10;

     // Find edges having out of whack stats.

     cout << Date( ) << ": searching edges" << endl;
     cout << "\nThere are " << ToStringAddCommas( hb.E( ) ) << " edges." << endl;
     int calls = 0;
     int calls6 = 0;
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    int n = hb.Bases(e);
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( inv[e] < e ) continue;

          // Do preliminary tests of counts.

          int minc = 1000000000, maxc = 0;
          for ( int s = 0; s < ns; s++ )
          {    minc = Min( minc, count[s][e] );
               maxc = Max( maxc, count[s][e] );    }
          if ( minc > max_look || maxc < min_high ) continue;

          // Exclude edges starting or ending with long homopolymers.

          const basevector& E = hb.O(e);
          int mp = 0;
          for ( int j = 1; j < E.isize( ); j++ )
          {    if ( E[j] != E[0] ) break;
               mp = Max( mp, j + 1 );    }
          for ( int j = E.isize( ) - 2; j >= 0; j-- )
          {    if ( E[j] != E[ E.isize( ) - 1 ] ) break;
               mp = Max( mp, E.isize( ) - j );    }
          if ( mp > max_poly ) continue;

          // Proceed.

          int nk = hb.Kmers(e), sum = 0;
          for ( int s = 0; s < ns; s++ )
               sum += count[s][e];

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

          for ( int s = 0; s < ns; s++ )
          {    int c = count[s][e];
               if ( c > max_look ) continue;

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
               {    
                    #pragma omp critical
                    {    cout << "\n";
                         PRINT8( e, nk, sum, goods, s, c, mean, p );    
                         calls++;
                         if ( p < 0.000001 ) calls6++;    }    }    }    }

     // Print summary stats.

     PRINT(calls);
     PRINT(calls6);
     cout << TimeSince(clock) << " used" << endl;    }
