///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/long/large/bng/BngAlign.h"

// Align a BNG map R to a DISCOVAR assembly line A.  We assume that position r_m on
// R maps to position a_m on A.  This uses a shortest path algorithm.  Diagnostic
// output goes to "out".
//
// ndirect = number of positions in "direct" correspondence in the alignment
// see PRINT_FANCY to understand what this means
// not entirely clear what it means

double BngAlign( const vec<int>& R, const vec<int>& A, const int r_m, const int a_m,
     ostringstream& out, vec<int>& p,
     int& rstart, int& rstop, int& ndirect, int& rsum, int& asum, const int mode )
{
     p.clear( );
     int nr = R.size( ), na = A.size( );
     int N = 2 + nr*na; // number of vertices
     PRINT3_TO( out, nr, na, N );

     // delta(x,y) = penalty for reading x for y              = |x-y|
     // skip(y)    = penalty for skipping segment of length y = y/3
     //
     // VERTEX                                     INDEX
     // 1. begin, end                              0, 1
     // 2. r_i x a_j                               2 + (i * na) + j
     //
     // EDGE                                       WEIGHT 
     // 1. begin --> r_0 x a_j, j < a_m            0
     //    begin --> r_i x a_0, i < r_m            0
     // 2. r_nr-1 x a_j --> end, j > a_m           delta( r_nr-1, a_j )
     //    r_i x a_na-1 --> end, i > r_m           delta( r_i, a_na-1 )
     //                          *******           (note seems not to be enforced)
     // 3. r_i x a_j --> r_i+1 x a_j+1             delta( r_i, a_j )
     // 4. r_i x a_j --> r_i+1 x a_j+2             delta( r_i, a_j + a_j+1 )
     //                                            + skip( Min( a_j, a_j+1 ) )
     // and similar to 4 for skips of two and three elements.
     // 5. r_i x a_j --> r_i+2 x a_j+1             delta( r_i + r_i+1, a_j )
     //                                            + skip( Min( r_i, r_i+1 ) )
     // and similar for skip of two elements.

     auto pos = [=]( int i, int j ){ return 2 + (i * na) + j; };
     auto r_ind = [=]( int pos ){ return (pos - 2) / na; };
     auto a_ind = [=]( int pos ){ return (pos - 2) % na; };

     auto delta = [=]( int x, int y ){ return Abs(x-y); };
     auto skip = [=]( int x ){ return x / 3.0; };

     auto vert_name = [=]( int v )
     {    if ( v == 0 ) return String("begin");
          if ( v == 1 ) return String("end");
          int i = r_ind(v);
          int j = v - 2 - (i*na);
          return "r" + ToString(i) + ".a" + ToString(j);    };

     // Set up a function f that defines weighted edges in a graph G.

     auto f = [&]( const int v )
     {    vec< pair<int,double> > x;
          if ( v == 0 ) // edge type 1
          {    for ( int j = 0; j < a_m; j++ )
                    x.push( pos(0, j), 0 );
               for ( int i = 0; i < r_m; i++ )
                    x.push( pos(i, 0), 0 );    }
          else if ( v >= pos( nr-1, a_m ) ) // edge type 2
          {    int j = v - pos( nr-1, 0 );
               x.push( 1, delta( R[nr-1], A[j] ) );    }
          else if ( v > 1 && (v-1) % na == 0 ) // also edge type 2
          {    int i = (v-1)/na - 1;
               if ( mode == 1 || i > r_m )
               x.push( 1, delta( R[i], A[na-1] ) );    }
          else // edge types 3, 4, 5;
          {    int i = r_ind(v), j = a_ind(v);
               if ( i < nr - 1 )
               {    if ( j < na - 1 ) x.push( pos(i+1, j+1), delta( R[i], A[j] ) );

                    for ( int m = 2; m <= Min( 4, na - j /* - 1 */ ); m++ )
                    {    
                         // Let sk be the sum of skip applied to A[j],...,A[j+m-1],
                         // but excluding the largest element, and let sum be the
                         // sum of all the elements.
     
                         int M = -1, Mi = -1;
                         for ( int l = 0; l < m; l++ )
                         {    if ( A[j+l] > M )
                              {    M = A[j+l];
                                   Mi = l;    }    }
                         double sk = 0, sum = 0;
                         for ( int l = 0; l < m; l++ )
                         {    sum += A[j+l];
                              if ( l != Mi ) sk += skip( A[j+l] );    }
     
                         // Create edge.

                         int target = pos( i+1, j+m );
                         if ( j+m == A.isize( ) ) 
                         {    if ( j >= a_m ) target = 1;
                              else continue;    }
                         x.push( target, delta( R[i], sum ) + sk );    }    }

               if ( i < nr - 2 && j < na - 1 )
               {    x.push( pos(i+2, j+1), delta( R[i] + R[i+1], A[j] )
                         + skip( Min( R[i], R[i+1] ) ) );    }    

               if ( i < nr - 3 && j < na - 1 )
               {    x.push( pos(i+3, j+1), delta( R[i] + R[i+1] + R[i+2], A[j] )
                         + skip( Min( R[i], R[i+1] ) )
                         + skip( Min( R[i+1], R[i+2] ) )
                              );    }    

                         }

          return x;    };

     // Now define the graph G.

     digraphE_V1<double> G( N, f );

     // Find shortest path from begin to end.

     double d = G.ShortestPath( 0, 1, p );
     double div = 0;
     if ( p.isize( ) < 3 )
     {    out << "|p| = " << p.size( ) << endl;
          rstart = -1;
          rstop = -1;
          return 1000000000;    }
     // ForceAssertGe( p.isize( ), 3 );

     int xa1 = ( p[1] - 2 ) % na;
     int xa2 = ( p[ p.isize( ) - 2 ] - 2 ) % na;

     int xr1 = ( p[1] - 2 ) / na;
     int xr2 = ( p[ p.isize( ) - 2 ] - 2 ) / na;

     for ( int j = xa1; j <= xa2; j++ )
          div += A[j];
     d /= div;

     // Print result.

     Bool PRINT_SIMPLE = False;
     Bool PRINT_FANCY = False;
     if ( mode == 2 ) PRINT_FANCY = True;
     if (PRINT_SIMPLE)
     {    for ( int j = 0; j < p.isize( ); j++ )
          {    if ( j > 0 ) out << "--> ";
               out << vert_name( p[j] ) << "\n";    }    }
     if (PRINT_FANCY)
     {    for ( int j = 0; j < p.isize( ); j++ )
          {    if ( j == 0 ) out << "begin\n";
               else
               {    out << "--> ";
                    int k;
                    for ( k = j + 1; k < p.isize( ); k++ )
                         if ( p[k] - p[k-1] != na + 1 ) break;
                    if ( k-1 > j )
                    {    out << vert_name( p[j] ) << " ... " 
                              << vert_name( p[k-1] ) << "\n";    }
                    else out << vert_name( p[j] ) << "\n";
                    j = k - 1;    }    }
          out << "d = " << setprecision(3) << 100*d << "%" << endl;
          rsum = 0, asum = 0;
          out << "R: ";
          for ( int l = xr1; l <= xr2; l++ )
          {    if ( l > xr1 ) out << ",";
               out << R[l];
               rsum += R[l];    }
          out << " (sum=" << rsum << ")\n";
          out << "A: ";
          for ( int l = xa1; l <= xa2; l++ )
          {    if ( l > xa1 ) out << ",";
               out << A[l];
               asum += A[l];    }
          out << " (sum=" << asum << ")\n";    }

     // Compute ndirect.

     ndirect = 0;
     for ( int j = 1; j < p.isize( ); j++ )
     {    int k;
          for ( k = j + 1; k < p.isize( ); k++ )
               if ( p[k] - p[k-1] != na + 1 ) break;
          if ( k-1 > j ) ndirect += k - j;
          j = k - 1;    }

     // Fix p.

     vec<int> q;
     for ( int j = 1; j < p.isize( ) - 1; j++ )
          q.push_back( p[j] - 2 );
     p = q;

     // Set return values.

     rstart = xr1;
     rstop = xr2;
     return d;    }
