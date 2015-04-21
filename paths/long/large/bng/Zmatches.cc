///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/long/large/bng/Zmatches.h"

// Build alignment seeds.

void Zmatches( const vec<vec<double>>& X, const vec<int>& S,
     vec< quad<int,int,int,int> >& zmatches, const Bool FILTER_DUPS,
     const int div, const int w, const int cutoff )
{

#include "paths/long/large/bng/SuperScore.h"

     zmatches.clear( );

     ForceAssert( w == 3 || w == 4 );

     cout << Date( ) << ": building M" << endl;
     double eclock = WallClockTime( );

     int cutoff_sum = cutoff;

     vec<vec<int>> Z( X.size( ) );
     for ( int s = 0; s < S.isize( ); s++ )
     {    int i = S[s];
          if ( X[i].isize( ) >= w )
          {    Z[i].resize( X[i].isize( ) - w + 1, 0 );
               for ( int j = 0; j < Z[i].isize( ); j++ )
               for ( int k = 0; k < w; k++ )
                    Z[i][j] += X[i][j+k];    }    }
     vec< vec<int> > V;
     vec< pair<int,int> > P;
     for ( int64_t i = 0; i < Z.jsize( ); i++ )
     for ( int j = 0; j < Z[i].isize( ); j++ )
     {    if ( X[i].isize( ) >= w )
          {    vec<int> m;
               m.push_back( Z[i][j] / cutoff_sum );
               for ( int k = 0; k < w; k++ )
                    m.push_back( X[i][j+k] / cutoff );
               V.push_back(m);    
               P.push( i, j );    }    }
     cout << Date( ) << ": sorting" << endl;
     ParallelSortSync( V, P );
     int mcount = 0;
     cout << Date( ) << ": finding orbits" << endl;
     for ( int64_t i = 0; i < V.jsize( ); i++ )
     {    int64_t j = V.NextDiff(i);
          mcount++;
          i = j - 1;    }
     PRINT3( Z.size( ), V.size( ), mcount );

     vec<vec<int>> adds;
     if ( w == 3 )
     {    for ( int j1 = -1; j1 <= +1; j1++ )
          for ( int j2 = -1; j2 <= +1; j2++ )
          for ( int j3 = -1; j3 <= +1; j3++ )
          for ( int j4 = -1; j4 <= +1; j4++ )
          {    vec<int> x;
               x.push_back( j1, j2, j3, j4 );
               adds.push_back(x);    }    }
     if ( w == 4 )
     {    for ( int j1 = -1; j1 <= +1; j1++ )
          for ( int j2 = -1; j2 <= +1; j2++ )
          for ( int j3 = -1; j3 <= +1; j3++ )
          for ( int j4 = -1; j4 <= +1; j4++ )
          for ( int j5 = -1; j5 <= +1; j5++ )
          {    vec<int> x;
               x.push_back( j1, j2, j3, j4, j5 );
               adds.push_back(x);    }    }
     
     const int64_t batches = 500;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int bi = 0; bi < batches; bi++ )
     {    const int64_t u1 = ( V.size( ) * bi ) / batches;
          const int64_t u2 = ( V.size( ) * (bi+1) ) / batches;
          vec< quad<int,int,int,int> > zmatchesi;
          vec<int> x;
          for ( int i = u1; i < u2; i++ )
          {    int i1 = P[i].first, j1 = P[i].second;
               if ( div >= 0 && i1 < div ) continue;
               vec< quad<int,int,int,int> > zmatchesj;
               for ( int j = 0; j < adds.jsize( ); j++ )
               {    x = V[i];
                    for ( int l = 0; l < 5; l++ )
                         x[l] += adds[j][l];
                    int64_t low = LowerBound( V, x ), high;
                    if ( low == V.jsize( ) ) continue;
                    for ( high = low + 1; high < V.jsize( ); high++ )
                         if ( V[high] != V[low] ) break;    
                    for ( int64_t b = low; b < high; b++ )
                    {    int i2 = P[b].first, j2 = P[b].second;
                         if ( i1 == i2 ) continue;

                         if ( div >= 0 && i2 >= div ) continue;

                         // Skip duplicates, and don't align a read to its reverse.

                         if (FILTER_DUPS)
                         {    if ( i1/2 == i2/2 ) continue;
                              if ( i2 < i1 ) continue;
                              int ri1 = ( i1 % 2 == 0 ? i1 + 1 : i1 - 1 );
                              int ri2 = ( i2 % 2 == 0 ? i2 + 1 : i2 - 1 );
                              if ( make_pair( ri1, ri2 ) < make_pair( i1, i2 ) ) 
                                   continue;
                              if ( make_pair( ri2, ri1 ) < make_pair( i1, i2 ) ) 
                                   continue;    }

                         // Check proximity requirements.

                         if ( Abs( Z[i1][j1] - Z[i2][j2] ) > cutoff_sum ) continue;
                         Bool off = False;
                         for ( int p = 0; p < w; p++ )
                         {    if ( Abs( X[i1][j1+p] - X[i2][j2+p] ) > cutoff )
                              {    off = True;
                                   break;    }    }
                         if (off) continue;

                         zmatchesj.push( i1, i2, j1, j2 );    }    }
               UniqueSort(zmatchesj);

               vec<Bool> to_delete( zmatchesj.size( ), False );
               for ( int m = 0; m < zmatchesj.isize( ); m++ )
               {    int i1 = zmatchesj[m].first, i2 = zmatchesj[m].second;
                    int j1 = zmatchesj[m].third, j2 = zmatchesj[m].fourth;
                    double super = gsuper( i1, i2, j1, j2 );
                    if ( super < -10 ) to_delete[m] = True;    }
               EraseIf( zmatchesj, to_delete );

               zmatchesi.append(zmatchesj);    }
          #pragma omp critical
          {    zmatches.append(zmatchesi);    }    }
     cout << Date( ) << ": sorting matches" << endl;
     ParallelUniqueSort(zmatches);
     int64_t rcount = 0;
     for ( int64_t i = 0; i < zmatches.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < zmatches.jsize( ); j++ )
          {    if ( zmatches[j].first != zmatches[i].first
                    || zmatches[j].second != zmatches[i].second )
               {    break;    }    }
          rcount++;
          i = j - 1;    }
     cout << ToStringAddCommas( zmatches.size( ) ) << " matches ";
     cout << "between " << ToStringAddCommas(rcount) << " pairs of reads" << endl;
     cout << "ratio = " << setprecision(3) 
          << double( zmatches.size( ) ) / rcount << endl;
     cout << TimeSince(eclock) << " used\n" << endl;    }
