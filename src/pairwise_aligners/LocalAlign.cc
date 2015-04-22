// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// LocalAlign( S, T, a )
//
// Find the best local alignment of basevector "S" and "T" and put it in align "a".
//
// This is computed relative to fixed constant scores for "match", "mismatch",
// and "gap".
//
// If the best alignment has length 0, the alignment will have zero blocks, but
// be otherwise meaningless.

#include "Basevector.h"
#include "math/Functions.h"
#include "pairwise_aligners/LocalAlign.h"
#include "PackAlign.h"
#include "ShortVector.h"

int LocalAlign( const basevector& S, const basevector& T, align& a ,
		const int match, const int mismatch, const int gap)
{    
     int n = S.size( ), N = T.size( );

     avector<char> s, t;
     s.resize(n);
     for ( int i = 0; i < n; i++ )
          s(i) = S[i];
     t.resize(N);
     for ( int i = 0; i < N; i++ )
          t(i) = T[i];

     vec< vec<int> > f(n+1);
     for ( int i = 0; i <= n; i++ )
          f[i].resize_and_set(N+1, 0);

     vec< vec<unsigned char> > from( n+1, vec<unsigned char>(N+1) );

     int Mf = 0, best_i = 0, best_j = 0;
     for ( int i = 1; i <= n; i++ )
     {    int *fim = &f[i-1][0], *fi = &f[i][0];
          unsigned char *fri = &from[i][0];
          int a, b, c, d = 0;
          for ( int j = 1; j <= N; j++ )
          {    a = fim[j-1] + ( s(i-1) == t(j-1) ? match : mismatch );
               b = fim[j] + gap;
               c = fi[j-1] + gap;
               fi[j] = std::max( {a, b, c, d} );
               if ( a >= b && a >= c && a >= d ) fri[j] = 'a';
               else if ( b >= a && b >= c && b >= d ) fri[j] = 'b';
               else if ( c >= a && c >= b && c >= d ) fri[j] = 'c';
               else if ( d >= a && d >= b && d >= c ) fri[j] = 'd';
               if ( fi[j] > Mf ) 
               {    Mf = fi[j];
                    best_i = i;
                    best_j = j;    }    }    }

     int i = best_i, j = best_j;
     vec< pair<int,int> > matches;
     while(1)
     {    if ( i == 0 || j == 0 ) break;
          if ( from[i][j] == 'a' )
          {    --i;
               --j;    
               matches.push_back( make_pair(i,j) );    }
          else if ( from[i][j] == 'b' ) --i;
          else if ( from[i][j] == 'c' ) --j;
          else break;    }
     matches.ReverseMe( );

     if ( matches.empty( ) )
     {    a.SetNblocks(0);
          return 0;    }

     int pos1 = matches[0].first, pos2 = matches[0].second;
     int p1 = pos1, p2 = pos2;
     avector<int> gaps(0), lengths(0);
     gaps.Append(0);
     for ( int k = 0; k < matches.isize( ); k++ )
     {    int l;
          for ( l = k + 1; l < matches.isize( ); l++ )
          {    if ( matches[l].first != matches[l-1].first + 1
                    || matches[l].second != matches[l-1].second + 1 )
               {    break;    }    }
          if ( k > 0 )
          {    int g = ( matches[k].second - p2 ) - ( matches[k].first - p1 );
               gaps.Append(g);    }
          lengths.Append(l-k);
          p1 = matches[k].first + l-k;
          p2 = matches[k].second + l-k;
          k = l - 1;    }
     a.Set( pos1, pos2, gaps, lengths );
     return Mf;    }
