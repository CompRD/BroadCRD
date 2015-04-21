///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/LongReadTools.h"
#include "paths/MapReadsToHyper.h"

void MapReadsToHyper( const vecbasevector& bases, const vecqualvector& quals, 
     const HyperBasevector& hb, vec< vec< vec<int> > >& rpaths, 
     vec< vec< pair<int,int> > >& aligns, vec< vec<int> >& sum, const int minq )
{
     // Hash the edges.

     vecbasevector edges;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          edges.push_back( hb.EdgeObject(i) );
     const int L = 12;
     vec< vec< pair<int,int> > > Ulocs( IPow( 4, L ) );
     for ( size_t i = 0; i < edges.size( ); i++ )
     {    for ( int j = 0; j <= edges[i].isize( ) - L; j++ )
          {    int n = KmerId( edges[i], L, j );
               Ulocs[n].push( i, j );    }    }

     // Align the reads to the graph.  Seed on 12-mers.  Do not allow indels.

     rpaths.clear( ), rpaths.resize( bases.size( ) );
     aligns.clear( ), aligns.resize( bases.size( ) );
     sum.clear( ), sum.resize( bases.size( ) );
     const int Qthresh = 20;
     #pragma omp parallel for
     for ( size_t i = 0; i < bases.size( ); i++ )
     {    const basevector& b = bases[i];
          const qualvector& q = quals[i];
          for ( int rpos = 0; rpos <= b.isize( ) - L; rpos++ )
          {    int n = KmerId( b, L, rpos );
               for ( int j = 0; j < Ulocs[n].isize( ); j++ )
               {    int u = Ulocs[n][j].first;
                    int upos = Ulocs[n][j].second;
                    int offset = rpos - upos;
                    aligns[i].push( u, offset );    }    }
          UniqueSort(aligns[i]);
          int naligns = aligns[i].size( );

          // Chain the aligns.
               
          const int K = hb.K( );
          vec<int> to_left, to_right;
          hb.ToLeft(to_left), hb.ToRight(to_right);
          vec< vec<int> > from(naligns), to(naligns);
          for ( int j1 = 0; j1 < naligns; j1++ )
          {    for ( int j2 = 0; j2 < naligns; j2++ )
               {    int u1 = aligns[i][j1].first, u2 = aligns[i][j2].first;
                    if ( to_right[u1] != to_left[u2] ) continue;
                    int o1 = aligns[i][j1].second, o2 = aligns[i][j2].second;
                    if ( o1 + edges[u1].isize( ) - o2 != K - 1 ) continue;
                    from[j1].push_back(j2);
                    to[j2].push_back(j1);    }    }    
          for ( int z = 0; z < naligns; z++ )
          {    Sort( from[z] );
               Sort( to[z] );    }
          digraph G( from, to );

          // Find all paths through the graph.  For now we ignore paths that
          // don't extend the full length of the read.

          G.AllPaths( -1, -1, rpaths[i] );
          vec<Bool> to_delete( rpaths[i].size( ), False );
          for ( int j = 0; j < rpaths[i].isize( ); j++ )
          {    const vec<int>& p = rpaths[i][j];
               const pair<int,int>& a1 = aligns[i][ p.front( ) ];
               const pair<int,int>& a2 = aligns[i][ p.back( ) ];
               if ( a1.second > 0 ) to_delete[j] = True;
               if ( a2.second + edges[a2.first].isize( ) < b.isize( ) )
                    to_delete[j] = True;    }
          EraseIf( rpaths[i], to_delete );

          // Remove superfluous unipaths from beginnings and ends of rpaths.

          for ( int j = 0; j < rpaths[i].isize( ); j++ )
          {    vec<int>& p = rpaths[i][j];
               int l1, l2;
               for ( l1 = 0; l1 < p.isize( ); l1++ )
                    if ( aligns[i][ p[l1] ].second > 0 ) break;
               for ( l2 = p.isize( ) - 1; l2 >= 0; l2-- )
               {    if ( aligns[i][ p[l2] ].second 
                         + edges[ aligns[i][ p[l2] ].first ].isize( ) 
                         < b.isize( ) )
                    {    break;    }    }
               vec<int> p2;
               for ( int m = l1 - 1; m <= l2 + 1; m++ )
                    p2.push_back( p[m] );
               p = p2;    }
          UniqueSort(rpaths[i]);

          // Score each path by computing its sum of quality scores at
          // mismatches.
                    
          sum[i].resize( rpaths[i].size( ), 0 );
          for ( int j = 0; j < rpaths[i].isize( ); j++ )
          {    const vec<int>& p = rpaths[i][j];
               const pair<int,int>& a1 = aligns[i][ p[0] ];
               int start = -a1.second;
               basevector x( edges[a1.first], start,
                    edges[a1.first].isize( ) - start );
               for ( int z = 1; z < p.isize( ); z++ )
               {    basevector y( edges[ aligns[i][ p[z] ].first ], K-1, 
                         edges[ aligns[i][ p[z] ].first ].isize( ) - (K-1) );
                    x = Cat( x, y );    }
               for ( int z = 0; z < b.isize( ); z++ )
                    if ( b[z] != x[z] && quals[i][z] >= minq ) 
                         sum[i][j] += quals[i][z];    }
          SortSync( sum[i], rpaths[i] );
          int j;
          for ( j = 0; j < rpaths[i].isize( ); j++ )
               if ( sum[i][j] > sum[i][0] + Qthresh ) break;
          sum[i].resize(j), rpaths[i].resize(j);    }    }
