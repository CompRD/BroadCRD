// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "CoreTools.h"
#include "math/Functions.h"
#include "graph/OddCycle.h"

// MinimalOddCycle.  Given a graph, represented by a symmetric square matrix of
// booleans, find a minimal odd cycle, returning the vertices which participate in 
// it (3 or 5 or 7 or ...), or no vertices, if there is no odd cycle.  Note that 
// there may be more than one odd cycle, but only one is returned, and there is 
// nothing canonical about the choice.  
//
// Note: the definition of "minimal" we use is that there is no odd cycle 
// properly contained within the one we are calling minimal.  There could be
// another smaller odd cycle, using different vertices.

void MinimalOddCycle( const vec< vec<Bool> >& edge, vec<int>& odd_cycle )
{    odd_cycle.clear( );
     int n = edge.size( );
     static vec<Bool> keep;
     keep.resize_and_set( n, False );
     static vec< vec<int> > neighbors;
     if ( neighbors.isize( ) < n ) neighbors.resize(n);
     for ( int i = 0; i < n; i++ )
          neighbors[i].clear( );
     for ( int i = 0; i < n; i++ )
     {    for ( int j = 0; j < n; j++ )
          {    if ( edge[i][j] ) neighbors[i].push_back(j);    }    }
     for ( int z = 0; z < n; z++ )
     {    
          // What all code in here does: if vertex z and all vertices set in
          // keep are removed, does there exist an odd cycle?  If so, set
          // keep[z] = True.

          int first_neg = 0;
          static vec<int> zeros, status;
          zeros.clear( );
          status.resize_and_set( n, -1 );

          // Meaning of status:
          //  1 means colored white
          //  2 means colored black
          //  0 means adjacent to a colored vertex
          // -1 means neither colored nor adjacent to colored.

          while(1)
          {    
               if ( zeros.nonempty( ) )
               {    int i = zeros.back( );
                    
                    // We know that vertex i is adjacent to a colored vertex.
                    // Find such a vertex, then color vertex i, then set status
                    // to 0 for vertices adjacent to i but uncolored.
                    // Run time: O(# of vertices adjacent to vertex i).

                    if ( !keep[i] && i != z )
                    {    for ( int x = 0; x < neighbors[i].isize( ); x++ )
                         {    int j = neighbors[i][x];
                              if ( !keep[j] && j != z
                                   && ( status[j] == 1 || status[j] == 2 ) )
                              {    if ( status[j] == 1 ) status[i] = 2;
                                   else status[i] = 1;
                                   zeros.resize( zeros.size( ) - 1 );
                                   for ( int y = 0; y < neighbors[i].isize( ); y++ )
                                   {    int k = neighbors[i][y];
                                        if ( !keep[k] && k != z && status[k] == -1 )
                                        {    status[k] = 0;    
                                             zeros.push_back(k);   }    }    
                                   break;    }    }    }
                    continue;    }
               int i;
               for ( i = first_neg; i < n; i++ )
                    if ( status[i] == -1 ) break;
               if ( i == n ) break;
               status[i] = 1;
               first_neg = i + 1;
               if ( i != z && !keep[i] )
               {    
                    // Set status to 0 for vertices adjacent to i but uncolored.

                    for ( int x = 0; x < neighbors[i].isize( ); x++ )
                    {    int k = neighbors[i][x];
                         if ( !keep[k] && k != z && status[k] == -1 ) 
                         {    status[k] = 0;    
                              zeros.push_back(k);    }    }    }    }
          for ( int i1 = 0; i1 < n; i1++ )
          {    if ( keep[z] ) break;
               if ( keep[i1] || i1 == z ) continue;
               for ( int x = 0; x < neighbors[i1].isize( ); x++ )
               {    int i2 = neighbors[i1][x];
                    if ( status[i1] == status[i2] && !keep[i2] && i2 != z )
                    {    keep[z] = True;    
                         break;    }    }    }    }
     if ( Sum(keep) == 0 ) return;
     for ( int i = 0; i < n; i++ )
          if ( !keep[i] ) odd_cycle.push_back(i);    }

// OddCycleVertices.  Given a graph, represented by a symmetric square matrix of
// booleans, find a non-canonically-chosen collection of vertices, which if removed 
// from the graph, eliminates all odd cycles.  It is obtained in the following 
// iterative fashion:
// 1. Start with an empty deletion list D.
// 2. Add all vertices appearing in a 3-cycle to D and delete them from the graph.
// 3. If there is no odd cycle in the graph, return D.
// 4. Pick a minimal odd cycle C "at random".  Add the vertices in C to D.  Pick one
//    of the vertices in C "at random", and delete it from the graph.
// 5. Go to #3.

void OddCycleVertices( const vec< vec<Bool> >& edge, vec<int>& to_delete )
{    static vec< vec<Bool> > e;
     e = edge;
     to_delete.clear( );
     int n = e.size( );
     for ( int i1 = 0; i1 < n; i1++ )
     {    for ( int i2 = i1+1; i2 < n; i2++ )
          {    if ( !edge[i1][i2] ) continue;
               for ( int i3 = i2+1; i3 < n; i3++ )
               {    if ( !edge[i2][i3] ) continue;
                    if ( !edge[i3][i1] ) continue;
                    to_delete.push_back( i1, i2, i3 );
                    for ( int j = 0; j < n; j++ )
                    {    e[i1][j] = e[j][i1] = False;
                         e[i2][j] = e[j][i2] = False;
                         e[i3][j] = e[j][i3] = False;    }    }    }    }
     while(1)
     {    static vec<int> odd_cycle;
          MinimalOddCycle( e, odd_cycle );
          if ( odd_cycle.empty( ) ) break;
          to_delete.append(odd_cycle);
          int k = odd_cycle[0];
          for ( int i = 0; i < n; i++ )
               e[i][k] = e[k][i] = False;    }
     UniqueSort(to_delete);    }
