/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   File: DigraphPaths

   Contains a set of functions and classes that are useful for path and
   distance finding within digraphE graphs.
   
   Instantiated for DigraphE<int>.
*/

#ifndef DIGRAPH_PATHS_H
#define DIGRAPH_PATHS_H

#include "graph/Digraph.h"

/**
   Function: GetMinDistanceInfinity

   Returns the very large value used to represents no path between vertices
   in the MinDistance functions. It is set to half the  maximum value of E
   for digraph<E>. The graph is passed only in order to discover its type.
*/
template<class E> E GetMinDistanceInfinity(const digraphE<E>& g) {
  // half maximum value of E so that infinity + infinity doesn't overflow
  return numeric_limits<E>::max() / 2;
}


/**
   Class: MinDistanceDatabase
   
   Database of the minimum distances between all vertex pairs in a DigraphE.

   Important - the graph must not have multiple edges between the same two
   vertices. It should be a simplified graph where all longer parallel edges
   have been removed leaving only the shortest edge.

   This information is stored in a fairly efficent way that will not explode in
   memory as the vertex count increases - so long as the graph is split into
   components. Interrogation is very fast (constant time). Building the database
   from a graph is also fast - O(CV^3) where V is the number of vertices in the
   component containing the most vertices, and C is the number of components.

   No path between a vertex pair is represented by the large value returned
   by the member function GetInfinity.
   
   Uses <MinDistance>
*/
template<class E> class MinDistanceDatabase {
public:
  /// Constructor: Builds database for DigraphE g
  MinDistanceDatabase(const digraphE<E>& g) {
    MinDistance(g, vertexIndex, componentIndex, minDistance);
    // infinity represents no path between vertices
    infinity = GetMinDistanceInfinity(g);
  }

  /// Constructor: Default
  MinDistanceDatabase() {}

  /// Returns the minimum distance through the graph from  vertex v to vertex w
  /// Returns large value <infinity> if no path exists between vertices
  inline E GetDistance(const vrtx_t v, const vrtx_t w) const {
    if (componentIndex[v] == componentIndex[w])
      return minDistance[componentIndex[v]][vertexIndex[v]][vertexIndex[w]];
    else
      return infinity;
  }

  /// Returns True if a path exists from vertex v to vertex w
  inline Bool IsConnected(const vrtx_t v,const vrtx_t w) const {
    return (GetDistance(v,w) != infinity);
  }

  /// Returns True if vertex v and vertex w are in the same component
  /// This does not imply that a path exists from v to w
  inline Bool IsSameComponent(const vrtx_t v,const vrtx_t w) const {
    return (componentIndex[v] == componentIndex[w]);
  }

  /// Returns the value used in the database to represent 'no path exists'
  inline E GetInfinity() const {
    return infinity;
  }

private:
  E infinity;
  vec<int> vertexIndex;
  vec<int> componentIndex;
  vec<vec<vec<E> > > minDistance;
};


/**
   Function: MinDistance - component friendly version

   Calculates the minimum distances between all vertex pairs in a DigraphE
   Modified version of Floyd-Warshall algorithm optimized for graphs that
   contain components. Fast, O(CV^3) where V is the number of vertices in the
   component containing the most vertices, and C is the number of components.

   Important - the graph must not have multiple edges between the same two
   vertices. It should be a simplified graph where all longer parallel edges
   have been removed leaving only the shortest edge.

   Returns the following:
   minDist        - the minimum distances between each vertex pair, indexed by
                    component, start vertex, end vertex.
   vertexIndex    - maps actual vertex id in the graph to the vertex index used
                    in minDist
   componentIndex - maps actual vertex id in the graph to the component index
                    used in minDist

   Example of how to determine the minimum path distance from v to w assuming
   that v and w are in the same component:
     distance = minDistance[componentIndex[v]][vertexIndex[v]][vertexIndex[w]];

   Exampe of how to determine if v and w are in the same component:
     componentIndex[v] == componentIndex[w]

   No path between a vertex pair is represented by the large value returned
   by the function GetMinDistanceInfinity.

   Uses <MinDistance> - basic version
*/
template<class E> void MinDistance(const digraphE<E>& g,
				   vec<int>& vertexIndex,
				   vec<int>& componentIndex,
				   vec<vec<vec<E> > >& minDist);

/**
   Function: MinDistance - basic version

   Calculates the minimum distances between the list of vertices passed into
   the function for the DigraphE. The list of vertices must contain either all
   vertices in the graph, or all the vertices of a single component. For a
   graph with many components it is much faster and more memory efficient to
   use other version of MinDistance or the MinDistanceDatabase class.

   Direct implimentation of the Floyd-Warshall algorithm.
   Has complexity O(V^3) where V is the number of vertices.

   Important - the graph must not have multiple edges between the same two
   vertices. It should be a simplified graph where all longer parallel edges
   have been removed leaving only the shortest edge.

   Returns the following:
   minDist - the minimum distances between each vertex pair, indexed by
             start vertex and end vertex.

   Example of how to determine the minimum path distance from v to w:
     distance = minDistance[v][w];

   No path between a vertex pair is represented by the large value returned
   by the function GetMinDistanceInfinity.

*/
template<class E> void MinDistance(const digraphE<E>& g,
				   const vec<vrtx_t>& vertex,
				   vec<vec<E> >& minDist);


/**
   Function: SelfLoopMinDistance

   Calculates the minimum distance, if any, for all self loops in the graph.
   A self loop is a path from a vertex back to itself.

   Important - the graph must not have multiple edges between the same two
   vertices. It should be a simplified graph where all longer parallel edges
   have been removed leaving only the shortest edge.

   Returns the following:
   selfLoopMinDist - the minimum self loop distance for each vertex in the graph.

   No self path for a vertex is represented by the large value returned
   by the function GetMinDistanceInfinity.
*/
template<class E> void SelfLoopMinDistance(const digraphE<E>& g,
					   vec<E>& selfLoopMinDist);


/**
   Function: AllPathsBetweenSelectPairs

   Find all possible paths between each pair of vertices as specified in the
   argument <paired> - a vec of VertexPair. This is much faster than solving for
   each pair of vertices in turn using the AllPaths function.

   Important - the graph must not have multiple edges between the same two
   vertices. It should be a simplified graph where all longer parallel edges
   have been removed leaving only the shortest edge.

   The set of vertex pairs to use is passed as a vec of VertexPair, where the first
   vertex of the pair is the starting point, and the 2nd vertex is the end point.
   The maximium allowed path length for each vertex pair is given in max_length, where
   the order matches the vertex pair order in paired. Set all values to zero for an
   exhaustive search with no upper limit on the path lengths returned.

   The resulting vec of paths found is returned in allpaths, again the order matches
   the vertex pair order in paired.

   Returns true if no problems were encoutered.

   If all paths cannot be found within 100 seconds (default) then a timeout will occur.
   If skip_on_timeout is set to true then this function will skip the remaining
   pairs of vertices of that max_length size on the particular component and start
   again on the next component.

   If maxpaths >= 0, fail, returning False if greater than that many paths or partial
   paths are found during the search.
   If allow_self_loop and v = w, then loops from v to v are allowed (but v cannot
   appear more than twice).
*/


template<class E> Bool AllPathsBetweenSelectVertices(const digraphE<E>& g,
					     const vec<VertexPair>& paired,
					     vec< vec< vec<vrtx_t> > > & allpaths,
					     const vec<int>& max_length,
					     const int maxpaths = -1,
					     const Bool allow_self_loop = False,
					     const int timeout = 100,
					     const Bool skip_on_timeout = False);


/**
   Function: AllPathsAssisted

   Find all possible paths from v to w, using any already know paths (allpaths) and
   minimum distance information (minDistDb) to accelerate the search.

   Important - the graph must not have multiple edges between the same two
   vertices. It should be a simplified graph where all longer parallel edges
   have been removed leaving only the shortest edge.

   See the implimentation of AllPathsBetweenSelectPairs for an example of how to use
   this function.

   If max_length > 0 then enforce this upper limit on the path length.
   If maxpaths >= 0, fail, returning False if greater than that many paths or partial
   paths are found during the search.
   If allow_self_loop and v = w, then loops from v to v are allowed (but v cannot
   appear more than twice).
*/
template<class E> Bool AllPathsAssisted(const digraphE<E>& g,
				vrtx_t v, vrtx_t w, vec< vec<vrtx_t> >& paths,
				const vec<vec<vec< vrtx_t> > > & allpaths,
				const vec<VertexPair>& paired,
				const MinDistanceDatabase<int>& minDistDb, 
				const int max_length = -1, const int maxpaths = -1,
				const Bool allow_self_loop = False, int timeout = 0 );

/**
   Function: RemoveEqualParallelEdges

   Removes parallel edges that are equal in length, leaving only one between
   each vertex pair. Parallel edges of differing length are unaffected.
   The edges removed are a subset of those removed by RemoveLongerParallelEdges.
*/
template<class E> int RemoveEqualParallelEdges(digraphE<E>& g);

/**
   Function: RemoveLongerParallelEdges

   Removes parallel edges, leaving only the shortest between each vertex pair.
   The edges removed are a superset of those removed by RemoveEqualParallelEdges.
*/
template<class E> int RemoveLongerParallelEdges(digraphE<E>& g);

#endif
