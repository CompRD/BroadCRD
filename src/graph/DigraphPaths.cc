/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "graph/DigraphPaths.h"
#include "graph/Digraph.h"
#include "system/TimeUtils.h"
#include "Set.h"
#include <numeric>

// Min Distance - basic version
template<class E> void MinDistance(const digraphE<E>& g,
				   const vec<vrtx_t>& vertex,
				   vec<vec<E> >& minDist) {

  // infinity represents no path between vertices
  const E infinity = GetMinDistanceInfinity(g);

  // Initialize minimum distance array and 
  int n = vertex.isize();
  minDist.resize(n);

  for (int i = 0; i < n; i++) {
    vec<E> dummy(n);
    for (int j = 0; j < n; j++) {
      if (i == j)
	dummy[j] = 0;
      else if (Member(g.From(vertex[i]), vertex[j]))
	dummy[j] = g.EdgeObjectByIndexFrom(vertex[i], Position(g.From(vertex[i]),vertex[j]));
      else
	dummy[j] = infinity;
    }
    minDist[i] = dummy;
  }
  
  // Build up answer
  for (int k = 0; k < n; k++) 
    for (int i = 0; i < n; i++) 
      for (int j = 0; j < n; j++)
	minDist[i][j] = Min(minDist[i][j], minDist[i][k] + minDist[k][j]);
}


// MinDistance - component friendly version
template<class E> void MinDistance(const digraphE<E>& g,
				   vec<int>& vertexIndex,
				   vec<int>& componentIndex,
				   vec<vec<vec<E> > >& minDist) {

  vertexIndex.resize(g.N());
  componentIndex.resize(g.N());

  // Obtain set of components
  vec<vec<int> > vertices;
  g.Components(vertices);
  minDist.resize(vertices.isize());

  // Find minimum distances for each component
  for ( int c = 0; c < vertices.isize( ); c++ ) {

    // Build vertex and component index

    int n = vertices[c].isize();

    for (int i = 0; i < n; i++) {
      vertexIndex[vertices[c][i]] = i;
      componentIndex[vertices[c][i]] = c;
    }
    
    // Calculate minimum distances for this component
    MinDistance(g, vertices[c], minDist[c]);
  }
}

template<class E> void SelfLoopMinDistance(const digraphE<E>& g,
					   vec<E>& selfLoopMinDist) {

  MinDistanceDatabase<int> minDistDb(g);
  const E infinity = minDistDb.GetInfinity();
  
  selfLoopMinDist.resize(g.N());
  for (vrtx_t v = 0; v < g.N(); v++) {
    E minDist = infinity;
    for (int k = 0; k < g.From(v).isize(); k++ ) {
      vrtx_t w = g.From(v)[k];
      if (v == w)
	minDist = Min(minDist, g.EdgeObjectByIndexFrom(v, k));
      else
	minDist = Min(minDist, g.EdgeObjectByIndexFrom(v, k)
		      + minDistDb.GetDistance(w, v) );
    }
    selfLoopMinDist[v] = minDist;
  }
}



// Sort criterion functor used by AllPathsBetweenSelectVertices to determine
// the optimal search order for the given vertex pairs.
class VertexPairSortCriterion : public binary_function<VertexPair, VertexPair, bool> {
public:
  VertexPairSortCriterion(const MinDistanceDatabase<int>& minDb)
    : minDistDb(minDb) {};
  // Order by 2nd vertex, then by minimum distance between vertices in pairs
  bool operator() (const VertexPair& lhs, const VertexPair& rhs) const {
    if (lhs.second < rhs.second)
      return true;
    if (lhs.second > rhs.second)
      return false;
    if (minDistDb.GetDistance(lhs.first, lhs.second) 
	< minDistDb.GetDistance(rhs.first, rhs.second)) {
      return true;
    }
    return false;
  }
private:
  const MinDistanceDatabase<int>& minDistDb;
};


template<class E> Bool AllPathsBetweenSelectVertices(const digraphE<E>& g,
					     const vec<VertexPair>& paired,
					     vec< vec< vec<vrtx_t> > > & allpaths,
					     const vec<int>& max_length,
					     const int maxpaths,
					     const Bool allow_self_loop,
					     const int timeout,
					     const Bool skip_on_timeout) {

  allpaths.resize(paired.isize());
  
  equiv_rel equiv;
  g.ComponentRelation(equiv);
  vec<vrtx_t> reps;
  equiv.OrbitRepsAlt(reps);

  MinDistanceDatabase<int> minDistDb(g);
  const int infinity = minDistDb.GetInfinity();

  set<int> unique_max;
  for (int k = 0; k < paired.isize(); k++) 
    unique_max.insert(max_length[k]);

  // Determine vertex pair seach order.
  vec<int> index;
  VertexPairSortCriterion cmp(minDistDb);
  SortIndex(paired, index, cmp);

  for (set<int>::iterator pos = unique_max.begin(); pos != unique_max.end(); ++pos) {
    int this_max_length = *pos;
    int skip_component = -1;
    //cout << Date() << ": Maximum possible insert size: " << this_max_length << "\n";
    for (int k = 0 ; k <index.isize(); k++) {
      int indexValue = index[k];
      if (max_length[indexValue] == this_max_length) {
        vrtx_t i = paired[indexValue].first;
 	vrtx_t j = paired[indexValue].second;
	if (!skip_on_timeout || (equiv.ClassId(i) != skip_component)) {	
	  vec<vec<vrtx_t> > newpaths;
	  if (!AllPathsAssisted(g, i,j, newpaths, allpaths, paired, minDistDb,
				this_max_length, maxpaths, allow_self_loop, timeout )) {
	    int component = BinPosition( reps, equiv.ClassId(i) );
	    cout << Date() << ": c = " << component
		 << ", i = " << i << ", j = " << j << " : Timeout - path of size "
		 << minDistDb.GetDistance(i,j) << " within max length\n";
	    if (skip_on_timeout) {
	      skip_component = equiv.ClassId(i);
	      cout << Date() << ": Skipping component " << component
		   << " for insert size " << this_max_length << "\n";
	    }
	  } else {
	    allpaths[indexValue] = newpaths;
	  }
	}
      }
    }
  }
  return True;
}


// Holds partial paths used by AllPathsAssisted
struct PartialPath {
  vec<vrtx_t> first;
  set<vrtx_t> second;
  int minDist;

  PartialPath(const vec<vrtx_t>& firstIn, const set<vrtx_t>& secondIn)
    : first(firstIn), second(secondIn) {
    minDist = 0;
  };

  PartialPath(const vec<vrtx_t>& firstIn, const set<vrtx_t>& secondIn, const int dist) 
    : first(firstIn), second(secondIn) {
    minDist = dist;
  };

  friend bool operator<(const PartialPath& lhs, const PartialPath& rhs) {
    return (lhs.first < rhs.first);
  }

};

template<class E> Bool AllPathsAssisted(const digraphE<E>& g,
					vrtx_t v, vrtx_t w,
					vec< vec<vrtx_t> >& paths,
					const vec<vec<vec< vrtx_t> > > & allpaths,
					const vec<VertexPair>& paired,
					const MinDistanceDatabase<int>& minDistDb, 
					const int max_length, const int maxpaths,
					const Bool allow_self_loop, int timeout ) {
  paths.clear( );
  set< PartialPath > partials;
  
  // Check that a path exists for v to w within max_length
  const int infinity = minDistDb.GetInfinity();
  int minLength = minDistDb.GetDistance(v,w);
  if ((minLength == infinity) || (max_length > 0 && minLength > max_length ))
    return True;
  
  // Check assist data
  Bool allpaths_assist = !allpaths.empty();
  if (allpaths_assist)
    AssertEq(allpaths.isize(), paired.isize() );
  // Set timeout
  TimeoutTimer::SetAlarm(timeout);
  
  // Starting value for paths is v
  vec<vrtx_t> just;
  set<vrtx_t> justs;
  just.push_back(v), justs.insert(v);

  // Special case - path from v to w when v==w
  if ( v == w ) {
    paths.push_back(just);
    if (!allow_self_loop)
      return True;
  }

  // Now build the answer.
  partials.insert( PartialPath( just, justs, 0 ) );
  while( !partials.empty( ) && !TimeoutTimer::TimeoutOccurred() ) {
    PartialPath p = *partials.begin( );
    partials.erase( partials.begin( ) );
    
    // Find distance so far (if needed)
    int distance = p.minDist;
    
    vrtx_t x = p.first.back( );
    for ( int i = 0; i < g.From(x).isize( ); i++ ) {
      vrtx_t y = g.From(x)[i];
      int edgeDist = g.EdgeObjectByIndexFrom(x, i);
      int minDistToEnd = minDistDb.GetDistance(y,w);
      if ((minDistToEnd == infinity) || 
	  (max_length > 0 && minDistToEnd > max_length - distance - edgeDist ))
	continue;
      
      if ( y == w ) {
	vec<vrtx_t> pf = p.first;
	pf.push_back(y);
	paths.push_back(pf);
	if ( maxpaths > 0 && (int) paths.size( ) > maxpaths ) 
	  return False;
	continue; 
      }
      
      if ( Member( p.second, y ) ) continue;
      
      int pos = BinPosition(paired, VertexPair(y,w));
      if (allpaths_assist && pos != -1 && !allpaths[pos].empty()) {
	for (int i = 0; i < allpaths[pos].isize(); ++i) {
	  bool looping = False;
	  int last = allpaths[pos][i].isize() - 1;
	  for (int j = 1; j < last; ++j) {
	    if (Member( p.second, allpaths[pos][i][j] ) ) {
	      looping = True; 
	      break;
	    }
	  }
	  if (looping)
	    continue;
	  
	  int distanceToEnd = distance + edgeDist;
	  for (int j = 0; j < last; j++) {
	    vrtx_t pfirst = allpaths[pos][i][j];
	    vrtx_t pnext = allpaths[pos][i][j+1];
	    for (int k = 0; k < g.From(pfirst).isize(); k++) {
	      if (g.From(pfirst)[k] == pnext) {
		distanceToEnd += g.EdgeObjectByIndexFrom(pfirst,k);
		break;
	      }
	    }
	    if (distanceToEnd > max_length)
	      break;
	  }
	  if (distanceToEnd > max_length)
	    continue;
	  
	  vec<vrtx_t> pf = p.first;
	  pf.append(allpaths[pos][i]);
	  paths.push_back(pf);
	  if ( maxpaths > 0 && (int) paths.size( ) > maxpaths ) 
	    return False;
	}
      } else {
	PartialPath q = p;
	q.first.push_back(y);
	q.second.insert(y);
	q.minDist = distance + edgeDist;
	partials.insert(q);
      }
      if ( maxpaths > 0 && (int) partials.size( ) > maxpaths ) 
                         return False;
    }
  }

  if (TimeoutTimer::TimeoutOccurred())
    return False;

  return True;
}

template<class E> int RemoveEqualParallelEdges(digraphE<E>& g) {
  vec<edge_t> edges_to_delete;
  for (int v = 0; v < g.N(); v++) {
    for (int j1 = 0; j1 < g.From(v).isize( ); j1++) {
      for (int j2 = j1 + 1; j2 < g.From(v).isize(); j2++) {
	if (g.From(v)[j2] != g.From(v)[j1])
	  continue;
	if (g.EdgeObjectByIndexFrom(v, j1) != g.EdgeObjectByIndexFrom(v, j2))
	  continue;
	edges_to_delete.push_back(g.EdgeObjectIndexByIndexFrom(v, j2));
      }
    }
  }
  UniqueSort(edges_to_delete);
  g.DeleteEdges(edges_to_delete);
  
  return edges_to_delete.isize();
}

template<class E> int RemoveLongerParallelEdges(digraphE<E>& g) {
  vec<edge_t> edges_to_delete;
  for (vrtx_t v = 0; v < g.N(); v++) {
    for (int j1 = 0; j1 < g.From(v).isize(); j1++) {
      for (int j2 = j1 + 1; j2 < g.From(v).isize(); j2++) {
	if (g.From(v)[j2] != g.From(v)[j1])
	  continue;
	if (g.EdgeObjectByIndexFrom(v, j1) <= g.EdgeObjectByIndexFrom(v, j2))
	  edges_to_delete.push_back(g.EdgeObjectIndexByIndexFrom(v, j2));
	else
	  edges_to_delete.push_back(g.EdgeObjectIndexByIndexFrom(v, j1));
      }
    }
  }
  UniqueSort(edges_to_delete);
  g.DeleteEdges(edges_to_delete);

  return edges_to_delete.isize();
}


// Instantiate for DigraphE<int>:

template void MinDistance(const digraphE<int>& g,
			  const vec<vrtx_t>& vertex,
			  vec<vec<int> >& minDist);

template void MinDistance(const digraphE<int>& g,
			  vec<int>& vertexIndex, 
			  vec<int>& componentIndex,
			  vec<vec<vec<int> > >& minDist);

template void SelfLoopMinDistance(const digraphE<int>& g,
				  vec<int>& selfLoopMinDist);

template Bool AllPathsBetweenSelectVertices(const digraphE<int>& g,
					    const vec<VertexPair>& paired,
					    vec< vec< vec<vrtx_t> > > & allpaths,
					    const vec<int>& max_length,
					    const int maxpaths,
					    const Bool allow_self_loop,
					    const int timeout,
					    const Bool skip_on_timeout);

template Bool AllPathsAssisted(const digraphE<int>& g,
			       vrtx_t v, vrtx_t w, vec< vec<vrtx_t> >& paths,
			       const vec<vec<vec< vrtx_t> > > & allpaths,
			       const vec<VertexPair>& paired,
			       const MinDistanceDatabase<int>& minDistDb, 
			       int max_length, int maxpaths,
			       const Bool allow_self_loop, int timeout ) ;

template int RemoveEqualParallelEdges(digraphE<int>& g);

template int RemoveLongerParallelEdges(digraphE<int>& g);
