/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "FeudalMimic.h"
#include "Map.h"
#include "PrintAlignment.h"
#include "PairsManager.h"
#include "graph/Digraph.h"
#include "graph/DigraphPaths.h"
#include "math/Functions.h"
#include "paths/AlignPairsToHyper.h"
#include "paths/AlignSeqsToHyper.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/SeqOnHyper.h"

// Computes alignments of reads to the hyperkmer path directly using Perfectlookup

void AlignPairsToHyper( const HyperKmerPath& h, const HyperBasevector& hb,
			// note that hb must be HyperBasevector(h)
			const vecbasevector& reads, const PairsManager& pairs,
			const String& sub_dir,
			vec<CompressedSeqOnHyper>& csaligns,
			vec<vec<int> >& csaligns_index,
			vec< IndexPair >& pair_places, 
			Bool verbose,  const int max_insert_size,
                        vecbasevector* p_cached_edges,
                        vec<partial_perfect_align>* p_cached_ppaligns)
{    
  AlignSeqsToHyper(h, hb, reads, sub_dir, csaligns, verbose, p_cached_edges, 
       p_cached_ppaligns );

  BuildAlignmentIndex( csaligns, reads.size(), csaligns_index, verbose );

  FindPairPlacements( h, hb, csaligns, csaligns_index, reads.size(), pairs,
		      pair_places, verbose, sub_dir, max_insert_size );
}

// EvaluatePair: given alignments ap1, ap2 of the reads of a read pair P to a
// HyperKmerPath, determine if the insert could be placed contiguously on the
// HyperKmerPath, without stretching too much.  Return True in that case.

Bool EvaluatePair( 
     const SeqOnHyper& ap1,            // alignment of first read
     const SeqOnHyper& ap2,            // alignment of second read
     const int& sep, const int& stdev, // the read pair
     double dev_mult,                  // how much stretch is allowed in devs
     const HyperBasevector& hb,        // from the HyperKmerPath
     const vec<int>& to_right_vertex,  // index from edge to right vertex
     const vec<int>& to_left_vertex,   // index from edge to left vertex
     const digraphE<int>& G,           // mirrors hb but edges are lengths in kmers
     const vec<VertexPair>& paired,
     const vec<Bool>& have_path,
     const vec<int>& min_path,
     const vec<int>& max_path,
     const vec<int>& a,
     const vec<int>& b )
{
  int K = hb.K( );
  if ( ap1.Rc1( ) || !ap2.Rc1( ) )
    return False;
  int N1 = ap1.N( ), N2 = ap2.N( );
  int e1 = ap1.Id2( N1 - 1 ), e2 = ap2.Id2(0);
  int sep1 = hb.EdgeLength(e1) - ap1.Pos2( N1 - 1 ), sep2 = ap2.pos2(0);
  int calc_sep = sep1 + sep2 - K + 1;
  int w1 = to_right_vertex[e1], v1 = to_left_vertex[e2];
  static vec<int> D;
  int min_sep 
       = int( round( double(sep) - dev_mult * double(stdev) - double(calc_sep) ) );
  int max_sep 
       = int( round( double(sep) + dev_mult * double(stdev) - double(calc_sep) ) );
  if ( e1 == e2 && ap1.Pos2( N1 - 1 ) <= ap2.pos2(0) ) {
    int sepdiff = ap2.pos2(0) - ap1.Pos2( N1 - 1 ) - sep;
    if ( double( Abs(sepdiff) ) <= double(stdev) * dev_mult )
      return True;
  }

  int pos = BinPosition(paired, VertexPair(w1,v1));
  if ( pos < 0 || pos >= have_path.isize()  ||  !have_path[pos] ) return False;
  if ( min_path[pos] > max_sep ) return False;
  if ( max_path[pos] < min_sep ) return False;
  if ( min_sep <= min_path[pos] && min_path[pos] <= max_sep ) return True;
  if ( min_sep <= max_path[pos] && max_path[pos] <= max_sep ) return True;
  
  int A = a[pos], B = b[pos];
  if ( B > 0 ) {
    if ( A <= min_sep && B <= max_sep - min_sep ) return True;
    if ( min_sep <= A && A <= max_sep ) return True;
  }
  G.Distance( w1, v1, max_sep, D );
  for ( int u = 0; u < D.isize( ); u++ )
    D[u] += calc_sep - sep;
  for ( int u = 0; u < D.isize( ); u++ ) {
    double x = D[u] / double(stdev);
    if ( x >= -dev_mult && x <= dev_mult ) return True;
  }
  
  return False;
}

// EvaluatePair: given alignments cs1, cs2 of the reads of a read pair P to a
// HyperKmerPath, determine if the insert could be placed contiguously on the
// on a single edge of the HyperKmerPath, without stretching too much.
// Return True in that case.
// Special case of the more general EvalutePair function. Optimized for single
// edge condition using CompressedSeqOnHyper information only.

Bool EvaluatePairOnSingleEdge( 
     const CompressedSeqOnHyper& cs1,     // alignment of first read
     const CompressedSeqOnHyper& cs2,     // alignment of second read
     const int& sep, const int& stdev,    // the read pair
     double dev_mult )                    // how much stretch is allowed in devs
{
  if ( cs1.Rc1( ) || !cs2.Rc1( ) )
    return False;
  if (!cs1.SingleEdge() || !cs2.SingleEdge())
    return False;

  const PartialSeqOnHyper& p_cs1 = cs1.Part0();
  const PartialSeqOnHyper& p_cs2 = cs2.Part0();
  int e1 = p_cs1.Id2();
  int e2 = p_cs2.Id2();
  if ( e1 == e2 && p_cs1.Pos2() <= p_cs2.pos2() ) { 
    int sepdiff = p_cs2.pos2() - p_cs1.Pos2() - sep;
    if ( double( Abs(sepdiff) ) <= double(sep) * dev_mult )
      return True;
  }
  return False;  
}

// Find placements of read pairs that are not too stretched. 

void FindPairPlacements( const HyperKmerPath& h, const HyperBasevector& hb, 
     const vec<CompressedSeqOnHyper>& csaligns, const vec<vec<int> >& csaligns_index,
     int nreads, const PairsManager& pairs, vec< IndexPair >& pair_places, 
     Bool verbose, const String& sub_dir, const int max_insert_size,
     const int timeout, Bool skip_on_timeout )
{
  if ( verbose ) {
    cout << Date( ) << ": Finding pair placements using HyperKmerPath alignments" << endl;
    cout << Date( ) << ": graph has " << h.ConnectedComponents()
	 << " components, " << hb.EdgeObjectCount() << " edges, and " << h.N()
	 << " vertices" << endl;
    PrintMemUsage( );
  }

  vec<int> to_left_vertex, to_right_vertex;
  h.ToLeft(to_left_vertex);
  h.ToRight(to_right_vertex);
  vec<int> edgelengths(h.EdgeObjectCount());
  for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
    edgelengths[i] = h.EdgeObject(i).KmerCount( ) ;
  digraphE<int> G( h, edgelengths );
  int N = h.N( );

  Destroy(edgelengths);

  // Create a graph H that is the same as h, except that we remove 
  // parallel edges of equal length.
  digraphE<int> H = G;
  RemoveEqualParallelEdges(H);

  // Create a graph L that is the same as H, except that we remove 
  // parallel edges until only the shortest one remains.
  digraphE<int> L = H;
  RemoveLongerParallelEdges(L);

  equiv_rel equiv;
  L.ComponentRelation(equiv);

  pair_places.clear( );

  // We precompute a bunch of stuff to help EvaluatePair go faster (two parts).

  // For each pair of vertices v, w, determine if there is a directed path 
  // from v to w.  Find the minimum length of such a path.  If there are a >= 0 
  // and b > 0 such that all the elements of {a+bt: t >= 0} occur as path 
  // lengths, find a and b, with b as small as possible and subject to that, a as 
  // small as possible.  Find the maximum length of a path if there is a max.  

  if ( verbose ) {
    if (max_insert_size == -1)
      cout << Date( ) << ": finding vertex pairs for all inserts" << endl;
    else
      cout << Date( ) << ": finding vertex pairs for inserts up to " 
	   << max_insert_size << " bases" << endl;
    PrintMemUsage( );
  }

  typedef map<VertexPair, int> VertexPairIntMap;
  VertexPairIntMap maxPairDist;
  int maxDist = 0;
  for ( size_t i = 0; i < pairs.nPairs( ); i++ ) {
    const int sep = pairs.sep(i);
    const int stdev = pairs.sd(i);
    if (max_insert_size != -1 && sep > max_insert_size)
      continue;
    int pair_sep = sep + stdev * 6;
    maxDist = max( maxDist, pair_sep );
    int id1x = pairs.ID1(i), id2x = pairs.ID2(i);
    for ( int i1 = 0; i1 < csaligns_index[id1x].isize( ); i1++ ) {
      for ( int i2 = 0; i2 < csaligns_index[id2x].isize( ); i2++ ) {
	int j1 = csaligns_index[id1x][i1];
	int j2 = csaligns_index[id2x][i2];
	if ( csaligns[j1].Rc1( ) == csaligns[j2].Rc1( ) )
	  continue;
	if ( csaligns[j1].Rc1( ) )
	  swap( j1, j2 );
	if ( csaligns[j1].Rc1( ) || !csaligns[j2].Rc1( ) )
	  continue;

	// Find vertex to the left of the 2nd read
	int e2 = csaligns[j2].Part0().Id2();
	int v1 = to_left_vertex[e2];

	// Assume 1st read is contained on a single edge (might not be true)
	// Find vertex to the right of the 1st read
	int e1 = csaligns[j1].Part0().Id2();
	int w1 = to_right_vertex[e1];

	// Determine if both reads are aligned to the same component
 	if (!equiv.Equiv(w1,v1))
	  continue;

	// Now check that 1st read is contained on a single edge, update values if not
	if (!csaligns[j1].SingleEdge()) {
	  SeqOnHyper ap1;
	  csaligns[j1].DecompressInto(ap1, hb);
	  e1 = ap1.Id2( ap1.N( ) - 1 );
	  w1 = to_right_vertex[e1];
	}

	VertexPair this_pair(w1, v1);
	int current_value = maxPairDist[this_pair];
	if (current_value < pair_sep)
	  maxPairDist[this_pair] = pair_sep;

      }
    }
  }
  
  vec<VertexPair> potentialPairs(keys(maxPairDist));
  vec<int> maxPairSep(values(maxPairDist));
  maxPairDist.clear();
  int potPairSize = potentialPairs.isize();
  
  // Find minimum loop from a vertex to itself.

  if ( verbose ) {
    cout << Date( ) << ": find minimum loops" << endl;
    PrintMemUsage( );
  }

  vec<int> min_loop( N, 0 );
  int selfLoopInfinity = GetMinDistanceInfinity(L);
  SelfLoopMinDistance(L, min_loop);
  
  if ( verbose ) {
    cout << Date( ) << ": found " << N - min_loop.CountValue(selfLoopInfinity)
	 << " loops" << endl;
    cout << Date( ) << ": finding paths between " << potPairSize << " vertex pairs" << endl;
    PrintMemUsage( );
  }
  
  vec< vec< vec<int> > > allpaths(potPairSize);

  AllPathsBetweenSelectVertices(L,potentialPairs, allpaths, maxPairSep, -1 , True,
				timeout, skip_on_timeout);

  longlong checksum1 = 0;
  longlong checksum2 = 0;
  for (int i = 0; i < allpaths.isize(); ++i) {
    checksum1 += allpaths[i].isize();
    for (int j = 0; j < allpaths[i].isize(); ++j)
      checksum2 += allpaths[i][j].isize();
  }
  
  Destroy(maxPairSep);

  // Finally find the distances.
  if ( verbose ) {
    cout << Date( ) << ": find distances between vertex pairs" << endl;
    PrintMemUsage( );
  }

  const int infinity = 1000000000;
  vec<Bool> have_path(potPairSize);
  vec<int> min_path(potPairSize, 0), max_path(potPairSize, infinity);
  vec<int> a(potPairSize, 0), b(potPairSize, 0);
  
  for (int k = 0; k < potPairSize; k++) {
    int v = potentialPairs[k].first;
    int w = potentialPairs[k].second;
      vec< vec<int> >& paths = allpaths[k];
      have_path[k] = paths.nonempty( );    

      int m = infinity, M = 0;
      for ( int i = 0; i < paths.isize( ); i++ ) {
	// Set d to the minimum length of paths[i], and D to the 
	// maximum.  The length is variable because there can be more
	// than one edge between two given vertices.
	
	int d = 0, D = 0;
	for ( int j = 1; j < paths[i].isize( ); j++ ) {
	  int x = paths[i][j-1], y = paths[i][j];
	  d += H.MinEdge( x, y );
	  D += H.MaxEdge( x, y );    
	  if ( d > maxDist )
	    break;
	}
	if ( d > maxDist )
	  continue;
	m = Min( m, d );
	M = Max( M, d );
	for ( int j = 0; j < paths[i].isize( ); j++ ) {
	  int ml = min_loop[ paths[i][j] ];
	  if ( ml != selfLoopInfinity ) {
	    if ( b[k] == 0 || ml < b[k] || ( ml == b[k] && d < a[k] ) ) {
	      a[k] = d;
	      b[k] = ml;
	    }
	  }
	}
      }

      if ( m == infinity ) { 
	have_path[k] = false;
	continue;
      }
      if ( b[k] == 0 ) max_path[k] = M;
      min_path[k] = m;
 
  }

  // Go through the pairs.
  
  if ( verbose ) {
    cout << Date( ) << ": find valid pair placements" << endl;
    PrintMemUsage( );
  }
  
  for ( size_t i = 0; i < pairs.nPairs( ); i++ ) {
    const int sep = pairs.sep(i);
    const int stdev = pairs.sd(i);
    int id1x = pairs.ID1(i), id2x = pairs.ID2(i);
    if (verbose)
      cout << "\nalignments of pair " << i << "\n";
    int count = 0;
    int pair_places_size = pair_places.size( );
    Bool trusted = True;
    for ( int i1 = 0; i1 < csaligns_index[id1x].isize( ); i1++ ) {
      for ( int i2 = 0; i2 < csaligns_index[id2x].isize( ); i2++ ) {
	int j1 = csaligns_index[id1x][i1];
	int j2 = csaligns_index[id2x][i2];
	if ( csaligns[j1].Rc1( ) == csaligns[j2].Rc1( ) )
	  continue;
	if ( csaligns[j1].Rc1( ) )
	  swap( j1, j2 );

	// Find vertex to the left of the 2nd read
	int v1 = to_left_vertex[csaligns[j2].Part0().Id2()];

	// Assume 1st read is contained on a single edge (might not be true)
	// Find vertex to the right of the 1st read
	int w1 = to_right_vertex[csaligns[j1].Part0().Id2()];

	// Determine if both reads are aligned to the same component
 	if (!equiv.Equiv(w1,v1))
	  continue;

	SeqOnHyper ap1, ap2;
	// Does pair land on a single edge of hyperkmerpath (quick)
	Bool have3_0 = EvaluatePairOnSingleEdge( csaligns[j1], csaligns[j2], sep, stdev, 3.0);
	if (!have3_0) {
	  // Pair not on single edge so perfom more complete evaluation (slow)
	  csaligns[j1].DecompressInto(ap1, hb);
	  csaligns[j2].DecompressInto(ap2, hb);
	  have3_0 = EvaluatePair( ap1, ap2, sep, stdev, 3.0, hb, to_right_vertex, to_left_vertex,
				  G, potentialPairs, have_path, min_path, max_path, a, b );
	}
	if ( !have3_0 ) {
	  Bool have3_5 = EvaluatePair( ap1, ap2, sep, stdev, 3.5, hb, to_right_vertex,
				       to_left_vertex, G, potentialPairs, 
				       have_path, min_path, max_path, a, b ); 
	  if (have3_5) trusted = False;
	}
	if (have3_0) {
	  if (verbose) {
	    if ( count++ > 0 )
	      cout << "----------------------------------------------------\n";
	    ap1.Print( cout, to_left_vertex, to_right_vertex );
	    ap2.Print( cout, to_left_vertex, to_right_vertex );
	    cout << sep << " +/- " << stdev << "\n";
	    if ( !trusted )
	      cout << "[NOT TRUSTED.]\n";
	  }
	  pair_places.push_back( make_pair( j1, j2 ) );
	}
      }
    }
    if ( !trusted )
      pair_places.resize(pair_places_size);
  }
  
  if ( verbose ) {
    cout << Date( ) << ": Found " << pair_places.size()  << " pair placements" << endl;
    PrintMemUsage( );
  }
}


// Find the bases on the HyperBasevector edges which are covered by alignments of 
// read pairs that are not too stretched.

void FindWellCovered( int min_sep, int max_sep, Bool use_fw, Bool use_rc,
     const HyperKmerPath& h, const HyperBasevector& hb, 
     const vec<CompressedSeqOnHyper>& csaligns, const vec< IndexPair >& pair_places,
     int nreads, const PairsManager& pairs,
     vecbitvector& cov, Bool verbose )
{
  vecbasevector edges;
  for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
    edges.push_back_reserve( hb.EdgeObject(i) );
  Mimic( edges, cov );

  for ( int i = 0; i < pair_places.isize( ); i++ ) {
    SeqOnHyper ap1, ap2;
    csaligns[ pair_places[i].first ].DecompressInto(ap1, hb);
    csaligns[ pair_places[i].second ].DecompressInto(ap2, hb);
    const int sep = pairs.sep( pairs.getPairID( ap1.Id1( ) ) );
    if ( sep < min_sep || sep > max_sep )
      continue;
    int N1 = ap1.N( ), N2 = ap2.N( );
    if ( ( !ap1.Rc1( ) && use_fw ) || ( ap1.Rc1( ) && use_rc ) ) {
      for ( int u = 0; u < N1; u++ ) {
	int id2 = ap1.Id2(u);
	for ( int v = ap1.pos2(u); v < ap1.Pos2(u); v++ )
	  cov[id2].Set( v, True );
      }
    }
    if ( ( !ap2.Rc1( ) && use_fw ) || ( ap2.Rc1( ) && use_rc ) ) {
      for ( int u = 0; u < N2; u++ ) {
	int id2 = ap2.Id2(u);
	for ( int v = ap2.pos2(u); v < ap2.Pos2(u); v++ )
	  cov[id2].Set( v, True );
      }
    }
  }
}


// Deletes edges of HyperKmerPath that are not covered by bases in a good pair.

void FindPoorlyCovered( HyperKmerPath& h, HyperBasevector& hb, 
     const KmerBaseBroker& kbb, const vec<CompressedSeqOnHyper>& csaligns, 
     const vec< IndexPair >& pair_places, const vecbasevector& reads, 
     const PairsManager& pairs, Bool remove, int verbosity,
     const Bool NEW_POORLY_COVERED )
{
  int nreads = reads.size( );

     vec<int> to_right_vertex, to_left_vertex;
     h.ToRight(to_right_vertex), h.ToLeft(to_left_vertex);
     vecbitvector cov;
     // TODO: potentially dangerous truncation of index by to_delete
     vec<int> to_delete;
     vec<int> trim_left( h.EdgeObjectCount( ), 0 );
     vec<int> trim_right( h.EdgeObjectCount( ), 0 );

     // Find the parts of the HyperKmerPath that are insufficiently covered by the 
     // full set of reads.

     if (NEW_POORLY_COVERED) 
     {    for ( int pass = 1; pass <= 3; pass++ ) 
          {    int low, high;
               if ( pass == 1 ) low = 0, high = 1000000;
               else if ( pass == 2 ) low = 0, high = 1000;
               else low = 1000, high = 1000000;
               FindWellCovered( low, high, True, True, h, hb, csaligns, pair_places,
	            nreads, pairs, cov, verbosity >= 2 );
               if ( verbosity >= 1 ) 
	            cout << "\ntesting for uncovered, pass " << pass << "\n"; 
               for ( size_t u = 0; u < cov.size( ); u++ )
               {    for ( unsigned int v = 0; v < cov[u].size( ); v++ )
                    {    if ( !cov[u][v] ) 
                         {    unsigned int w;
	                      for ( w = v + 1; w < cov[u].size( ); w++ )
	                           if ( cov[u][w] ) break;
	                      if ( verbosity >= 1 )
                              {    cout << "uncovered: " << BaseAlpha(u) << "." 
                                        << v << "-" << w 
                                        << " (l = " << hb.EdgeLength(u) << ")";    }
                              int L = to_left_vertex[u], R = to_right_vertex[u];
                              int N = hb.EdgeLength(u);
	                      if ( ( h.Source(L) && w <= 20 )
		                   || ( h.Sink(R) && N - v <= 20 ) ) 
                              {    if (verbosity >= 1) cout << " (trimming)\n";
	                           if ( h.Source(L) && w <= 20 )
                                        trim_left[u] = Max( trim_left[u], static_cast<int>(w) );
		                   if ( h.Sink(R) && N - v <= 20 )
                                   {    trim_right[u] 
                                             = Max( trim_right[u], static_cast<int>(N-v) );    }    }
	                      else if ( N >= 10000 ) 
                              {    if ( verbosity >= 1 )
		                        cout << " (ignoring - large edge)\n";    }
	                      else if ( pass == 2 && w - v < 10 ) 
                              {    if (verbosity >= 1) cout << " (ignoring)\n";    }
	                      else 
                              {    if ( verbosity >= 1 ) cout << "\n";
	                           to_delete.push_back(u);    }
	                      v = w - 1;    }    }    }
               if ( verbosity >= 1 ) cout << "\n";    }

  } else {
    FindWellCovered( 0, 1000000000, True, True, h, hb, csaligns, pair_places, 
		     nreads, pairs, cov, verbosity >= 2 );
    if ( verbosity >= 1 )
      cout << "\n"; 
    for ( size_t u = 0; u < cov.size( ); u++ ) {
      for ( unsigned int v = 0; v < cov[u].size( ); v++ ) {
	if ( !cov[u][v] ) {
	  unsigned int w;
	  for ( w = v + 1; w < cov[u].size( ); w++ )
	    if ( cov[u][w] )
	      break;
	  if ( verbosity >= 1 )
	    cout << "uncovered: " << BaseAlpha(u) << "." << v << "-" << w 
		 << " (l = " << hb.EdgeLength(u) << ")";
	  if ( ( h.Source( to_left_vertex[u] ) && w <= 5 )
	       || ( h.Sink( to_right_vertex[u] ) && hb.EdgeLength(u) - v <= 5 ) ) {
	    if ( verbosity >= 1 )
	      cout << " (ignoring)\n";
	  } else if ( hb.EdgeLength(u) >= 10000 ) {
	    if ( verbosity >= 1 ) 
	      cout << " (ignoring - large edge)\n";
	  } else {
	    if ( verbosity >= 1 )
	      cout << "\n";
	    to_delete.push_back(u);
	  }
	  v = w - 1;
	}
      }
    }
    if ( verbosity >= 1 )
      cout << "\n";
  }

     // Clean the HyperKmerPath.

     if (remove) 
     {    for ( int u = 0; u < h.EdgeObjectCount( ); u++ )
          {    if ( trim_left[u] > 0 || trim_right[u] > 0 )
               {    if ( trim_left[u] + trim_right[u] >= h.EdgeLength(u) )
                    {    to_delete.push_back(u);
                         continue;    }
                    else
                    {    KmerPath& p = h.EdgeObjectMutable(u);
                         KmerPathLoc start = p.Begin( ), stop = p.End( );
                         start += trim_left[u], stop -= trim_right[u];
                         KmerPath pnew;
                         p.CopySubpath( start, stop, pnew );
                         p = pnew;    }    }    }
          UniqueSort(to_delete);
          h.DeleteEdges(to_delete);
          h.RemoveUnneededVertices( );
          h.RemoveDeadEdgeObjects( );
          hb = HyperBasevector( h, kbb );    }    }
