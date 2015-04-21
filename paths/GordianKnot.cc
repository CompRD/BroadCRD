///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// For class/module documentation, see GordianKnot.h
#include "PrettyPrintTable.h"
#include "graph/Digraph.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/AlignPairsToHyper.h"

// This takes care of all #includes for GordianKnot member variables' data types
#include "paths/GordianKnot.h" 






/* RemoveUnsupportedAdjacencies
 *
 *
 * BACKGROUND
 *
 * Every vertex in the HyperKmerPath contains a set of incoming edges and a set
 * of outgoing edges.  Let's call these sets A and B respectively, and the
 * number of incoming/outgoing edges nA and nB.
 *
 * By definition, a vertex represents a (K-1)-mer overlap between the edges in A
 * and the edges in B.  (From here we ignore vertices with nA = 0 or nB = 0,
 * because they do not imply any overlaps.)  We refer to each of these (K-1)-mer
 * overlaps as an "adjacency".
 *
 * Each vertex implies an adjacency between every edge in A and every edge in B,
 * for a total of nA*nB adjacencies.  Denote these as A0-B0, A0-B1, etc.  Some
 * of these adjacencies are genomic, but some of them have been spuriously
 * introduced at some stage during the graph-building process.
 *
 * Here we attempt to determine which adjacencies Ai-Bj are true by looking for
 * supporting evidence in the fragment reads and their pairs.  If there are at
 * least N_READS_TO_SUPPORT reads *or* at least N_PAIRS_TO_SUPPORT pairs that
 * fall on both Ai and Bj, the Ai-Bj adjacency is considered supported.
 *
 *
 * ALGORITHM STEPS
 *
 * 1. Place the short-insert pairs on the HyperBasevector via AlignPairsToHyper.
 * 2. For every vertex v in the HyperKmerPath, find all the short-insert pairs
 *    that bridge the vertex.  Determine which adjacencies Ai-Bj are supported.
 * 3. Split the vertex v and re-route its edges A-B in order to remove
 *    unsupported adjacencies while preserving all supported adjacencies.
 *
 *
 * FUNCTION INPUTS
 *
 * -- N_READS_TO_SUPPORT: Used in algorithm step 2, above.  An adjacency is
 *    considered "supported" if there are at least this many reads that fall on
 *    both edges.
 *
 * -- N_PAIRS_TO_SUPPORT: Used in algorithm step 2, above.  An adjacency is
 *    considered "supported" if there are at least this many pairs of reads
 *    (with the reads placed properly within the pair) that fall on both edges.
 *
 * You can set N_READS_TO_SUPPORT or N_PAIRS_TO_SUPPORT to -1; if you do this,
 * no number of reads/pairs will support an adjacency.  If you set both of these
 * to -1, the function will abort rather than break up the entire graph into
 * disconnected edges.
 *
 ******************************************************************************/
void
GordianKnot::RemoveUnsupportedAdjacencies( const int N_READS_TO_SUPPORT, const int N_PAIRS_TO_SUPPORT )
{
  cout << Date() << ": BEGIN GordianKnot::RemoveUnsupportedAdjacencies!" << endl;
  if ( N_READS_TO_SUPPORT == -1 && N_PAIRS_TO_SUPPORT == -1 ) {
    cout << "WARNING: You cannot run with N_READS_TO_SUPPORT = N_PAIRS_TO_SUPPORT = -1.\nIt would break up the entire graph into disconnected edges.\nI am aborting the function instead.";
    return;
  }
  
  
  // This fills _edge_coverage and _C_hbv.
  getCoverageByFragReads( );
  
  
  cout << Date() << ": Examining each vertex and its adjacencies" << endl;
  int total_adjs = 0, total_supported_adjs = 0;
  int total_vs_split = 0, total_vs_added = 0;
  
  if ( _verbose ) {
    String info;
    if ( N_READS_TO_SUPPORT == -1 ) info = ToString( N_PAIRS_TO_SUPPORT ) + " pairs";
    else if (N_PAIRS_TO_SUPPORT == -1 ) info = ToString( N_READS_TO_SUPPORT ) + " read(s)";
    else info = ToString( N_READS_TO_SUPPORT ) + " read(s) *or* " + ToString( N_PAIRS_TO_SUPPORT ) + " pair(s)";
    cout << "*** An adjacency is considered supported if it is verified by at least " << info << endl;
  }
  
  // Loop over every vertex v (that existed at the beginning of this function.)
  for ( int v = 0; v < _n_vertices; v++ ) {
    
    // Find the sets A and B and the numbers nA and nB.
    // (See this function's documentation for an explanation of variable names.)
    vec<int> A = _hkp.ToEdgeObj(v), B = _hkp.FromEdgeObj(v);
    int nA = A.size(), nB = B.size();
    if ( nA == 0 || nB == 0 ) continue;
    
    // Create a digraph-style matrix of supports for the adjacency Ai-Bj.
    matrix<Bool> supported( nA+nB, nA+nB, False );
    
    
    if ( N_READS_TO_SUPPORT != -1 ) {
      
      // Get the list of all reads falling on each edge.
      vec< vec<int> > A_reads, B_reads;
      for ( int i = 0; i < nA; i++ )
	A_reads.push_back( _edge_to_reads[ A[i] ] );
      for ( int i = 0; i < nB; i++ )
	B_reads.push_back( _edge_to_reads[ B[i] ] );
      
      // Mark which adjacencies are supported by reads.
      for ( int i = 0; i < nA; i++ )
	for ( int j = 0; j < nB; j++ ) {
	  vec<int> common_reads = Intersection( A_reads[i], B_reads[j] );
	  if ( common_reads.isize() >= N_READS_TO_SUPPORT )
	    supported[i][nA+j] = True;
	}
    }
    
    
    if ( N_PAIRS_TO_SUPPORT != -1 ) {
      
      // Get the list of pairs falling on each edge.  This is more restrictive
      // than the list of reads, even if the reads are all paired:
      // AlignPairsToHyper will only place a read pair on the HKP if both reads
      // in the pair are placed with appropriate separation and orientation.
      vec< vec<int> > A_pairs, B_pairs;
      for ( int i = 0; i < nA; i++ )
	A_pairs.push_back( _edge_to_pairs[ A[i] ] );
      for ( int i = 0; i < nB; i++ )
	B_pairs.push_back( _edge_to_pairs[ B[i] ] );
      
      // Mark which adjacencies are supported by paired reads.
      for ( int i = 0; i < nA; i++ )
	for ( int j = 0; j < nB; j++ ) {
	  vec<int> common_pairs = Intersection( A_pairs[i], B_pairs[j] );
	  if ( common_pairs.isize() >= N_PAIRS_TO_SUPPORT )
	    supported[i][nA+j] = True;
	}
    }
    
    
    // Add to tallies.
    total_adjs += nA*nB;
    for ( int i = 0; i < nA; i++ )
      for ( int j = 0; j < nB; j++ )
	if ( supported[i][nA+j] ) total_supported_adjs++;
    
    
    // Split the vertex v into as many separate vertices as necessary, in order
    // to preserve all supported adjacencies while removing as many unsupported
    // adjacencies as possible.
    int n_new_vertices = SplitVertex( v, supported );
    if ( n_new_vertices > 0 ) total_vs_split++;
    total_vs_added += n_new_vertices;
  }
  
  
  // Clean up the graph by merging any edges that do not branch.
  cleanupHKP();
  
  // Final report.
  cout << "\tMarked " << total_supported_adjs << " of " << total_adjs << " edge-edge adjacencies as supported." << endl;
  cout << "\tACTION: Split " << total_vs_split << " vertices, adding " << total_vs_added << " new ones (before cleanup)." << endl;
  cout << Date() << ": END GordianKnot::RemoveUnsupportedAdjacencies!" << endl;
}







/* RemoveNongenomicAdjacencies
 *
 * This function is exactly the same as RemoveUnsupportedAdjacencies (see
 * documentation, above) except that it uses the reference genome to determine
 * which adjacencies are supported.   An adjacency Ai-Bj is considered supported
 * if the bases from Ai and the bases from Bj both align to the reference with a
 * sufficiently low error rate, and the alignments are adjacent.
 *
 * This function is obviously "cheating".  Any assembly produced using
 * RemoveNongenomicAdjacencies should be used for evaluation only.
 *
 *
 * FUNCTION INPUTS
 *
 * -- genome: The reference genome.
 *
 *
 * WARNING: This function can be slow because it calls AlignHyperKmerPath.
 *
 ******************************************************************************/
void
GordianKnot::RemoveNongenomicAdjacencies( const String & GENOME_HEAD )
{
  cout << Date() << ": BEGIN GordianKnot::RemoveNongenomicAdjacencies!" << endl;
  cout << "\tWARNING: The function RemoveNongenomicAdjacencies is a cheating step.\n\tThis output assembly is not truly de novo and should be used for evaluation only." << endl;
  
  // Convert the HyperKmerPath to a HyperBasevector and align it to the
  // reference genome.  This is time-consuming.
  cout << Date( ) << ": Aligning HyperKmerPath to reference" << endl;
  vec<look_align> aligns;
  vec< vec<int> > aligns_index;
  AlignHyperKmerPath( _hkp, _kbb, GENOME_HEAD, _temp_dir, aligns, aligns_index );
  
  cout << Date() << ": Examining each vertex and its adjacencies" << endl;
  int total_adjs = 0, total_supported_adjs = 0;
  int total_vs_split = 0, total_vs_added = 0;
  
  // Loop over every vertex v (that existed at the beginning of this function.)
  for ( int v = 0; v < _n_vertices; v++ ) {
    
    // Find the sets A and B and the numbers nA and nB.
    // (See this function's documentation for an explanation of variable names.)
    vec<int> A = _hkp.ToEdgeObj(v), B = _hkp.FromEdgeObj(v);
    int nA = A.size(), nB = B.size();
    if ( nA == 0 || nB == 0 ) continue;
    
    // Create a digraph-style matrix of supports for the adjacency Ai-Bj.
    matrix<Bool> supported( nA+nB, nA+nB, False );
    
    // Get the list of all alignments-to-reference for each edge.
    vec< vec<look_align> > A_aligns( nA ), B_aligns( nB );
    for ( int i = 0; i < nA; i++ ) {
      const vec<int> & A_index = aligns_index[ A[i] ];
      for ( int j = 0; j < A_index.isize(); j++ )
	A_aligns[i].push_back( aligns[ A_index[j] ] );
    }
    for ( int i = 0; i < nB; i++ ) {
      const vec<int> & B_index = aligns_index[ B[i] ];
      for ( int j = 0; j < B_index.isize(); j++ )
	B_aligns[i].push_back( aligns[ B_index[j] ] );
    }
    
    // Mark which adjacencies are supported by alignments to reference.
    // An adjacency Ai-Bj is considered supported if the bases from Ai and the
    // bases from Bj align to adjacent locations in the reference.
    for ( int i = 0; i < nA; i++ )
      for ( int j = 0; j < nB; j++ ) {
	bool found = false;
	
	// Loop over every alignment of these edges (if they align multiply.)
	for ( int a = 0; a < A_aligns[i].isize(); a++ ) {
	  const look_align & Ai = A_aligns[i][a];
	  for ( int b = 0; b < B_aligns[j].isize(); b++ ) {
	    const look_align & Bj = B_aligns[j][b];
	    
	    // Look for the characteristic (K-1)-base overlap between alignments
	    // which indicates adjacent unibases.
	    if ( Ai.Pos2() == Bj.pos2() + _K-1 )
	      found = true;
	    else if ( Bj.Pos2() == Ai.pos2() + _K-1 )
	      found = true;
	    
	    if ( found ) break;
	  }
	  
	  if ( found ) break;
	}
	
	
	if ( found )
	  supported[i][nA+j] = True;
      }
    
    
    // Add to tallies.
    total_adjs += nA*nB;
    for ( int i = 0; i < nA; i++ )
      for ( int j = 0; j < nB; j++ )
	if ( supported[i][nA+j] ) total_supported_adjs++;
    
    
    // Split the vertex v into as many separate vertices as necessary, in order
    // to preserve all supported adjacencies while removing as many unsupported
    // adjacencies as possible.
    int n_new_vertices = SplitVertex( v, supported );
    if ( n_new_vertices > 0 ) total_vs_split++;
    total_vs_added += n_new_vertices;
  }
  
  // Clean up the graph by merging any edges that do not branch.
  cleanupHKP();
  
  // Final report.
  cout << "\tMarked " << total_supported_adjs << " of " << total_adjs << " edge-edge adjacencies as supported." << endl;
  cout << "\tACTION: Split " << total_vs_split << " vertices, adding " << total_vs_added << " new ones (before cleanup)." << endl;
  cout << Date() << ": END GordianKnot::RemoveNongenomicAdjacencies!" << endl;
}







/* Shave
 *
 *
 * BACKGROUND
 *
 * A "hair" in the HyperKmerPath is an edge E that connects two vertices v and w
 * (in either direction) such that:
 * -- v does not touch any other edges besides E.  Depending on the direction of
 *    E, v is either a source or a sink.
 * -- w has at least one other 'to' edge and one other 'from' edge besides E.
 *    Additionally, the other edges entering/exiting w (in the same direction as
 *    E) must be longer than E.
 *        v
 *        |
 *      E |
 *        |
 * -----> w ----->
 *
 *
 * ALGORITHM STEPS
 *
 * 1. Place the short-insert pairs on the HyperBasevector via AlignPairsToHyper.
 *    Calculate the coverage of each individual edge.
 * 2. Find all edges E that are hairs.  Find the coverage of E, and also the
 *    coverage of all other edges entering and exiting w.  If the former is less
 *    than MIN_HAIR_COV times the latter, remove E.
 * 3. Merge edges in the graph, where possible.
 *
 *
 * FUNCTION INPUTS
 *
 * -- MIN_HAIR_COV: Used in the algorithm step 2, above.  A higher MIN_HAIR_COV
 *    leads to more aggressive shaving.
 *
 ******************************************************************************/
void
GordianKnot::Shave( const double MIN_HAIR_COV )
{
  cout << Date() << ": BEGIN GordianKnot::Shave!" << endl;
  
  
  // This fills _edge_coverage and _C_hbv.
  getCoverageByFragReads( );
  
  
  vec<int> weak_hairs;
  int n_hairs = 0;
  
  cout << Date() << ": Looking for hairs to shave..." << endl;
  if ( _verbose ) cout << "*** Let E be a hair-edge and w be the vertex connecting it to the graph.  E will be removed if its coverage is lower than MIN_HAIR_COV=" << MIN_HAIR_COV << " times the coverage of the other edges entering and exiting w." << endl;
  
  // Find all hairs in the graph.
  for ( int v = 0; v < _n_vertices; v++ ) {
    
    // By definition, every hair-edge E connects vertices v and w, where v is a
    // source/sink, and w connects to the larger graph.
    bool w_E_v;
    int E, w;
    if ( _hkp.FromSize(v) == 0 && _hkp.ToSize(v) == 1 ) {
      E = _hkp.EdgeObjectIndexByIndexTo( v, 0 );
      w = _hkp.To(v)[0];
      w_E_v = true;
    }
    else if ( _hkp.FromSize(v) == 1 && _hkp.ToSize(v) == 0 ) {
      E = _hkp.EdgeObjectIndexByIndexFrom( v, 0 );
      w = _hkp.From(v)[0];
      w_E_v = false;
    }
    else continue;
    
    // Require w to have at least one edge entering it and one edge exiting it,
    // aside from E.
    if ( w_E_v ) {
      if ( _hkp.ToSize(w) < 1 || _hkp.FromSize(w) < 2 ) continue;
    }
    else {
      if ( _hkp.ToSize(w) < 2 || _hkp.FromSize(w) < 1 ) continue;
    }
    
    int E_length = _hkp.EdgeLength(E);
    
    // Find the coverage of all edges (other than E) entering/exiting w.
    // Also require the other edges entering/exiting w, alongside E, to be
    // longer than E.
    bool found_shorter_edge = false;
    int total_cov = 0, total_len = 0;
    for ( int i = 0; i < _hkp.ToSize(w); i++ ) {
      int E1 = _hkp.EdgeObjectIndexByIndexTo( w, i );
      if ( E1 == E ) continue;
      total_cov += _edge_coverage[E1];
      total_len += _hkp.EdgeLength(E1) + _K-1;
      if ( !w_E_v && _hkp.EdgeLength(E1) < E_length )
	found_shorter_edge = true;
    }
    for ( int i = 0; i < _hkp.FromSize(w); i++ ) {
      int E2 = _hkp.EdgeObjectIndexByIndexFrom( w, i );
      if ( E2 == E ) continue;
      total_cov += _edge_coverage[E2];
      total_len += _hkp.EdgeLength(E2) + _K-1;
      if ( w_E_v && _hkp.EdgeLength(E2) < E_length )
	found_shorter_edge = true;
    }
    if ( found_shorter_edge ) continue;
    
    n_hairs++;
    
    // Calculate coverage averages.
    double cov_E = double( _edge_coverage[E] ) / ( E_length + _K-1 );
    double cov_w = double( total_cov ) / total_len;
    
    
    if ( _verbose )
      cout << "Edge " << E << " is a hair with loose end at v=" << v << " and connected end at w=" << w << ". cov_E = " << cov_E << ", cov_w = " << cov_w;
    
    // If E is not sufficiently covered, mark it for deletion.
    if ( cov_E < MIN_HAIR_COV * cov_w ) {
      weak_hairs.push_back( E );
      if ( _verbose ) cout << " -> REMOVED!";
    }
    
    if ( _verbose ) cout << endl;
  }
  
  
  _hkp.DeleteEdges( weak_hairs );
  cleanupHKP();
  
  // Final report.
  cout << "\tACTION: Found " << n_hairs << " hairs; deemed " << weak_hairs.size() << " to be weak and removed them." << endl;
  cout << Date() << ": END GordianKnot::Shave!" << endl;
}



/* PopBubbles
 *
 *
 * BACKGROUND
 *
 * A "bubble" in the HyperKmerPath is a pair (v,w) of graph vertices such that
 * v has exactly two edges leaving it, w has exactly two edges entering it, and
 * they're the same two edges.  Note that v may have any number of edges
 * entering it, and w may have any number of edges leaving it.
 *        _-_
 *  --> v     w -->
 *        -_-
 *
 * To "pop a bubble" is to choose one of the two edges in a bubble and detach
 * it from the graph.
 *
 *
 * ALGORITHM STEPS
 *
 * 1. Place the short-insert pairs on the HyperBasevector via AlignPairsToHyper.
 *    Calculate the coverage of each individual edge.
 * 2. Find bubbles.  In each bubble, choose the edge with lower coverage (C1).
 *    Also find the average coverage (C2) of all edges entering or exiting the
 *    bubble.  If C1 < MIN_BUBBLE_COV * C2, remove the edge, popping the bubble.
 * 3. Merge edges in the graph, where possible.
 *
 *
 * FUNCTION INPUTS
 *
 * -- MIN_BUBBLE_COV: Used in the algorithm step 2, above.  A higher
 *    MIN_BUBBLE_COV leads to more aggressive bubble-popping.
 *
 *
 * WARNINGS
 *
 * Note that bubble-popping is dangerous in diploid assemblies because it can
 * remove true SNPs.  It's a good idea to set MIN_BUBBLE_COV lower in this case.
 *
 ******************************************************************************/
void
GordianKnot::PopBubbles( const double MIN_BUBBLE_COV )
{
  cout << Date() << ": BEGIN GordianKnot::PopBubbles!" << endl;
  
  // This fills _edge_coverage and _C_hbv.
  getCoverageByFragReads( );
  
  cout << Date() << ": Looking for bubbles to pop..." << endl;
  if ( _verbose ) cout << "*** A bubble will only be popped if its lower-coverage edge has coverage less than MIN_BUBBLE_COV=" << MIN_BUBBLE_COV  << " times the 'local coverage' of the edges around the bubble." << endl;
  
  vec<int> bubble_edges;
  int n_bubbles = 0;
  
  // Find all bubbles in the graph.
  for ( int v = 0; v < _n_vertices; v++ ) {
    if ( _hkp.FromSize(v) != 2 ) continue;
    if ( _hkp.From(v)[0] != _hkp.From(v)[1] ) continue;
    int w = _hkp.From(v)[0];
    if ( _hkp.ToSize(w) != 2 ) continue;
    n_bubbles++;
    
    // Calculate the coverage of each edge in the bubble.
    int edge1 = _hkp.FromEdgeObj(v)[0], edge2 = _hkp.FromEdgeObj(v)[1];
    double cov1 = double( _edge_coverage[edge1] ) / ( _hkp.EdgeLength(edge1) + _K-1 );
    double cov2 = double( _edge_coverage[edge2] ) / ( _hkp.EdgeLength(edge2) + _K-1 );
    
    // Find the 'local coverage', which is the average coverage of all edges
    // entering/exiting the bubble (but not in the bubble.)
    vec<int> local_edges = _hkp.ToEdgeObj(v);
    local_edges.append( _hkp.FromEdgeObj(w) );
    int total_cov = 0, total_len = 0;
    for ( int i = 0; i < local_edges.isize(); i++ ) {
      total_cov += _edge_coverage [ local_edges[i] ];
      total_len += _hkp.EdgeLength( local_edges[i] ) + _K-1;
    }
    double local_cov = double( total_cov ) / total_len;
    
    // Find the less-covered edge, and apply the coverage threshold.
    // If it passes, remove the edge.  The graph around the bubble will be
    // cleaned up later by a call to RemoveUnneededVertices.
    double low_cov  = cov1 < cov2 ?  cov1 :  cov2;
    int    low_edge = cov1 < cov2 ? edge1 : edge2;
    double high_cov  = cov1 < cov2 ?  cov2 :  cov1;
    int    high_edge = cov1 < cov2 ? edge2 : edge1;
    
    if ( _verbose )
      cout << setiosflags(ios::fixed) << setprecision(1)
	   << "BUBBLE, vert " << v << " and " << w << ", local cov = "
	   << local_cov << "; low: " << low_edge
           << "(" << BaseAlpha(low_edge) << "), cov = " << low_cov << "; high: "
           << high_edge << "(" << BaseAlpha(high_edge) << "), cov = " << high_cov;
    
    if ( low_cov < MIN_BUBBLE_COV * local_cov ) {
      bubble_edges.push_back( low_edge );
      
      if ( _verbose ) cout << " -> POP!";
    }
    
    if ( _verbose ) cout << endl;
  }
  
  _hkp.DeleteEdges( bubble_edges );
  
  // Clean up the graph by merging any edges that do not branch.
  cleanupHKP();
  
  cout << "\tACTION: Found " << n_bubbles << " bubbles.  Popped " << bubble_edges.size() << " of them." << endl;
  cout << Date() << ": END GordianKnot::PopBubbles!" << endl;
}





/* UnravelLoops
 *
 *
 * BACKGROUND
 *
 * In this context, a "loop" in the HyperKmerPath is an edge E that loops back
 * on itself, going from a vertex v back to the same vertex v.  (We are not yet
 * worrying about multi-edge loops.)  Any non-loop edges entering/exiting v are
 * designated E1/E2; there may be more than one of either of these.  A loop is
 * called "unique" if it is the only loop at vertex v.
 *
 * To "unravel a loop" is to split the edge v into v1 and v2, and place a number
 * of copies (N) of the edge E between v1 and v2.  (The value of N must be
 * determined by the algorithm.)  This can only be done with unique loops.
 *       E
 *       _                     E1       E^N      E2
 *  E1  ( ) E2      becomes   ----> v1 -----> v2 ---->
 * ----> v ---->
 *
 *
 * ALGORITHM STEPS
 *
 * 1. Place the short-insert pairs on the HyperBasevector via AlignPairsToHyper.
 *    Calculate the coverage of each individual edge.
 * 2. Find all loops.  For each unique loop, find the coverage of the loop edge
 *    E, and compare it to the coverages of the other edges E1 entering v and
 *    the other edges E2 leaving v.  If cov(E) is close to a multiple of cov(E1)
 *    and cov(E2), unravel the loop, using the coverage ratio as N.
 * 3. Merge edges in the graph, where possible.
 *
 *
 * FUNCTION INPUTS
 *
 * -- LOOP_COV_TOLERANCE: Used in the algorithm step 2, above.  
 *
 *
 * WARNINGS
 *
 * Non-unique loops cannot be unraveled: If a vertex v contains two loop edges
 * E, UnravelLoops cannot to determine which edge would go first in the newly
 * created edge E^N.
 *
 ******************************************************************************/
void
GordianKnot::UnravelLoops( const double LOOP_COV_TOLERANCE )
{
  cout << Date() << ": BEGIN GordianKnot::UnravelLoops!" << endl;
  
  // This fills _edge_coverage and _C_hbv.
  getCoverageByFragReads( );
  
  // Find all self-looping edges in the HyperKmerPath.
  cout << Date() << ": Find self-looping edges" << endl;
  vec<int> self_loops = _hkp.SelfLoops();
  
  vec<int> to_left, to_right;
  _hkp.ToLeft ( to_left );
  _hkp.ToRight( to_right );
  
  // For each loop, find the edge index E, the vertex index v, and the indices
  // of the non-loop edges E1 entering v and E2 exiting v.
  cout << Date() << ": Finding which loops to unravel" << endl;
  vec< triple<int,int,int> > loops;
  for ( int i = 0; i < self_loops.isize(); i++ ) {
    int E = self_loops[i];
    int v = to_left[E];
    ForceAssertEq( v, to_right[E] );
    vec<int> E1 = _hkp.ToEdgeObj(v), E2 = _hkp.FromEdgeObj(v);
    
    if ( _verbose )
      cout << "Considering loop edge " << E << " on vertex " << v << "... " << flush;
    
    // Filter out non-unique loops (multiple loops at the same vertex.)
    vec<int> loop_edges = Intersection( E1, E2 );
    if ( !loop_edges.solo() ) {
      if ( _verbose ) cout << "not a unique loop." << endl;
      continue;
    }
    ForceAssertEq( loop_edges[0], E );
    
    // Find the average coverage of the edges E1 and E2, and the coverage of
    // the edge E.
    double E12_length = 0, E12_base_cov = 0;
    for ( int i = 0; i < E1.isize(); i++ ) {
      if ( E1[i] == E ) continue;
      E12_base_cov += _edge_coverage[ E1[i] ];
      E12_length += ( _hkp.EdgeLength( E1[i] ) + _K-1 );
    }
    for ( int i = 0; i < E2.isize(); i++ ) {
      if ( E2[i] == E ) continue;
      E12_base_cov += _edge_coverage[ E2[i] ];
      E12_length += ( _hkp.EdgeLength( E2[i] ) + _K-1 );
    }
    if ( E12_length == 0 ) {
      if ( _verbose ) cout << "vertex contains no non-looped edges." << endl;
      continue;
    }
    double E12_cov = E12_base_cov / E12_length;
    double E_cov = double( _edge_coverage[E] ) / double( _hkp.EdgeLength(E) + _K-1 );
    
    // We now have three coverages:
    // (1) _C_hbv: The average coverage of the HyperBasevector.
    // (2) E12_cov: The average coverage of all edges E1, E2.
    // (3) E_cov: The coverage of the loop edge E.
    // Require (2) to be within LOOP_COV_TOLERANCE of a nonzero multiple of (1).
    // Require (3) to be within LOOP_COV_TOLERANCE of a nonzero multiple of (2).
    int ipart = round( E12_cov / _C_hbv );
    if ( ipart == 0 ) {
      if ( _verbose ) cout << "coverage on E1,E2 edges is too low." << endl;
      continue;
    }
    double fpart = fabs( _C_hbv * ipart - E12_cov ) / _C_hbv;
    if ( fpart > LOOP_COV_TOLERANCE ) {
      if ( _verbose ) cout << "coverage on E1,E2 edges does not fall within LOOP_COV_TOLERANCE." << endl;
      continue;
    }
    
    ipart = round( E_cov / E12_cov );
    if ( ipart == 0 ) {
      if ( _verbose ) cout << "coverage on E is too low." << endl;
      continue;
    }
    fpart = fabs( E12_cov * ipart - E_cov ) / E12_cov;
    if ( fpart > LOOP_COV_TOLERANCE ) {
      if ( _verbose ) cout << "coverage on E does not fall within LOOP_COV_TOLERANCE." << endl;
      continue;
    }
    
    // We have an unravelable loop!
    loops.push_back( make_triple( v, E, ipart ) );
    if ( _verbose ) cout << "This is an unravellable loop!  Coverages: _C_hbv = " << _C_hbv << ", E12_cov = " << E12_cov << ", E_cov = " << E_cov << ".  Multiplicity = " << ipart << "." << endl;
  }
  
  
  // Make the new vertices that we'll need to unravel the loops.
  int n_loops = loops.isize();
  _hkp.AddVertices( n_loops );
  
  // Unravel each loop!
  cout << Date() << ": Unraveling " << n_loops << " loops" << endl;
  for ( int i = 0; i < n_loops; i++ ) {
    
    int v = loops[i].first, E = loops[i].second;
    ForceAssert( BinMember( self_loops, E ) );
    // v2 is the newly added vertex that will be used by the unraveled loop.
    int v2 = i + _n_vertices;
    
    // Replace the edge E with the new edge, E^N, consisting of N copies of E.
    KmerPath EN;
    for ( int j = 0; j < loops[i].third; j++ )
      EN.Append( _hkp.EdgeObject( E ) );
    EN.Canonicalize();
    _hkp.SetEdgeObject( E, EN );
    
    // Re-route E^N to v2, making it no longer a loop.
    _hkp.GiveEdgeNewToVx( E, v, v2 );
    
    // Re-route all E2 edges so they come from v2 instead of v.
    vec<int> E2 = _hkp.FromEdgeObj(v);
    for ( int j = 0; j < E2.isize(); j++ ) {
      if ( E2[j] == E ) continue;
      _hkp.GiveEdgeNewFromVx( E2[j], v, v2 );
    }
  }
  
  // Clean up the graph by merging any edges that do not branch.
  cleanupHKP();
  
  cout << "\tACTION: Found " << self_loops.size() << " loops; unraveled " << n_loops << " of them." << endl;
  cout << Date() << ": END GordianKnot::UnravelLoops!" << endl;
}






// Clean up the graph by merging any edges that do not branch.
// Note that this may cause an apparent decrease in reference coverage:
// If two adjacent edges in the HKP align to different locations on the
// reference, and you merge those edges, the resulting edge will have no
// complete alignment.
void
GordianKnot::cleanupHKP( )
{
  _hkp.RemoveUnneededVertices( );
  _hkp.RemoveDeadEdgeObjects( );
  
  // The preceding functions may change edge/vertex functions, so any data
  // structures based on this numbering must be re-loaded.
  _n_vertices = _hkp.N();
  _n_edges = _hkp.EdgeObjectCount();
  _edge_coverage.clear();
}





// Find the coverage of each edge in the HyperBasevector by fragment read pairs.
// The coverage is reported as a total number of bases covered (with each base
// counted once for each time it is covered) rather than as an average.
// This function is time-intensive because it calls AlignPairsToHyper.
void
GordianKnot::getCoverageByFragReads( )
{
  // Don't bother doing this if it's already been done.
  if ( _edge_coverage.nonempty() ) return;
  
  HyperBasevector hbv( _hkp, *_kbb );
  
  // Outputs of AlignPairsToHyper.
  vec<CompressedSeqOnHyper> aligns;
  vec< vec<int> > aligns_index;
  vec< pair<int,int> > pair_places;
  
  // Call AlignPairsToHyper.  Note that the pairs we're aligning are the
  // short-insert pairs, after error correction but before filling fragments.
  AlignPairsToHyper( _hkp, hbv, *_frag_reads, *_frag_pairs,
		     _temp_dir, aligns, aligns_index, pair_places,
		     false, -1, &_cached_edges, &_cached_ppaligns );
  
  
  // Loop over the alignments we've just found, and gather data.
  
  // _edge_to_reads: Index of edge ID to IDs of reads covering that edge.
  // _edge_coverage: Count of the number of bases covering each edge.
  _edge_to_reads.resize( _n_edges, vec<int>() );
  _edge_coverage.resize( _n_edges, 0 );
  for ( int i = 0; i < aligns.isize( ); i++ ) {
    SeqOnHyper a;
    aligns[i].DecompressInto( a, hbv );
    for ( int j = 0; j < a.N(); j++ ) {
      _edge_to_reads[ a.Id2(j) ].push_back( i );
      _edge_coverage[ a.Id2(j) ] += ( a.Pos2(j) - a.pos2(j) );
    }
  }
  
  for ( int i = 0; i < _n_edges; i++ )
    UniqueSort( _edge_to_reads[i] );
  
  // _edge_to_pairs: Index of edge ID to IDs of pairs covering that edge.
  _edge_to_pairs.resize( _n_edges, vec<int>() );
  for ( int i = 0; i < pair_places.isize(); i++ ) {
    SeqOnHyper a1, a2;
    aligns[ pair_places[i].first  ].DecompressInto( a1, hbv );
    aligns[ pair_places[i].second ].DecompressInto( a2, hbv );
    
    for ( int j = 0; j < a1.N(); j++ )
      _edge_to_pairs[ a1.Id2(j) ].push_back(i);
    for ( int j = 0; j < a2.N(); j++ )
      _edge_to_pairs[ a2.Id2(j) ].push_back(i);
  }
  
  for ( int i = 0; i < _n_edges; i++ )
    UniqueSort( _edge_to_pairs[i] );
  
  
  // Find _C_hbv, the average coverage of the HyperBasevector.
  // By counting base coverage rather than kmer coverage, we effectively give
  // extra weight to alignments that bridge graph adjacencies, which is good.
  _C_hbv = double( BigSum( _edge_coverage ) ) / hbv.TotalEdgeLength();
  ForceAssertGt( _C_hbv, 0 );
}
  






// Split the vertex v on _hkp into as many separate vertices as necessary, in
// order to preserve all supported adjacencies while removing as many
// unsupported adjacencies as possible.
// It may be impossible to remove all unsupported adjacencies:
// e.g., if A0-B0, A0-B1, A1-B0 are all supported, the vertex connecting these
// edges will also connect A1-B1 even if it is unsupported.
// Returns the number of vertices created.
// May create new vertices, but will not renumber existing vertices or edges.
int
GordianKnot::SplitVertex( const int v, const matrix<Bool> & supported )
{
  vec<int> A = _hkp.ToEdgeObj(v), B = _hkp.FromEdgeObj(v);
  int nA = A.size(), nB = B.size();
  
  // This digraph has the edges A,B as vertices and contains an edge for every
  // supported adjacency.  Each component in this digraph corresponds to a
  // vertex coming out of the split of v.
  digraph v_supports( supported );
  vec< vec<int> > v_components;
  v_supports.Components( v_components );
  int n_split_vertices = v_components.size();
  if ( n_split_vertices == 1 ) return 0; // can't split this vertex
  
  
  if ( _verbose ) {
    // Print a table illustrating supported adjacencies at this vertex.
    cout << "Splitting vertex " << v << " into " << n_split_vertices << " vertices.  Support of adjacencies is as follows: (X = supported)" << endl;
    vec< vec<String> > table( nA+1, vec<String>(nB+1) );
    for ( int i = 0; i < nA; i++ )
      table[i+1][0] = BaseAlpha( A[i] );
    for ( int j = 0; j < nB; j++ )
      table[0][j+1] = BaseAlpha( B[j] );
    for ( int i = 0; i < nA; i++ )
      for ( int j = 0; j < nB; j++ )
	if ( supported[i][nA+j] )
	  table[i+1][j+1] = "X";
	else
	  table[i+1][j+1] = ".";
    BeautifyAndPrintTable( table, cout );
    cout << endl;
  }
  
  // Add the necessary new vertices and re-route to them the edges in A,B.
  _hkp.AddVertices( n_split_vertices - 1 );
  for ( int i = 1; i < n_split_vertices; i++ ) {
    int new_vertex = _hkp.N() + i - n_split_vertices;
    
    for ( int j = 0; j < v_components[i].isize(); j++ ) {
      if ( v_components[i][j] < nA )
	_hkp.GiveEdgeNewToVx( A[ v_components[i][j] ], v, new_vertex );
      else
	_hkp.GiveEdgeNewFromVx( B[ v_components[i][j] - nA ], v, new_vertex );
    }
  }
  
  return n_split_vertices - 1;
}
