///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* GORDIAN KNOT: From the name of a legendary knot tied to a pole near the
 * temple of Zeus in Gordium.  It was prophesied that whoever loosed the knot
 * would become ruler of all Asia.  Alexander the Great solved the puzzle by
 * slicing through the knot and took it as a sign of Zeus's favor.  He then
 * proceeded to conquer much of the known world.
 *         -- Wiktionary
 *
 *
 *
 * The GordianKnot module is where we attempt to disentangle the knottiness in
 * an assembly produced via the RunAllPathsLG pipeline.  The wrapper module
 * is in CutTheGordianKnot.cc.  If we can do this, we can do anything!
 *
 * The GordianKnot class is simply a struct that contains a HyperKmerPath as
 * well as several additional (const) data structures used in its algorithms.
 * The class functions implement algorithms that modify the HyperKmerPath.  Each
 * algorithm detaches from the HyperKmerPath certain edges that are deemed
 * spurious.  The GetHKP function allows access to the modified HyperKmerPath,
 * with optional cleanup.
 *
 *
 * FUNCTIONS (see GordianKnot.cc for documentation of the algorithms)
 * -- PopBubbles
 *
 *
 * Josh Burton
 * November 2009
 *
 ******************************************************************************/

#ifndef GORDIAN_KNOT_H
#define GORDIAN_KNOT_H

#include "Basevector.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "math/Matrix.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/SeqOnHyper.h"





// The GordianKnot class is simply a struct that contains a HyperKmerPath as
// well as several additional data structures used in its algorithms.
// The class functions modify the HyperKmerPath without affecting the auxiliary
// structures.
class GordianKnot
{
public:
  
  GordianKnot( const HyperKmerPath & hkp,
	       const vecbasevector * frag_reads,
	       const PairsManager * frag_pairs,
	       const KmerBaseBroker * kbb,
	       const String temp_dir,
	       const bool verbose = true )
    : _hkp( hkp ),
      _n_edges( hkp.EdgeObjectCount() ),
      _n_vertices( hkp.N() ),
      _K( hkp.K() ),
      _frag_reads( frag_reads ),
      _frag_pairs( frag_pairs ),
      _kbb( kbb ),
      _temp_dir( temp_dir ),
      _verbose( verbose )
  { }
  
  // GetHKP: This is the only way to access the output of the GordianKnot.
  HyperKmerPath GetHKP( ) const { return _hkp; }
  
  /*****************************************************************************
   *
   *      ALGORITHMS
   *      These functions are where the magic happens.
   *      See GordianKnot.cc for a thorough documentation of these algorithms.
   *      Each of these functions calls cleanupHKP() at the end.
   *
   ****************************************************************************/
  void RemoveUnsupportedAdjacencies( const int N_READS_TO_SUPPORT, const int N_PAIRS_TO_SUPPORT );
  void Shave( const double MIN_HAIR_COV );
  void PopBubbles( const double MIN_BUBBLE_COV );
  void UnravelLoops( const double LOOP_COV_TOLERANCE );
  
  // The following algorithms "cheat" by using the reference genome to modify
  // the assembly.  Any assembly produced using these functions should be used
  // for evaluation only.
  void RemoveNongenomicAdjacencies( const String & GENOME_HEAD );
  
private:
  
  // Clean up the graph by merging any edges that do not branch.  Note that
  // this may change edge/vertex numbering, so any data structures based on that
  // numbering are invalidated and will be cleared by this function.
  void cleanupHKP( );
  
  // Fill the data structures _edge_coverage, _edge_to_reads, etc.
  void getCoverageByFragReads( );
  
  // Split the vertex v on _hkp into as many separate vertices as necessary, in
  // order to preserve all supported adjacencies while removing as many
  // unsupported adjacencies as possible.
  // It may be impossible to remove all unsupported adjacencies: e.g., if A0-B0,
  // A0-B1, A1-B0 are all supported, the vertex connecting these edges will also
  // connect A1-B1 even if it is unsupported.
  // May create new vertices, but will not renumber existing vertices or edges.
  int SplitVertex( const int v, const matrix<Bool> & supported );
  
  /*****************************************************************************
   *
   *      DATA STRUCTURES
   *
   ****************************************************************************/
  
  // The HyperKmerPath that we're trying to disentangle!
  HyperKmerPath _hkp;
  int _n_edges, _n_vertices;
  const int _K;
  
  // Large data structures.
  // These are all constant, and they are not owned by the GordianKnot class.
  const vecbasevector * _frag_reads;
  const PairsManager * _frag_pairs;
  const KmerBaseBroker * _kbb;
  
  // Derived data structures.  These are not filled with data until necessary.
  vec<int> _edge_coverage; // coverage of edges in HBV by _frag_reads
  vec< vec<int> > _edge_to_reads; // list of reads covering each edge
  vec< vec<int> > _edge_to_pairs; // list of read pairs covering each edge
  double _C_hbv; // Average coverage of the HyperBasevector
  
  // Cached data.
  vecbasevector _cached_edges;
  vec<partial_perfect_align> _cached_ppaligns;
  
  // Miscellany.
  const String _temp_dir;
  const bool _verbose;
};




#endif
