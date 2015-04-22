/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "PrintAlignment.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "dna/Bases.h"
#include "math/Functions.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/MergeScaffoldGraph.h"
#include "paths/ContigsManager.h"
#include "paths/ScaffoldGraphIntegrity.h"

/**
 * PrintInfo
 *
 * Helper function to print core info for an edge between two supers.
 */
void PrintInfoMS( ostream &out,
		  const int v1,
		  const int v2,
		  const CLinkBundle &edge,
		  const vec<superb> &supers )
{
  int s1 = v1 / 2;
  int len1 = supers[s1].TrueLength( );;
  bool rc1 = ( 0 != v1 % 2 );

  int s2 = v2 / 2;
  int len2 = supers[s2].TrueLength( );;
  bool rc2 = ( 0 != v2 % 2 );
  
  out << "s" << s1 << ( rc1 ? "[-]" : "[+]" ) << " (" << len1 << " bp)"
      << "  to  s" << s2 << ( rc2 ? "[-]" : "[+]" ) << " (" << len2 << " bp)"
      << "  :  " << edge.AsString( true )
      << "\n";
}

/**
 * class CMergeEvent
 *
 * Track merge events between contigs, and adjust alignments of reads
 * on merged contigs.
 */
class CMergeEvent {
  
public:
  
  CMergeEvent( ) : cg1_ (-1), rc1_ (false), cg2_ (-1), rc2_ (false) { }
  
  void Set( int cg1, bool rc1, int cg2, bool rc2, alignment &al ) {
    cg1_ = cg1;
    rc1_ = rc1;
    cg2_ = cg2;
    rc2_ = rc2;
    al_ = al;
  }

  bool IsEmpty( ) const { return cg1_ < 0; }
  int SignedContig1( ) const { return ( rc1_ ? -cg1_ - 1 : cg1_) ; }
  int SignedContig2( ) const { return ( rc2_ ? -cg2_ - 1 : cg2_ ); }
  String StrContig1( ) const { return ToString(cg1_) + (rc1_ ? "rc" : "fw"); }
  String StrContig2( ) const { return ToString(cg2_) + (rc2_ ? "rc" : "fw"); }
  int AlignLength( ) const { return al_.Pos1( ) - al_.pos1( ); }
  alignment Alignment( ) const { return al_; }   // return copy
  
  friend ostream &operator<< ( ostream &out, const CMergeEvent &event ) {
    out << "alignment of " << event.AlignLength( )
	<< " bp between " << event.StrContig1( )
	<< " and " << event.StrContig2( );
    return out;
  }
  
  
private:
  
  int cg1_;
  bool rc1_;
  int cg2_;
  bool rc2_;
  alignment al_;
  
};

/**
 * AlignToMerge
 *
 * Align to eventually merge two (singleton) supers. Supers may be
 * oriented fw or rc (orientations are captured by parity of idv and
 * idw), but idw is always assumed to follow idv. It returns a
 * CMergeEvent (which may be empty if AlignToMerge failed).
 */
CMergeEvent AlignToMerge( const CLinkBundle &bundle,
			  const vec<superb> &supers,
			  const vec<fastavector> &contigs,
			  const int idv,
			  const int idw,
			  vec<bool> &touched,
			  ostream &out,
			  const bool VERBOSE,
			  const int MIN_OVERLAP )
{
  CMergeEvent merger;   // default (empty) merge
  
  // HEURISTICS - stretch for gap.
  const double MAX_STRETCH = 6.0;

  // One of both the two contigs have already been merged.
  int s1 = idv / 2;
  int s2 = idw / 2;
  int cg1 = supers[s1].Tig( 0 );
  int cg2 = supers[s2].Tig( 0 );
  if ( touched[cg1] || touched[cg2] ) return merger;

  // Expected overlap.
  int sep = bundle.Sep( );
  int dev = bundle.Dev( );
  int expected_overlap = Min( -sep, (int)contigs[cg2].size( ) );
  if ( expected_overlap  < MIN_OVERLAP ) return merger;

  int radius = (int)( MAX_STRETCH * (double)dev );
  int min_overlap = Max( 0, expected_overlap - radius );
  int max_overlap = expected_overlap + radius;

  // Log.
  if ( VERBOSE ) PrintInfoMS( out, idv, idw, bundle, supers );
  
  // From fastavectors to basevectors.
  GenCharToRandomBaseMapper mapper;
  bool rc1 = ( 1 == idv % 2 );
  bvec b1( contigs[cg1].begin( ), contigs[cg1].end( ), mapper );
  if ( rc1 ) b1.ReverseComplement( );

  bool rc2 = ( 1 == idw % 2 );
  bvec b2( contigs[cg2].begin( ), contigs[cg2].end( ), mapper );
  if ( rc2 ) b2.ReverseComplement( );
       
  // AlignTwoBasevectors arguments.
  ofstream devnull( "/dev/null" );
  int minOv = 0;
  int maxOv = b1.size( ) + b2.size( );
  float maxErrRate = 0.2;
  ostream &a2bLog = devnull;
  int a2bRC = false;
  int mode = 5;
  int K = 24;
  int stretch = 12;
  int nstretch = 6;
  const qvec q1( 0 );
  const qvec q2( 0 );
  float maxScore = 1000000.0;
  float maxErr = 5000;
  ostream &badlog = devnull;
  bool avoidProm = False;
  int maxCliq8 = 1000;
  int maxAligns8 = 10000;
  int maxErr8 = 1000;
  int locMaxErr = 50;
  bool altMethod = false;
  int band = 0;
  int minMutmer = 0;
  float ambigThreshold = 3.0;
  double minPolyscoreRatio = 100.0;
  bool swMethod = False;
  
  if ( VERBOSE ) 
    out << "aligning  c" << cg1 << ( rc1 ? "rc" : "fw" )
	<< " (" << b1.size( )
	<< " bp) to c" << cg2 << ( rc2 ? "rc" : "fw" )
	<< " (" << b2.size( )
	<< " bp)  expected overlap = [" << min_overlap
	<< "," << max_overlap
	<< "]  ";
  
  align theAl;
  AlignTwoBasevectors( b1, b2, theAl, minOv, maxOv, maxErrRate, &a2bLog,
		       a2bRC, mode, K, stretch, nstretch, q1, q2,
		       maxScore, maxErr, badlog, avoidProm, maxCliq8,
		       maxAligns8, maxErr8, locMaxErr, altMethod, band,
		       minMutmer, ambigThreshold, minPolyscoreRatio,
		       swMethod );
  
  alignment al;
  al.Set( theAl, ActualErrors( b1, b2, theAl, 1, 1 ) );
  
  // Is align valid?
  bool is_valid = true;
  
  int al_len = al.Pos1( ) - al.pos1( );
  if ( al_len < 1 ) {
    if ( VERBOSE ) out << "Align not found!\n\n";
    is_valid = false;
  }  
  
  int pos1 = al.pos1( );
  int Pos1 = al.Pos1( );
  int pos2 = al.pos2( );
  int Pos2 = al.Pos2( );
  int len1 = b1.size( );
  int len2 = b2.size( );
  if ( is_valid ) {
    if ( pos2 > 0 || ( Pos1 < len1 && Pos2 < len2 ) ) {
      if ( VERBOSE ) out << "Invalid align found.\n\n";
      is_valid = false;
    }
  }

  if ( is_valid ) {
    if ( al_len < MIN_OVERLAP ) {
      if ( VERBOSE ) out << "Align found, but short (" << al_len << ").\n\n";
      is_valid = false;
    }
  }

  if ( is_valid ) {
    if ( al_len < min_overlap || max_overlap < al_len ) {
      if ( VERBOSE ) out << "Align found, wrong length (" << al_len << ").\n\n";
      is_valid = false;
    }
  }
  
  if ( ! is_valid ) return merger;

  // Align found.
  if ( VERBOSE ) {
    out  << "FOUND!  length = " << al_len << "\n";
    PrintVisualAlignment( True, out, b1, b2, al );
  }
  
  // Good match, set merger and return.
  touched[cg1] = true;
  touched[cg2] = true;
  merger.Set( cg1, rc1, cg2, rc2, al );

  return merger;
  
}
    
/**
 * MergeScaffoldGraph
 */
void MergeScaffoldGraph( vec<superb> &supers,
			 vec<fastavector> &contigs,
			 vec<alignlet> &aligns,
			 vec<int> &index,
			 ostream &out,
			 const PairsManager &pairs,
			 const bool VERBOSE,
			 const int MIN_OVERLAP )
{ 
  out << Date( ) << ": starting MergeScaffoldGraph" << endl;
  
  // Generate trivial supers (needed for BuildScaffoldGrap).
  out << Date( ) << ": generate trivial supers" << endl;

  // Build scaffold graph.
  out << Date( ) << ": building scaffold graph" << endl;
  digraphE<sepdev> unused;
  digraphE<CLinkBundle> graph;
  BuildScaffoldGraph( pairs, supers, aligns, index, unused, &graph );
  ForceAssert( ScaffoldGraphIntegrity( graph ) );
  
  // Maps to nagavigate the initial graph.
  int n_edges = graph.EdgeObjectCount( );
  int n_vertices = graph.N( );
  vec<int> to_left;
  vec<int> to_right;
  graph.ToLeft( to_left );
  graph.ToRight( to_right );
  
  // Sort bundles.
  out << Date( ) << ": sorting bundles" << endl;
  vec<const CLinkBundle*> bundles;
  vec<int> bids;
  bundles.reserve( n_edges/2 );
  bids.reserve( n_edges/2 );
  for (int ii=0; ii<n_edges; ii++) {
    // each edge (bundle) is duplicated for the fw and rc orientations;
    // just take one of the links for sorting and evaluation. --bruce
    if (to_left[ii] < to_right[ii]) {
      bundles.push_back( &( graph.EdgeObject( ii ) ) );
      bids.push_back( ii );
    }
  }
  pCLinkBundle_order_best sorter;
  SortSync( bundles, bids, sorter );

  if ( bundles.isize( ) < 1 ) {
    out << "no bundles found. Leaving now.\n" << endl;
    return;
  }

  // Progress log.
  int dotter = 100;
  out << Date( ) << ": done sorting" << ( VERBOSE ? "\n" : "" ) << endl;
  if ( ! VERBOSE )
    out << Date( ) << ": testing "
	<< bundles.isize( ) << " edges for overlap (.="
	<< dotter  << " edges):\n";
  
  // Keep track of merging events.
  vec<bool> touched( contigs.size( ), false );
  
  // Test each bundle for overlap (in pairs).
  int n_merges = 0;
  CMergeEvent event;
  vec<alignlet> dummy(0);
  vec< vec<int> > dummy_index(0);
  ContigsManager manager( contigs, aligns, index, dummy, dummy_index );
  for (int bundle_id=0; bundle_id<bundles.isize( ); bundle_id++) {
    if ( ! VERBOSE && bundle_id % dotter == 0 )
      Dot( out, bundle_id / dotter );
    
    // nodes
    int idv = to_left[ bids[bundle_id] ];
    int idw = to_right[ bids[bundle_id] ];
    
    // Check next bundle connects idw-rc to idv-rc.
    int idw_rc = ( idw % 2 == 0 ) ? idw + 1 : idw - 1;
    int idv_rc = ( idv % 2 == 0 ) ? idv + 1 : idv - 1;

    // supers
    //int sv = idv / 2, sw = idw / 2;

    // contigs
    //int cw = supers[sw], cv = supers[sv];
    
    const CLinkBundle *pbundle = bundles[bundle_id];

    // Align and merge.
    event = AlignToMerge( *pbundle, supers, contigs, idv, idw, touched,
			  out, VERBOSE, MIN_OVERLAP );
    if ( ! event.IsEmpty( ) ) {
      n_merges++;
      int id1 = event.SignedContig1( );
      int id2 = event.SignedContig2( );
      alignment al = event.Alignment( );
      if (VERBOSE) 
	out << endl << Date( ) << ": merging contigs " << id1 << " & " << id2 << endl;
      manager.MergeContigs( id1, id2, al );
    }
  }

  if ( ! VERBOSE ) out << "\n" << Date() << ": adjusting supers" << endl;
  int sn = 0;
  for (u_int s = 0; s < supers.size(); ++s) {
    int c = supers[s].Tig(0);
    int c_len = contigs[c].size();
    if (c_len > 0 ) {
      supers[sn].Clear();
      supers[sn].PlaceFirstTig(c, c_len);
      ++sn;
    }
  }
  supers.resize(sn);

  // Done.
  out << Date( ) << ": MergeScaffoldGraph done ("
      << n_merges << " merges done)\n"
      << endl;
  
}
