///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "TokenizeString.h"
#include "graph/Digraph.h"
#include "paths/assisted/CDisplayProxDigraph.h"
#include "paths/assisted/CProx.h"
#include "util/RunCommand.h"

/**
 * CDisplayProxDigraph
 * Constructor
 */
CDisplayProxDigraph::CDisplayProxDigraph( const proxgraph &graph,
					  const vecbvec &contigs ) :
  graph_ ( graph ),
  contigs_ ( contigs )
{
  this->Setup( );
  this->DefaultArgs( );
}

/**
 * CDisplayProxDigraph
 * GoDisplay
 */
void CDisplayProxDigraph::GoDisplay( ) const
{
  while( this->FetchCommand ( cin ) )
    ;
}

/**
 * CDisplayProxDigraph
 * Setup
 * private
 */
void CDisplayProxDigraph::Setup( )
{
  graph_.ToLeft( to_left_ );
  graph_.ToRight( to_right_ );
}

/**
 * CDisplayProxDigraph
 * DefaultArgs
 * private
 *
 * IMPORTANT: if you want to add an argument, remember to update the
 * following:
 *
 *  - Add a default value in here (and document it!)
 *  - Make sure PrintArgs( ) lists it as a modifiable arg
 *  - List it in the help message in ExecCommand( )
 *  - Allow ExecCommand( ) to set it to a different value
 *
 */
void CDisplayProxDigraph::DefaultArgs( )
{
  // Filter ref-only edges from a vertex v. The idea is to show
  //  ref-only edges between long contigs, or just the ref-only edge
  //  between v and its closest neighbor.
  SMART_REF_FILTER_ = True;

  // Show an edgeonly  if its connected contigs are longer than this.
  MIN_CLEN_ = 1000;
}

/**
 * CDisplayProxDigraph
 * PrintArgs
 * private
 */
void CDisplayProxDigraph::PrintArgs( ) const
{
  cout << "\n"
       << " === CURRENT VALUES FOR ARGUMENTS ===\n"
       << "\n"
       << " SMART_REF_FILTER: " << ToString( SMART_REF_FILTER_ ) << "\n"
       << "         MIN_CLEN: " << ToString( MIN_CLEN_ ) << "\n"
       << "\n";
}

/**
 * CDisplayProxDigraph
 * FetchCommand
 * private
 */
bool CDisplayProxDigraph::FetchCommand( istream &in ) const
{
  cout << "> " << flush;
  
  // Wait for command.
  String cmnd;
  getline( cin, cmnd );
  if ( ! cin ) {
    cout << endl;
    return false;
  }

  // "q" (quit) is treated separately.
  if ( cmnd == "q" )
    return false;
  
  // Try to run command.
  if ( ! this->ExecCommand( cmnd ) )
    cout << cmnd << ": malformed command. Type h for syntax help\n" << endl;
  
  // Done.
  return true;
}

/**
 * CDisplayProxDigraph
 * private
 * ExecCommand
 */
bool CDisplayProxDigraph::ExecCommand( const String &cmnd ) const
{
  // Help.
  if ( cmnd == "h" ) {
    cout << "\n"
	 << " === BASIC COMMANDS ===\n"
	 << "\n"
	 << "            q: quit\n"
	 << "           bs: basic stats\n"
	 << "         args: print current values for args\n"
	 << "        c<id>: show this component\n"
	 << "        v<id>: show component containing this vertex\n"
	 << "        e<id>: show component containing this edge\n"
	 << "    e<id>r<R>: show digraph at this edge, with given radius\n"
	 << "    v<id>r<R>: show digraph at this vertex, with given radius\n"
	 << "e<id1> e<id2>: show digraph between these two edges\n"
	 << "v<id1> v<id2>: show digraph between these two vertices\n"
	 << "\n"
	 << " === OPTIONAL ARGUMENTS ===\n"
	 << "\n"
	 << "set SMART_REF_FILTER (true or false - default is true):\n"
	 << "  apply a smart filter to remove some of the ref-only edges\n"
	 << "set MIN_CLEN (int - default is 1000):\n"
	 << "  only show edges between long contigs\n"
	 << "\n";
    return true;
  }
  
  // Command string's tokens.
  vec<String> toks;
  Tokenize( cmnd, toks );

  // An optional argument (set).
  if ( toks.size( ) > 0 && toks[0] == "set" ) {
    if ( toks.size( ) == 3 && toks[1] == "SMART_REF_FILTER" ) {
      bool ok = ( toks[2] == "true" || toks[2] == "false" );
      SMART_REF_FILTER_ = ( toks[2] == "true" ? true : false );
      return ok;
    }
    if ( toks.size( ) == 3 && toks[1] == "MIN_CLEN" ) {
      MIN_CLEN_ = toks[2].Int( );
      return true;
    }
    return false;
  }

  // Basic stats.
  if ( toks.size( ) == 1 && toks[0] == "bs" ) {
    this->BasicStats( );
    return true;
  }
  
  // Print args.
  if ( toks.size( ) == 1 && toks[0] == "args" ) {
    this->PrintArgs( );
    return true;
  }

  // Display whole component.
  if ( toks.size( ) == 1 && toks[0].Contains( "c", 0 ) ) {
      String str_cid = toks[0].After( "c" );
      if ( ! str_cid.IsInt( ) ) return false;
      int c_id = str_cid.Int( );
      this->DisplayComponent( &c_id );
      return true;
  }

  // Display digraph around a given edge (or whole component containing it).
  if ( toks.size( ) == 1 && toks[0].Contains( "e", 0 ) ) {
    String str_edge = toks[0].After( "e" );
    if ( ! str_edge.Contains( "r" ) ) {
      if ( ! str_edge.IsInt( ) ) return false;
      int e_id = str_edge.Int( );
      int v_id = to_left_[e_id];
      int w_id = to_right_[e_id];
      this->DisplayComponent( 0, &v_id );
      return true;
    }
    else {
      String strEdge = str_edge.Before( "r" );
      String strRad = str_edge.After( "r" );
      if ( ! ( strEdge.IsInt( ) || strRad.IsInt( ) ) ) return false;
      int left = to_left_[strEdge.Int( )];
      int right = 1 + strRad.Int( );
      this->DisplayRadius( left, right );
      return true;
    }
  }

  // Display digraph around a given vertex (or whole component containing it).
  if ( toks.size( ) == 1 && toks[0].Contains( "v", 0 ) ) {
    String str_vid = toks[0].After( "v" );
    if ( ! str_vid.Contains( "r" ) ) {
      if ( ! str_vid.IsInt( ) ) return false;
      int v_id = str_vid.Int( );
      this->DisplayComponent( 0, &v_id );
      return true;
    }
    else {
      String strVtx = str_vid.Before( "r" );
      String strRad = str_vid.After( "r" );
      if ( ! ( strVtx.IsInt( ) || strRad.IsInt( ) ) ) return false;
      DisplayRadius( strVtx.Int( ), strRad.Int( ) );
      return true;
    }
  }
  
  // Display component between two vertices or edges.
  if ( toks.size( ) == 2 ) {
    bool vcase = ( toks[0].Contains( "v", 0 ) && toks[1].Contains( "v", 0 ) );
    bool ecase = ( toks[0].Contains( "e", 0 ) && toks[1].Contains( "e", 0 ) );
    if ( ! ( vcase || ecase ) ) return false;

    if ( vcase ) {
      String str1 = toks[0].After( "v" );
      String str2 = toks[1].After( "v" );
      if ( ! str1.IsInt( ) ) return false;
      if ( ! str2.IsInt( ) ) return false;
      this->DisplayRegion( str1.Int( ), str2.Int( ) );
      return true;
    }
    if ( ecase ) {
      String str1 = toks[0].After( "e" );
      String str2 = toks[1].After( "e" );
      if ( ! str1.IsInt( ) ) return false;
      if ( ! str2.IsInt( ) ) return false;
      int left = to_left_[str1.Int( )];
      int right = to_left_[str2.Int( )];
      this->DisplayRegion( left, right );
      return true;
    }
  }

  // We should never get here.
  return false;
  
}

/**
 * CDisplayProxDigraph
 * AllPredecessors
 * private
 *
 * Recursive function to find all predecessors of a vertex (these are
 * added to data).
 */
void CDisplayProxDigraph::AllPredecessors( vec<int> &data,
					   int &level,
					   int rad,
					   int vtx ) const
{
  level++;
  data.push_back( vtx );
  if ( level > rad ) return;
  
  vec<int> direct = graph_.To( vtx );
  for (int ii=0; ii<direct.isize( ); ii++)
    this->AllPredecessors( data, level, rad, direct[ii] );
}

/**
 * CDisplayProxDigraph
 * AllSuccessors
 * private
 * 
 * Recursive function to find all successors of a vertex (these are
 * added to data).
 */
void CDisplayProxDigraph::AllSuccessors( vec<int> &data,
					 int &level,
					 int rad,
					 int vtx ) const
{
  level++;
  data.push_back( vtx );
  if ( level > rad ) return;
  
  vec<int> direct = graph_.From( vtx );
  for (int ii=0; ii<direct.isize( ); ii++)
    this->AllSuccessors( data, level, rad, direct[ii] );
}

/**
 * CDisplayProxDigraph
 * ContigLengthFilter
 * private
 *
 * Keep only edges between contigs longer than or equal to MIN_CLEN_.
 */
void CDisplayProxDigraph::ContigLengthFilter( vec<int> &edges ) const
{
  vec<int> keepers;
  keepers.reserve( edges.size( ) );
  for (int ii=0; ii<edges.isize( ); ii++) {
    int edge_id = edges[ii];
    const CProx &prox = graph_.EdgeObject( edges[ii] );
    int vlen = contigs_[ to_left_[edge_id] / 2 ].size( );
    int wlen = contigs_[ to_right_[edge_id] / 2 ].size( );
    if ( vlen >= MIN_CLEN_ && wlen >= MIN_CLEN_ )
      keepers.push_back( edges[ii] );
  }

  swap( edges, keepers );
}

/**
 * CDisplayProxDigraph
 * SmartRefFilter
 * private
 *
 * Filter edges from a single vertex by removing all but a few of the
 * ref-only ones (edges are assumed to have the same to_left_, but
 * this is not checked). Remark: edges will be overwritten.
 */
void CDisplayProxDigraph::SmartRefFilter( vec<int> &edges ) const
{
  // HEURISTICS.
  const int MIN_REF_FILTER_CLEN = 2000;   // defines a long contig

  // This should never happen.
  if ( edges.size( ) < 1 ) return;
  
  // Tag ref and ref-only edges.
  vec<bool> rtag( edges.size( ), false );
  for (int ii=0; ii<edges.isize( ); ii++)
    if ( graph_.EdgeObject( edges[ii] ).RefGap( ) != INT_MAX )
      rtag[ii] = true;

  // Find the shortest ref edge.
  int best_id = -1;
  for (int ii=0; ii<edges.isize( ); ii++) {
    const CProx &prox = graph_.EdgeObject( edges[ii] );
    if ( prox.RefGap( ) == INT_MAX ) continue;
    if ( best_id < 0 ) {
      best_id = ii;
      continue;
    }
    int best_gap = graph_.EdgeObject( edges[best_id] ).RefGap( );
    int this_gap = prox.RefGap( );
    if ( this_gap < best_gap ) {
      best_id = ii;
      continue;
    }
  }

  // The keepers.
  vec<int> keepers;
  keepers.reserve( edges.size( ) );
  for (int ii=0; ii<edges.isize( ); ii++) {
    int edge_id = edges[ii];
    const CProx &prox = graph_.EdgeObject( edges[ii] );

    // Keep edge if it is supported by links.
    if ( prox.NLinksWinner( ) > 0 ) {
      keepers.push_back( edges[ii] );
      continue;
    }

    // Ref-only edge, keep it if it connects long contigs.
    int vlen = contigs_[ to_left_[edge_id] / 2 ].size( );
    int wlen = contigs_[ to_right_[edge_id] / 2 ].size( );
    if ( vlen >= MIN_REF_FILTER_CLEN && wlen >= MIN_REF_FILTER_CLEN ) {
      keepers.push_back( edges[ii] );
      continue;
    }

    // Keep the ref-only edge implying the smallest gap.
    if ( ii == best_id ) {
      keepers.push_back( edges[ii] );
      continue;
    }
  }

  // Done.
  swap( keepers, edges );
}

/**
 * CDisplayProxDigraph
 * BasicStats
 * private
 */
void CDisplayProxDigraph::BasicStats( ) const
{
  // How many "components" to print.
  const uint max_to_print = 12;

  // Find all components.
  vec< vec<int> > comp;
  graph_.Components( comp );

  // Map number_of_vertices to component_id
  vec< pair<uint,uint> > nv2id;
  nv2id.reserve( comp.size( ) );
  for (uint ii=0; ii<comp.size( ); ii++)
    nv2id.push_back( make_pair( comp[ii].size( ), ii ) );
  sort( nv2id.rbegin( ), nv2id.rend( ) );

  // Print info
  cout << "\nBASIC STATS:\n"
       << "  n_components: " << comp.size( ) << "\n"
       << "    n_vertices: " << graph_.N( ) << "\n"
       << "\n"
       << "COMPONENTS WITH THE MOST VERTICES:\n";
  for (uint ii=0; ii<Min( max_to_print, (uint)nv2id.size( ) ); ii++)
    cout << "  " << nv2id[ii].second << "\t" << nv2id[ii].first << "\n";
  if ( nv2id.size( ) > max_to_print ) {
    uint nleft = (uint)nv2id.size( ) - max_to_print;
    cout << "  (" << nleft << " other components)\n";
  }
  cout << endl;
  
}

/**
 * CDisplayProxDigraph
 * DisplayRadius
 * private
 *
 * Warning! If radius is too large, this can explode.
 */
void CDisplayProxDigraph::DisplayRadius( int vtx, int rad ) const
{
  // Find all predecessors and all successors for vtx. No memory is reserved!
  vec<int> allvert;
  int level = 0;
  this->AllPredecessors( allvert, level, rad, vtx );

  level = 0;
  this->AllSuccessors( allvert, level, rad, vtx );

  // Sort unique.
  sort( allvert.begin( ), allvert.end( ) );
  allvert.erase( unique( allvert.begin( ), allvert.end( ) ), allvert.end( ) );

  // Early exit if there are no edges.
  if ( allvert.size( ) < 1 ) {
    cout << "No edges around " << vtx << "\n" << endl;
    return;
  }
  cout << endl;
  
  // Generate plot.
  this->DotAndPlot( allvert );
  
}

/**
 * CDisplayProxDigraph
 * DisplayComponent
 * private
 *
 * Display the component c_id (if c_id is given), or the component
 * that contains the given vertex v_id (if v_id is given).
 */
void CDisplayProxDigraph::DisplayComponent( int *c_id, int *v_id ) const
{
  ForceAssert( c_id || v_id );

  // Find all components.
  vec< vec<int> > comp;
  graph_.Components( comp );
  
  // Find component to print.
  int comp_id = -1;
  if ( c_id )
    comp_id = *c_id;
  else {
    ForceAssert( v_id );
    for (int ii=0; ii<comp.isize( ); ii++) {
      if ( binary_search( comp[ii].begin( ), comp[ii].end( ), *v_id ) ) {
	comp_id = ii;
	break;
      }
    }
    if ( comp_id < 0 ) {
      cout << "Error: " << *v_id << " does not exist\n" << endl;
      return;
    }
  }
  cout << endl;

  // Generate plot.
  DotAndPlot( comp[comp_id] );
  
}

/**
 * CDisplayProxDigraph
 * DisplayRegion
 * private
 */
void CDisplayProxDigraph::DisplayRegion( int v1, int v2 ) const
{
  // Find paths between vertices (each path is a set of vertices).
  vec< vec<int> > paths;
  graph_.AllPaths( v1, v2, paths );

  vec<int> allvert;
  for (int ii=0; ii<paths.isize( ); ii++)
    for (int jj=0; jj<paths[ii].isize( ); jj++)
      allvert.push_back( paths[ii][jj] );
  sort( allvert.begin( ), allvert.end( ) );
  allvert.erase( unique( allvert.begin( ), allvert.end( ) ), allvert.end( ) );

  // Early exit if there are no edges.
  if ( allvert.size( ) < 1 ) {
    cout << "No edges between v" << v1 << " and v" << v2 << "\n" << endl;
    return;
  }
  cout << endl;

  // Generate plot.
  DotAndPlot( allvert );
  
}

/**
 * CDisplayProxDigraph
 * DotAndPlot
 * private
 *
 * Generate dot file and view (gv) ps for the given graph. NB:
 * vertices must be sorted, but for efficiency we do not check
 * this.
 */
void CDisplayProxDigraph::DotAndPlot( const vec<int> &vertices ) const
{
  // Open dot stream.
  temp_file tmpf = temp_file( "/tmp/DisplayCProxDigraphE.dot_XXXXXX" );
  ofstream dotout( tmpf.c_str( ) );
  dotout << "digraph G {\n"
	 << "\n"
	 << "rankdir=LR;\n"
	 << "node [width=0.1, height=0.1, fontsize=10, shape=ellipse];\n"
	 << "edge [fontsize=12];\n"
	 << "\n";

  // List of edges and vertices to print.
  vec<int> e_ids;
  vec<int> v_ids;

  // Select edges to print.
  for (int ii=0; ii<vertices.isize( ); ii++) {
    int vid = vertices[ii];
    vec<int> all_edges = graph_.FromEdgeObj( vid );
    
    // Only keep edges between the input set of vertices.
    vec<int> edges;
    edges.reserve( all_edges.size( ) );
    for (int jj=0; jj<all_edges.isize( ); jj++) {
      int wid = to_right_[ all_edges[jj] ];
      if ( binary_search( vertices.begin( ), vertices.end( ), wid  ) )
	edges.push_back( all_edges[jj] );
    }
    
    // Discard edges between short contigs.
    this->ContigLengthFilter( edges );

    // Run SmartRefFilter, if required.
    if ( SMART_REF_FILTER_ ) this->SmartRefFilter( edges );

    // Add edges to list.
    copy( edges.begin( ), edges.end( ), back_inserter( e_ids ) );
  }
  
  // Select vertices to print.
  v_ids.reserve( 2 * e_ids.size( ) );
  for (int ii=0; ii<e_ids.isize( ); ii++) {
    v_ids.push_back( to_left_[ e_ids[ii] ] );
    v_ids.push_back( to_right_[ e_ids[ii] ] );
  }
  sort( v_ids.begin( ), v_ids.end( ) );
  v_ids.erase( unique( v_ids.begin( ), v_ids.end( ) ), v_ids.end( ) );

  // Print edges.
  for (int ii=0; ii<e_ids.isize( ); ii++) {
    int eid = e_ids[ii];
    int vid = to_left_[eid];
    int wid = to_right_[eid];
    const CProx &prox = graph_.EdgeObject( eid );
    String text = ToString(prox.NLinksWinner());
    CLinkBundle clb = prox.Bundle();
    if (clb.Weight() > 0)
      text = ToString(clb.Weight());
    text += "/" + ToString(prox.EstimatedGap().first);
    String color = ( prox.RefGap( ) == INT_MAX ) ? "black" : "red";
    String label
      = "label = \"" + text + "\", "
      + String( "color = " + color + ", " )
      + String( "style = bold, " );
    dotout << vid << " -> " << wid << " [" << label << "]\n";
  }
  dotout << "\n";

  // Print vertices.
  for (int ii=0; ii<v_ids.isize( ); ii++) {
    int vid = v_ids[ii];
    bool rc = ( vid % 2 != 0 );
    int vlen = contigs_[ vid / 2 ].size( );
    String color;
    if ( vlen < 100 ) color = "gray";
    else if ( vlen < 1000 ) color = "black";
    else if ( vlen < 10000 ) color = "red";
    else color = "magenta";
    String label
      = "label = \"" + ToString( vid / 2 ) + ( rc ? "[-]" : "[+]") + ToString(vlen) + "\", " 
      + String( "style = bold, " )
      + String( "color = " ) + color;
    dotout << vid << " [" << label << "];\n"; 
  }
  dotout << "\n";

  // Close dot stream.
  dotout << "}\n";
  dotout.close( );

  // Make a ps file out of the dot file, and view it (remove it when done).
  String tmpf_eps = tmpf + ".eps";
  String command( "dot -Tps " + tmpf + " -o " + tmpf_eps );
  RunCommandWithLog( command, "/dev/null" );
  
  command = "( gv " + tmpf_eps + "; rm " + tmpf_eps + " ) &";
  RunCommandWithLog( command, "/dev/null" );
  
}

