/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/reporting/CLinkBundle.h"
#include "util/RunCommand.h"
#include "feudal/BinaryStream.h"

void DisplaySG( const String&in_file, ostream &out );

/**
 * DisplayScaffoldGraph
 *
 * Interactive tool to display portions of a given scaffold graph,
 * given in input as a digraphE<CLinkBundle>.
 *
 * Remark: dot must be in your path.
 */
int main( int argc, char *argv[] )
{
  
  BeginCommandArguments;
  CommandArgument_String( SCAFFOLD_GRAPH );
  CommandArgument_Bool_OrDefault( VERBOSE, True )
  EndCommandArguments;

  if ( StringOfOutput( "which dot" ) == "no" ) {
    cout << "Fatal error: dot is not installed on this system." << endl;
    return 0;
  }
  
  ofstream devnull ( "/dev/null" );
  ostream &out = VERBOSE ? cout : devnull;

  DisplaySG( SCAFFOLD_GRAPH, out );
  
  cout << Date( ) << ": done" << endl;
}



/*************************
 * FUNCTIONS START HERE  *
 *************************/



/**
 * BasicStats
 *
 * Basic stats on the given graph.
 */
void BasicStats( const digraphE<CLinkBundle> &graph, ostream &out )
{
  // How many "components" to print.
  const uint max_to_print = 15;
  
  // Find all components.
  vec< vec<int> > comp;
  graph.Components( comp );

  // Map number_of_vertices to component_id
  vec< pair<uint,uint> > nv2id;
  nv2id.reserve( comp.size( ) );
  for (uint ii=0; ii<comp.size( ); ii++)
    nv2id.push_back( make_pair( comp[ii].size( ), ii ) );
  sort( nv2id.rbegin( ), nv2id.rend( ) );

  // Print info
  out << "\nBASIC STATS:\n"
      << "  n_components: " << comp.size( ) << "\n"
      << "    n_vertices: " << graph.N( ) << "\n"
      << "\n"
      << "COMPONENTS WITH THE MOST VERTICES:\n";
  for (uint ii=0; ii<Min( max_to_print, (uint)nv2id.size( ) ); ii++)
    out << "  " << nv2id[ii].second << "\t" << nv2id[ii].first << "\n";
  if ( nv2id.size( ) > max_to_print ) {
    uint nleft = (uint)nv2id.size( ) - max_to_print;
    out << "  (" << nleft << " other components)\n";
  }
  out << endl;

}

/**
 * DotAndPlot
 *
 * Generate dot file and view (gv) ps of the graph (only selected
 * vertices).
 */
void DotAndPlot( const digraphE<CLinkBundle> &graph,
		 const vec<int> &vselect,
		 ostream &out )
{
  ForceAssert( is_sorted( vselect.begin( ), vselect.end ( ) ) );

  // Dot the sub-digraphE.
  temp_file tmpf = temp_file( "/tmp/DisplayScaffoldGraph.dot_XXXXXX" );
  ofstream dotout( tmpf.c_str( ) );

  // Plot the graph.
  dotout << "digraph G {\n\n"
   	 << "rankdir=LR;\n"
    	 << "node [width=0.1,height=0.1,fontsize=10,shape=point];\n"
     	 << "edge [fontsize=12];\n"
     	 << "rankdir=LR;\n"
     	 << "labeljust=l;\n"
     	 << "color=white;\n";

  for (int vpos=0; vpos<vselect.isize( ); vpos++) {
    int v1 = vselect[vpos];
    bool rc1 = 1 == v1 % 2;
    int id1 = v1 / 2;
    String vlabel1 = ToString( id1 ) + ( rc1 ? "rc" : "fw" );
    dotout << ToString(v1) + " [label=\"" + vlabel1 + "\",shape=circle];" << endl;
    for (int jj=0; jj<graph.From( v1 ).isize( ); jj++) {
      int v2 = graph.From( v1 )[jj];
      if ( ! binary_search( vselect.begin( ), vselect.end( ), v2 ) )
	continue;
      vec<int> edge_ids = graph.EdgesBetween( v1, v2 );
      

      bool rc2 = 1 == v2 % 2;
      int id2 = v2 / 2;
      String vlabel2 = ToString( id2 ) + ( rc2 ? "rc" : "fw" );
      
      String vlabel = vlabel1 + ">" + vlabel2;
      
      vec<String> elabels;
      for (int kk=0; kk<edge_ids.isize( ); kk++) {
	String ename = "E_" + ToString(v1) + "_" + ToString(v2) + "_" + ToString(kk);
	String elabel = graph.EdgeObject( edge_ids[kk] ).AsString( );
	dotout << ename << " [label=\"" + elabel + "\",fontsize=8,shape=box];" << endl;
	dotout << v1 << " -> " << ename << "[arrowhead=none];" << endl;
	dotout << ename << " -> " << v2 << ";" << endl;
      }
    }
  }
  dotout << "\n}" << endl;
  
  dotout.close( );
  
  // Make a ps file out of the dot file, and view it (remove it when done).
  String tmpf_eps = tmpf + ".eps";
  String command( "dot -Tps " + tmpf + " -o " + tmpf_eps );
  RunCommandWithLog( command, "/dev/null" );
  
  command = "( gv " + tmpf_eps + "; rm " + tmpf_eps + " ) &";
  RunCommandWithLog( command, "/dev/null" );
  
}

/**
 * AllPredecessors
 *
 * Recursive function to find all predecessors of a vertex (these are
 * added to data).
 */
void AllPredecessors( const digraphE<CLinkBundle> &graph,
		      vec<int> &data,
		      int &level,
		      int rad,
		      int vtx )
{
  level++;
  data.push_back( vtx );
  if ( level > rad ) return;
  
  vec<int> direct = graph.To( vtx );
  for (int ii=0; ii<direct.isize( ); ii++)
    AllPredecessors( graph, data, level, rad, direct[ii] );
}

/**
 * AllSuccessors
 *
 * Recursive function to find all successors of a vertex (these are
 * added to data).
 */
void AllSuccessors( const digraphE<CLinkBundle> &graph,
		    vec<int> &data,
		    int &level,
		    int rad,
		    int vtx )
{
  level++;
  data.push_back( vtx );
  if ( level > rad ) return;
  
  vec<int> direct = graph.From( vtx );
  for (int ii=0; ii<direct.isize( ); ii++)
    AllSuccessors( graph, data, level, rad, direct[ii] );
}

/**
 * DisplayRadius
 *
 * Warning! If radius is too large, this can explode.
 */
void DisplayRadius( const digraphE<CLinkBundle> &graph,
		    int vtx,
		    int rad,
		    ostream &out )
{
  // Find all predecessors and all successors for vtx. No memory is reserved!
  vec<int> allvert;
  int level = 0;
  AllPredecessors( graph, allvert, level, rad, vtx );
  
  level = 0;
  AllSuccessors( graph, allvert, level, rad, vtx );

  // Sort unique.
  sort( allvert.begin( ), allvert.end( ) );
  allvert.erase( unique( allvert.begin( ), allvert.end( ) ), allvert.end( ) );
  
  // Early exit if there are no edges.
  if ( allvert.size( ) < 1 ) {
    out << "No edges around " << vtx << "\n" << endl;
    return;
  }
  out << endl;
  
  // Plot.
  DotAndPlot( graph, allvert, out );
  
}

/**
 * DisplayComponent
 *
 * Display the component c_id (if c_id is given), or the component
 * that contains the given vertex v_id (if v_id is given). Print the
 * whole graph if *c_id = -1.
 */
void DisplayComponent( const digraphE<CLinkBundle> &graph,
		       ostream &out,
		       int *c_id = 0,
		       int *v_id = 0 )
{
  ForceAssert( c_id || v_id );

  // Find all components.
  vec< vec<int> > comp;
  graph.Components( comp );
  
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
      cout << "Error: " << *v_id << " does not appear in any component" << endl;
      return;
    }
  }

  // Either print the whole graph, or a select component.
  vec<int> allvert;
  if ( comp_id < 0 ) {
    allvert.reserve( graph.N( ) );
    for (int ii=0; ii<graph.N( ); ii++)
      allvert.push_back( ii );
  }
  else {
    allvert.reserve( comp[comp_id].isize( ) );
    for (int ii=0; ii<comp[comp_id].isize( ); ii++)
      allvert.push_back( comp[comp_id][ii] );
    sort( allvert.begin( ), allvert.end( ) );
    allvert.erase( unique( allvert.begin( ), allvert.end( ) ), allvert.end( ) );
  }  

  // Plot.
  DotAndPlot( graph, allvert, out );
  
}

/**
 * DisplayRegion
 */
void DisplayRegion( const digraphE<CLinkBundle> &graph,
		    int v1,
		    int v2,
		    ostream &out )
{
  // Find paths between vertices (each path is a set of vertices).
  vec< vec<int> > paths;
  graph.AllPaths( v1, v2, paths );

  vec<int> allvert;
  for (int ii=0; ii<paths.isize( ); ii++)
    for (int jj=0; jj<paths[ii].isize( ); jj++)
      allvert.push_back( paths[ii][jj] );
  sort( allvert.begin( ), allvert.end( ) );
  allvert.erase( unique( allvert.begin( ), allvert.end( ) ), allvert.end( ) );
  
  // Early exit if there are no edges.
  if ( allvert.size( ) < 1 ) {
    out << "No edges between " << v1 << " and " << v2 << "\n" << endl;
    return;
  }
  out << endl;

  // Plot.
  DotAndPlot( graph, allvert, out );
  
}

/**
 * ExecCmnd
 */
bool ExecCmnd( const digraphE<CLinkBundle> &graph,
	       const String &cmnd,
	       ostream &out )
{
  // Help.
  if ( cmnd == "h" ) {
    cout << "\n"
	 << "            q: quit\n"
	 << "           bs: basic stats\n"
	 << "          all: print the full graph\n"
	 << "        c<id>: show this component\n"
	 << "        v<id>: show component containing this vertex\n"
	 << "    v<id>r<R>: show digraph at this vertex, with given radius\n"
	 << "v<id1> v<id2>: show digraph between these two vertices\n"
	 << "\n";
    return true;
  }
  
  // Command string's tokens.
  vec<String> toks;
  Tokenize( cmnd, toks );
  
  // Basic stats.
  if ( toks.size( ) == 1 && toks[0] == "bs" ) {
    BasicStats( graph, out );
    return true;
  }
  
  // Display whole component.
  if ( toks.size( ) == 1 && toks[0] == "all" ) {
    int c_id = -1;
    DisplayComponent( graph, out, &c_id );
    return true;
  }
  
  // Display whole component.
  if ( toks.size( ) == 1 && toks[0].Contains( "c", 0 ) ) {
    String str_cid = toks[0].After( "c" );
    if ( ! str_cid.IsInt( ) ) return false;
    int c_id = str_cid.Int( );
    DisplayComponent( graph, out, &c_id );
    return true;
  }

  // Display digraph around a given vertex (or whole component containing it).
  if ( toks.size( ) == 1 && toks[0].Contains( "v", 0 ) ) {
    String str_vid = toks[0].After( "v" );
    if ( ! str_vid.Contains( "r" ) ) {
      if ( ! str_vid.IsInt( ) ) return false;
      int v_id = str_vid.Int( );
      DisplayComponent( graph, out, 0, &v_id );
      return true;
    }
    else {
      String strVtx = str_vid.Before( "r" );
      String strRad = str_vid.After( "r" );
      if ( ! ( strVtx.IsInt( ) || strRad.IsInt( ) ) ) return false;
      
      DisplayRadius( graph, strVtx.Int( ), strRad.Int( ), out );
      return true;
    }
  }
  
  // Display component between two vertices.
  if ( toks.size( ) == 2 ) {
    if ( ! toks[0].Contains( "v", 0 ) ) return false;
    if ( ! toks[1].Contains( "v", 0 ) ) return false;
    String str1 = toks[0].After( "v" );
    String str2 = toks[1].After( "v" );
    if ( ! str1.IsInt( ) ) return false;
    if ( ! str2.IsInt( ) ) return false;

    DisplayRegion( graph, str1.Int( ), str2.Int( ), out );
    return true;
  }

  // We should never get here.
  return false;
  
}

/**
 * FetchCmnd
 */
bool FetchCmnd( const digraphE<CLinkBundle> &graph, ostream &out, istream &in )
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
  if ( ! ExecCmnd( graph, cmnd, out ) )
    cout << cmnd << ": malformed command. Type h for syntax help" << endl;
  
  // Done.
  return true;

}

/**
 * DisplaySG
 */
void DisplaySG( const String&in_file, ostream &out )
{
  if ( ! IsRegularFile( in_file ) )
    FatalErr( in_file << " not found.\n" << endl );
  
  out << Date( ) << ": loading digraphE" << endl;
  digraphE<CLinkBundle> graph;
  BinaryReader::readFile( in_file, &graph );
  
  while( FetchCmnd ( graph, out, cin ) )
    ;
}
