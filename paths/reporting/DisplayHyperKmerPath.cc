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
#include "Basevector.h"
#include "Fastavector.h"
#include "TokenizeString.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperFastavector.h"
#include "paths/HyperKmerPath.h"
#include "util/RunCommand.h"
#include "feudal/BinaryStream.h"

template<class HYPER_T>
void DisplayHyper( const String&, const String&, ostream &out );

template
void DisplayHyper<HyperKmerPath> ( const String&, const String&, ostream& );

template
void DisplayHyper<HyperBasevector> ( const String&, const String&, ostream& );

template
void DisplayHyper<HyperFastavector> ( const String&, const String&, ostream& );

/**
 * DisplayHyperKmerPath
 *
 * Interactive tool to display portions of a given Hyper object
 * (HyperKmerPath, HyperBasevector, or HyperFastavector).  The
 * argument HYPER_T specifies which type of Hyper this is:
 *
 *  HYPER_T = kp       [HyperKmerPath]
 *  HYPER_T = bv       [HyperBasevector]
 *  HYPER_T = fv       [HyperFastavector]
 *
 * Remark: dot must be in your path.
 *
 * HYPER: the Hyper object (full path name)
 * HYPER_T: specify type of Hyper object (see above)
 * ELABELS: full path name to an optional extra_labels file (a vec<String>)
 * VERBOSE: if true, show (in the terminal) the ids of the vertices
 */
int main( int argc, char *argv[] )
{
  
  BeginCommandArguments;
  CommandArgument_String( HYPER );
  CommandArgument_String_OrDefault( HYPER_T, "kp" );
  CommandArgument_String_OrDefault( ELABELS, "" );
  CommandArgument_Bool_OrDefault( VERBOSE, True );
  EndCommandArguments;
    
  if ( StringOfOutput( "which dot" ) == "no" ) {
    cout << "Fatal error: dot is not installed on this system." << endl;
    return 0;
  }
  
  ofstream devnull ( "/dev/null" );
  ostream &out = VERBOSE ? cout : devnull;

  if ( HYPER_T == "kp" )
    DisplayHyper<HyperKmerPath>( HYPER, ELABELS, out );
  else if ( HYPER_T == "bv" )
    DisplayHyper<HyperBasevector>( HYPER, ELABELS, out );
  else if ( HYPER_T == "fv" )
    DisplayHyper<HyperFastavector>( HYPER, ELABELS, out );
  else {
    cout << "Fatal error: " << HYPER_T << " is not a supported.\n" << endl;
    return 0;
  }

  cout << Date( ) << ": done" << endl;

}



/*************************
 * FUNCTIONS START HERE  *
 *************************/



/**
 * BasicStats
 *
 * Basic stats on the given hyper.
 */
template<class HYPER_T>
void BasicStats( const HYPER_T &hyper, ostream &out )
{
  // How many "components" to print.
  const uint max_to_print = 25;

  // Find all components.
  vec< vec<int> > comp;
  hyper.Components( comp );

  // Map number_of_vertices to component_id
  vec< pair<uint,uint> > nv2id;
  nv2id.reserve( comp.size( ) );
  for (uint ii=0; ii<comp.size( ); ii++)
    nv2id.push_back( make_pair( comp[ii].size( ), ii ) );
  sort( nv2id.rbegin( ), nv2id.rend( ) );

  // Print info
  out << "\nBASIC STATS:\n"
      << "  n_components: " << comp.size( ) << "\n"
      << "    n_vertices: " << hyper.N( ) << "\n"
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
 * Generate dot file and view (gv) ps for the given hyper.
 */
template<class HYPER_T>
void DotAndPlot( const HYPER_T &hyper,
		 const vec<String> &elabels,
		 const vec<int> &vertices,
		 ostream &out )
{
  // Dot the sub-digraphE.
  temp_file tmpf = temp_file( "/tmp/DisplayHyperKmerPath.dot_XXXXXX" );
  ofstream dotout( tmpf.c_str( ) );
  const vec<String> *plab = elabels.size( ) < 1 ? 0 : &elabels;
  hyper.PrintSummaryDOT0w( dotout, True, True, True, 0, 0, plab, 0, &vertices );
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
template<class HYPER_T>
void AllPredecessors( const HYPER_T &hyper,
		      vec<int> &data,
		      int &level,
		      int rad,
		      int vtx )
{
  level++;
  data.push_back( vtx );
  if ( level > rad ) return;
  
  vec<int> direct = hyper.To( vtx );
  for (int ii=0; ii<direct.isize( ); ii++)
    AllPredecessors( hyper, data, level, rad, direct[ii] );
}

/**
 * AllSuccessors
 *
 * Recursive function to find all successors of a vertex (these are
 * added to data).
 */
template<class HYPER_T>
void AllSuccessors( const HYPER_T &hyper,
		    vec<int> &data,
		    int &level,
		    int rad,
		    int vtx )
{
  level++;
  data.push_back( vtx );
  if ( level > rad ) return;
  
  vec<int> direct = hyper.From( vtx );
  for (int ii=0; ii<direct.isize( ); ii++)
    AllSuccessors( hyper, data, level, rad, direct[ii] );
}

/**
 * DisplayRadius
 *
 * Warning! If radius is too large, this can explode.
 */
template<class HYPER_T>
void DisplayRadius( const HYPER_T &hyper,
		    const vec<String> &elabels,
		    int vtx,
		    int rad,
		    ostream &out )
{
  // Find all predecessors and all successors for vtx. No memory is reserved!
  vec<int> allvert;
  int level = 0;
  AllPredecessors( hyper, allvert, level, rad, vtx );

  level = 0;
  AllSuccessors( hyper, allvert, level, rad, vtx );

  // Sort unique.
  sort( allvert.begin( ), allvert.end( ) );
  allvert.erase( unique( allvert.begin( ), allvert.end( ) ), allvert.end( ) );

  // Early exit if there are no edges.
  if ( allvert.size( ) < 1 ) {
    out << "No edges around " << vtx << "\n" << endl;
    return;
  }
  out << endl;
  
  // Generate plot.
  DotAndPlot( hyper, elabels, allvert, out );
  
}

/**
 * DisplayComponent
 *
 * Display the component c_id (if c_id is given), or the component
 * that contains the given vertex v_id (if v_id is given).
 */
template<class HYPER_T>
void DisplayComponent( const HYPER_T &hyper,
		       const vec<String> &elabels,
		       ostream &out,
		       int *c_id = 0,
		       int *v_id = 0 )
{
  ForceAssert( c_id || v_id );

  // Find all components.
  vec< vec<int> > comp;
  hyper.Components( comp );
  
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
      out << "Error: " << *v_id << " does not exist\n" << endl;
      return;
    }
  }
  out << endl;

  // Generate plot.
  DotAndPlot( hyper, elabels, comp[comp_id], out );
  
}

/**
 * DisplayRegion
 */
template<class HYPER_T>
void DisplayRegion( const HYPER_T &hyper,
		    const vec<String> &elabels,
		    int v1,
		    int v2,
		    ostream &out )
{
  // Find paths between vertices (each path is a set of vertices).
  vec< vec<int> > paths;
  hyper.AllPaths( v1, v2, paths );

  vec<int> allvert;
  for (int ii=0; ii<paths.isize( ); ii++)
    for (int jj=0; jj<paths[ii].isize( ); jj++)
      allvert.push_back( paths[ii][jj] );
  sort( allvert.begin( ), allvert.end( ) );
  allvert.erase( unique( allvert.begin( ), allvert.end( ) ), allvert.end( ) );

  // Early exit if there are no edges.
  if ( allvert.size( ) < 1 ) {
    out << "No edges between v" << v1 << " and v" << v2 << "\n" << endl;
    return;
  }
  out << endl;

  // Generate plot.
  DotAndPlot( hyper, elabels, allvert, out );
  
}

/**
 * ExecCmnd
 */
template<class HYPER_T>
bool ExecCmnd( const HYPER_T &hyper,
	       const vec<String> &elabels,
	       const vec<int> &to_left,
	       const vec<int> &to_right,
	       const String &cmnd,
	       ostream &out )
{
  // Help.
  if ( cmnd == "h" ) {
    cout << "\n"
	 << "            q: quit\n"
	 << "           bs: basic stats\n"
	 << "        c<id>: show this component\n"
	 << "        v<id>: show component containing this vertex\n"
	 << "        e<id>: show component containing this edge\n"
	 << "    e<id>r<R>: show digraph at this edge, with given radius\n"
	 << "    v<id>r<R>: show digraph at this vertex, with given radius\n"
	 << "e<id1> e<id2>: show digraph between these two edges\n"
	 << "v<id1> v<id2>: show digraph between these two vertices\n"
	 << "\n";
    return true;
  }
  
  // Command string's tokens.
  vec<String> toks;
  Tokenize( cmnd, toks );

  // Basic stats.
  if ( toks.size( ) == 1 && toks[0] == "bs" ) {
    BasicStats( hyper, out );
    return true;
  }
  
  // Display whole component.
  if ( toks.size( ) == 1 && toks[0].Contains( "c", 0 ) ) {
      String str_cid = toks[0].After( "c" );
      if ( ! str_cid.IsInt( ) ) return false;
      int c_id = str_cid.Int( );
      DisplayComponent( hyper, elabels, out, &c_id );
      return true;
  }

  // Display digraph around a given edge (or whole component containing it).
  if ( toks.size( ) == 1 && toks[0].Contains( "e", 0 ) ) {
    String str_edge = toks[0].After( "e" );
    if ( ! str_edge.Contains( "r" ) ) {
      if ( ! str_edge.IsInt( ) ) return false;
      int e_id = str_edge.Int( );
      int v_id = to_left[e_id];
      int w_id = to_right[e_id];
      DisplayComponent( hyper, elabels, out, 0, &v_id );
      return true;
    }
    else {
      String strEdge = str_edge.Before( "r" );
      String strRad = str_edge.After( "r" );
      if ( ! ( strEdge.IsInt( ) || strRad.IsInt( ) ) ) return false;
      int left = to_left[strEdge.Int( )];
      int right = 1 + strRad.Int( );
      DisplayRadius( hyper, elabels, left, right, out );
      return true;
    }
  }

  // Display digraph around a given vertex (or whole component containing it).
  if ( toks.size( ) == 1 && toks[0].Contains( "v", 0 ) ) {
    String str_vid = toks[0].After( "v" );
    if ( ! str_vid.Contains( "r" ) ) {
      if ( ! str_vid.IsInt( ) ) return false;
      int v_id = str_vid.Int( );
      DisplayComponent( hyper, elabels, out, 0, &v_id );
      return true;
    }
    else {
      String strVtx = str_vid.Before( "r" );
      String strRad = str_vid.After( "r" );
      if ( ! ( strVtx.IsInt( ) || strRad.IsInt( ) ) ) return false;
      DisplayRadius( hyper, elabels, strVtx.Int( ), strRad.Int( ), out );
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
      DisplayRegion( hyper, elabels, str1.Int( ), str2.Int( ), out );
      return true;
    }
    if ( ecase ) {
      String str1 = toks[0].After( "e" );
      String str2 = toks[1].After( "e" );
      if ( ! str1.IsInt( ) ) return false;
      if ( ! str2.IsInt( ) ) return false;
      int left = to_left[str1.Int( )];
      int right = to_left[str2.Int( )];
      DisplayRegion( hyper, elabels, left, right, out );
      return true;
    }
  }

  // We should never get here.
  return false;
  
}

/**
 * FetchCmnd
 */
template<class HYPER_T>
bool FetchCmnd( const HYPER_T &hyper, 
		const vec<String> &elabels,
		const vec<int> &to_left,
		const vec<int> &to_right,
		ostream &out,
		istream &in )
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
  if ( ! ExecCmnd( hyper, elabels, to_left, to_right, cmnd, out ) )
    cout << cmnd << ": malformed command. Type h for syntax help\n" << endl;
  
  // Done.
  return true;

}

/**
 * DisplayHyper
 */
template<class HYPER_T>
void DisplayHyper( const String &in_file,
		   const String &elabels_file,
		   ostream &out )
{
  if ( ! IsRegularFile( in_file ) )
    FatalErr( in_file << " not found.\n" << endl );

  out << Date( ) << ": loading hyper object" << endl;
  HYPER_T hyper;
  BinaryReader::readFile( in_file, &hyper );
  
  vec<String> elabels;
  if ( elabels_file != "" ) {
    out << Date( ) << ": loadind extra labels" << endl;
    READX( elabels_file, elabels );
  }

  vec<int> to_left;
  vec<int> to_right;
  hyper.ToLeft( to_left );
  hyper.ToRight( to_right );

  while( FetchCmnd ( hyper, elabels, to_left, to_right, out, cin ) )
    ;
}
