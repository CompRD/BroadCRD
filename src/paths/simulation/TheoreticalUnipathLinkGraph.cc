/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "feudal/BinaryStream.h"
#include "paths/Unipath.h"
#include "paths/UnipathNhood.h"
#include "paths/simulation/Placement.h"

/**
 * TheoreticalUnipathLinkGraph
 *
 * Build the theoretical UnipathLinkGraph of the genomic unipaths. Two
 * vertices (unipaths) are linked by an edge iff they satisfy a few
 * simple criteria (copy number <= MAX_CN; the separation between them
 * is <= MAX_SEP; and both unipaths are >= MIN_LEN kmer long).
 *
 * K: kmer size
 * BASE_DIR: full path name for input/output
 * GENOME: base name for genome files
 * MAX_CN: max copy number
 * MAX_SEP: max sepearation between to unipaths to "link" them
 * MIN_LEN: min length (in kmers) for unipaths
 * WRITE: if false do not save graphs
 * PRINT_LINKS: if true, print readabe version as well
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( BASE_DIR );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault( GENOME, "genome" );
  CommandArgument_UnsignedInt_OrDefault( MAX_CN, 5 );
  CommandArgument_Int_OrDefault( MAX_SEP, 4000 );
  CommandArgument_Int_OrDefault( MIN_LEN, 10 );
  CommandArgument_Bool_OrDefault( WRITE, True );
  CommandArgument_Bool_OrDefault( PRINT_LINKS, True );
  EndCommandArguments;
  
  // File names.
  String strK = ToString( K );

  String placs_file = BASE_DIR + "/" + READS + ".unipaths.k" + strK + ".locs";
  String unipaths_file = BASE_DIR + "/" + READS + ".unipaths.k" + strK;
  String unipathsdb_file = BASE_DIR + "/" + READS + ".unipathsdb.k" + strK;

  String out_base = BASE_DIR + "/" + GENOME;
  String   graph_file = out_base + ".unipathgraph.k" + strK;
  String F_graph_file = out_base + ".unipathgraph_fp.k" + strK;
  String   links_file = out_base + ".unipathgraph_links.k" + strK;
  String F_links_file = out_base + ".unipathgraph_fp_links.k" + strK;
  
  // Load.
  cout << Date( ) << ": loading unipath placements" << endl;
  VecPlacementVec placs;
  placs.ReadAll( placs_file );

  cout << Date( ) << ": map unipaths to their rc's" << endl;
  vec<int> to_rc;
  {
    vec<tagged_rpint> unipathsdb;
    BinaryReader::readFile( unipathsdb_file, &unipathsdb );
    vecKmerPath unipaths(unipaths_file);
    UnipathInvolution( unipaths, unipathsdb, to_rc );
  }
  int n_unipaths = to_rc.size( );
  
  // Turn placements into a flat file.
  cout << Date( ) << ": digesting placements" << endl;
  unsigned int n_placs = 0;
  for (VecPlacementVec::size_type ii=0; ii<placs.size( ); ii++)
    n_placs += placs[ii].size( );
  
  vec< pair<placement,int> > flat2pos;
  flat2pos.reserve( n_placs );
  for (VecPlacementVec::size_type ii=0; ii<placs.size( ); ii++)
    for (PlacementVec::size_type jj=0; jj<placs[ii].size( ); jj++)
      flat2pos.push_back( make_pair( placs[ii][jj], ii ) );

  sort( flat2pos.begin( ), flat2pos.end( ) );

  // Collect edges.
  cout << Date( ) << ": collecting edges" << endl;
  vec< Tedge<int> >      edges;
  vec< Tedge<double> > F_edges;

  for (uint flat_id=0; flat_id<flat2pos.size( ); flat_id++) {
    const int unipath_id = flat2pos[flat_id].second;
    const placement &plac = flat2pos[flat_id].first;
    if ( plac.Rc( ) ) continue;
    if ( plac.Pos( ) - plac.pos( ) < MIN_LEN ) continue;
    if ( placs[unipath_id].size( ) > MAX_CN ) continue;
    
    for (uint jj=flat_id+1; jj<flat2pos.size( ); jj++) {
      const int next_unipath_id = flat2pos[jj].second;
      const placement &next_plac = flat2pos[jj].first;
      const int sep = next_plac.pos( ) - plac.Pos( );
      if ( sep > MAX_SEP ) break;
      if ( next_plac.Rc( ) ) continue;
      if ( next_plac.Pos( ) - next_plac.pos( ) < MIN_LEN ) continue;
      if ( placs[next_unipath_id].size( ) > MAX_CN ) continue;
      
      edges.push( unipath_id, next_unipath_id, sep, 0 );
      edges.push( to_rc[unipath_id], to_rc[next_unipath_id], sep, 0 );
      F_edges.push( unipath_id, next_unipath_id, sep, 0 );
      F_edges.push( to_rc[unipath_id], to_rc[next_unipath_id], sep, 0 );
    }
  }
  
  sort( edges.begin( ), edges.end( ) );
  edges.erase( unique( edges.begin( ), edges.end( ) ), edges.end( ) );
  sort( F_edges.begin( ), F_edges.end( ) );
  F_edges.erase( unique( F_edges.begin( ), F_edges.end( ) ), F_edges.end( ) );

  // Build graphs.
  cout << Date( ) << ": generating plain graphs" << endl;
  digraphE< sepdev>   graph;
  digraphE<fsepdev> F_graph;
  BuildGraphFromEdges(   edges, n_unipaths,   graph );
  BuildGraphFromEdges( F_edges, n_unipaths, F_graph );
  
  // Save graphs.
  if ( WRITE ) {
    cout << Date( ) << ": saving plain graphs" << endl;
    BinaryWriter::writeFile(   graph_file,   graph );
    BinaryWriter::writeFile( F_graph_file, F_graph );
  }

  // Print linking info.
  if ( PRINT_LINKS ) {
    cout << Date( ) << ": saving links info" << endl;
    ofstream   links(   links_file.c_str( ) );
    ofstream F_links( F_links_file.c_str( ) );
    PrintCommandPretty(   links );
    PrintCommandPretty( F_links );

    for (int vv=0; vv<graph.N( ); vv++) {
      for (int ii=0; ii<graph.From( vv ).isize( ); ii++) {
	int ww = graph.From( vv )[ii];
	sepdev    sep =   graph.EdgeObjectByIndexFrom( vv, ii );
	fsepdev F_sep = F_graph.EdgeObjectByIndexFrom( vv, ii );
	links << vv << " --- "
	      << sep.Sep( ) << " +/- "
	      << sep.Dev( ) << " --> "
	      << ww << "\n";
	F_links << vv << " --- "
		<< sep.Sep( ) << " +/- "
		<< sep.Dev( ) << " --> "
		<< ww << "\n";
      }
    }
    links << "\n";
    F_links << "\n";

    links.close( );
    F_links.close( );
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;
}
