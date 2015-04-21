/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Superb.h"
#include "paths/assisted/Backbone.h"
#include "paths/assisted/Proxgraph.h"
#include "paths/RegapSupers.h"
#include "util/RunCommand.h"

/**
 * BuildBackboneScaffolds
 *
 * Load a proximity graph and generate a backbone scaffold.
 *
 * INPUT:
 *   <HEAD_IN>.CProx.digraphE
 *   <HEAD_IN>.fastb
 *
 * OUTPUT:
 *   <HEAD_OUT>.CProx.digraphE
 *   <HEAD_OUT>.contigs.fastb (contigs)
 *   <HEAD_OUT>.superb (supers)
 *   <HEAD_OUT>.orig_ids (map to original ids - signed for orient.)
 *   <HEAD_OUT>.interim_* (if DUMP_INTERIM is true)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( HEAD_IN );
  CommandArgument_String( HEAD_OUT );
  CommandArgument_Bool_OrDefault( DUMP_INTERIM, false );
  CommandArgument_Int_OrDefault( VERBOSITY, 0 );
  EndCommandArguments;

  // Dir and file names.
  String in_proxgraph_file = HEAD_IN + ".CProx.digraphE";
  String in_bases_file = HEAD_IN + ".fastb";
  String in_map_file = HEAD_IN + ".orig.ids";

  String out_proxgraph_file = HEAD_OUT + ".CProx.digraphE";
  String out_bases_file = HEAD_OUT + ".contigs.fastb";
  String out_supers_file = HEAD_OUT + ".superb";
  String out_map_ids_file = HEAD_OUT + ".orig.ids";
  
  vec<String> needed;
  needed.push_back( in_proxgraph_file );
  needed.push_back( in_bases_file );
  needed.push_back( in_map_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  if ( HEAD_IN == HEAD_OUT ) {
    cout << "Fatal error, HEAD_IN and HEAD_OUT must be different.\n" << endl;
    return 1;
  }

  Mkpath( Dirname( HEAD_OUT ) );

  // Load 
  cout << Date( ) << ": loading proxgraph" << endl;
  proxgraph graphE;
  BinaryReader::readFile( in_proxgraph_file, &graphE );
  
  cout << Date( ) << ": loading unibases" << endl;
  vecbvec bases( in_bases_file );

  cout << Date( ) << ": loading orig_ids map" << endl;
  READ( in_map_file, vec<int>, orig_ids );

  // Build initial backbone.
  cout << Date( ) << ": building initial backbone" << endl;
  backbone in_backbone;
  BackboneFromCProx( graphE, in_backbone );
  
  cout << Date( ) << ": removing reference-only links with contrary jump evidence" << endl;
  RemoveContraryRefLinks(in_backbone);

  backbone bblin;
  cout << Date( ) << ": linearizing backbone" << endl;
  BackboneLinearize( in_backbone, bblin );
  
  if ( DUMP_INTERIM ) {
    String out_file = HEAD_OUT + ".interim0.lin";
    PrintBackbone( bblin, out_file, &bases );
  }
  
  // Iterate: triangularize, linearize.
  cout << Date( ) << ": running iterations (one dot per iteration)" << endl;
  int iter = 0;
  while ( 1 ) {
    Dot( cout, ++iter );
    String tri_file = HEAD_OUT + ".interim" + ToString( iter ) + ".tri";
    String lin_file = HEAD_OUT + ".interim" + ToString( iter ) + ".lin";

    backbone bbtri;
    int n_events = BackboneTriangularize( bblin, bbtri );
    if ( DUMP_INTERIM ) PrintBackbone( bbtri, tri_file, &bases );
    
    BackboneLinearize( bbtri, bblin );
    if ( DUMP_INTERIM ) PrintBackbone ( bblin, lin_file, &bases );

    if ( n_events < 1 ) break;
  }
  cout << endl;

  // Convert to proximity graph, and save it.
  cout << Date( ) << ": converting backbone to proxomity graph" << endl;
  proxgraph graphE_out;
  BackboneToCProx( bblin, graphE, graphE_out );
  BinaryWriter::writeFile( out_proxgraph_file, graphE_out );
  
  // Involution maps.
  vec<int> VertRC;
  vec<int> EdgeRC;
  ForceAssert( IsSymmetric( bblin, &VertRC, &EdgeRC ) );

  // Generate and save output contigs and supers.
  cout << Date( ) << ": generating contigs and supers" << endl;
  int n_Vert = bblin.N( );
  int n_supers = 0;
  int n_contigs = 0;
  vec<bool> selected( n_Vert, false );
  for (int ii=0; ii<n_Vert; ii++) {
    const vec<int> &Vert = bblin.Vert( ii );
    if ( Vert.size( ) < 1 ) continue;
    int iiRC = VertRC[ii];
    ForceAssert( ii != iiRC );
    if ( iiRC < ii ) continue;
    
    n_supers += 1;
    n_contigs += Vert.isize( );
    selected[ii] = true;
  }

  vecbvec out_bases;
  vec<int> out_map_ids;
  vec<superb> out_supers;
  out_bases.reserve( n_contigs );
  out_map_ids.reserve( n_contigs );
  out_supers.reserve( n_supers );

  for (int ii=0; ii<n_Vert; ii++) {
    if ( ! selected[ii] ) continue;

    const vec<int> &Vert = bblin.Vert( ii );
    for (int jj=0; jj<Vert.isize( ); jj++) {
      int vert_id = Vert[jj] / 2;
      bool vert_rc = ( Vert[jj] % 2 == 1 );
      int orig_id = orig_ids[vert_id];
      int signed_orig_id = ( vert_rc ? - 1 - orig_id : orig_id );
      const bvec &orig_bases = bases[vert_id];
      bvec rc_orig_bases;
      if ( vert_rc ) {
	rc_orig_bases = orig_bases;
	rc_orig_bases.ReverseComplement( );
      }
      const bvec &out_contig = vert_rc ? rc_orig_bases : orig_bases;
      int out_id = out_bases.size( );
      int out_len = out_contig.size( );

      out_bases.push_back( out_contig );
      out_map_ids.push_back( signed_orig_id );
      if ( jj == 0 ) {
	superb new_super;
	new_super.PlaceFirstTig( out_id, out_len );
	out_supers.push_back( new_super );
      }
      else {
	vec<int> edge_ids = graphE_out.EdgesBetween( Vert[jj-1], Vert[jj] );
	ForceAssert( edge_ids.size( ) == 1 );
	const CProx &edge = graphE_out.EdgeObject( edge_ids[0] );

	pair<int,int> estim_gap = edge.EstimatedGap( );
	int gap = estim_gap.first;
	int dev = estim_gap.second;
	int sid = out_supers.size( ) - 1;
	out_supers[sid].AppendTig( out_id, out_len, gap, dev );
      }
    }
  }

  // Save contigs, supers, and map.
  cout << Date( ) << ": saving contigs, supers, and map" << endl;
  out_bases.WriteAll( out_bases_file );
  WriteSuperbs( out_supers_file, out_supers );
  WRITE( out_map_ids_file, out_map_ids );

  // Done.
  cout << Date( ) << ": BuildBackboneScaffolds done" << endl;
  
}

