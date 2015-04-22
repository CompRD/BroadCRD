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
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "paths/Alignlet.h"
#include "paths/reporting/CSuperLinks.h"

/**
 * PrintRawSeparations
 *
 * Print all implied separations (gap sizes) between pairs of contigs
 * (output is sent to cout).
 *
 * READS: it loads <READS>.pairs
 * ALIGNS: it loads <ALIGNS>.{qltoutlet,index}
 * CONTIGS: it loads <CONTIGS>.fasta
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( ALIGNS );
  CommandArgument_String( CONTIGS );
  EndCommandArguments;
  
  String contigs_file = CONTIGS + ".fasta";
  String pairs_file = READS + ".pairs";
  String aligns_file = ALIGNS + ".qltoutlet";
  String index_file = ALIGNS + ".index";

  // Load.
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  int n_tigs = contigs.size( );

  PairsManager pairs( pairs_file );

  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );
  
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  // Build singleton supers, and a linker.
  vec<superb> supers;
  supers.resize( n_tigs );
  for (int tig=0; tig<n_tigs; tig++) {
    supers[tig].SetNtigs( 1 );
    supers[tig].SetTig( 0, tig );
    supers[tig].SetLen( 0, contigs[tig].size( ) );
  }
  
  CSuperLinks linker( &pairs, &supers, &aligns, &index );
  
  // Legend.
  cout << "cg1 "
       << "cg2 "
       << "rc2 "
       << "cluster_id "
       << "lib_id "
       << "gap\n"
       << "\n";

  // Loop over all contigs.
  for (int tig_id=0; tig_id<n_tigs; tig_id++) {
    float slop = -1.0;
    float stretch = 1.0;
    vec<COffset> all_offsets;
    linker.AllLinks( tig_id, all_offsets, slop, &stretch );

    for (int set_id=0; set_id<all_offsets.isize( ); set_id++) {
      const COffset &offset = all_offsets[set_id];
      int cg1 = offset.Super1( );
      int cg2 = offset.Super2( );
      int len1 = offset.Slen1( ); 
      int rc2 = offset.Rc2( );
      
      vec< vec<String> > table;
      vec<String> line;
      for (size_t cluster_id=0; cluster_id<offset.NClusters( ); cluster_id++) {
	table.clear( );
	const vec<SLink> &links = offset.Links( cluster_id );

	for (int link_id=0; link_id<links.isize( ); link_id++) {
	  line.clear( );
	  
	  int lib_id = pairs.libraryID( links[link_id].pair_id_ );
	  int gap = links[link_id].offset_.mu_ - len1;
	  
	  line.push_back( ToString( cg1 ) );
	  line.push_back( ToString( cg2 ) );
	  line.push_back( rc2 ? "-" : "+" );
	  line.push_back( ToString( cluster_id ) );
	  line.push_back( ToString( lib_id ) );
	  line.push_back( ToString( gap ) );
	  
	  table.push_back( line );
	}
	PrintTabular( cout, table, 2, "rrrrrr" );
	cout << "\n";
      }
    }
  }

  // Done.
  cout << "EOF" << endl;
  
}
