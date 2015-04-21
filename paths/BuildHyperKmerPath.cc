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
#include "paths/HyperFastavector.h"
#include "paths/HyperKmerPath.h"
#include "paths/Unipath.h"
#include "util/RunCommand.h"

/**
 * BuildHyperKmerPath
 *
 * Load various *.paths and *.unipaths files and generate the
 * corresponding HyperKmerPath.
 *
 * Input:
 *   <HEAD_IN>.paths{,db,_rc}.k<K>
 *   <HEAD_IN>.unipaths{,db}.k<K>
 *   <HEAD_IN>.unibases.k<K> (optional)
 *
 * Output:
 *   <HEAD_OUT>.HyperKemrPath.k<K>
 *   <HEAD_OUT>.HyperFastavector.k<K> (optional)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( HEAD_IN );
  CommandArgument_String_OrDefault( HEAD_OUT, HEAD_IN );
  CommandArgument_Bool_OrDefault( HYPER_FASTAVECTOR, True );
  EndCommandArguments;

  // File names.
  String strK = "k" + ToString( K );

  String paths_file      = HEAD_IN + ".paths." + strK;
  String pathsrc_file    = HEAD_IN + ".paths_rc." + strK;
  String pathsdb_file    = HEAD_IN + ".pathsdb." + strK;
  String unipaths_file   = HEAD_IN + ".unipaths." + strK;
  String unipathsdb_file = HEAD_IN + ".unipathsdb." + strK;
  String unibases_file   = HEAD_IN + ".unibases." + strK;
  
  String hyperKmerPath_file     = HEAD_OUT + ".HyperKmerPath." + strK;
  String hyperFastavector_file  = HEAD_OUT + ".HyperFastavector." + strK;
  
  vec<String> needed;
  needed.push_back( paths_file );
  needed.push_back( pathsrc_file );
  needed.push_back( pathsdb_file );
  needed.push_back( unipaths_file );
  needed.push_back( unipathsdb_file );
  if ( HYPER_FASTAVECTOR ) needed.push_back( unibases_file );

  // File(s) not found: fatal error.
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  // Load.
  cout << Date( ) << ": loading files" << endl;
  vecKmerPath paths( paths_file );
  vecKmerPath pathsrc( pathsrc_file );
  vecKmerPath unipaths( unipaths_file );
  BREAD2( pathsdb_file, vec<tagged_rpint>, pathsdb );
  BREAD2( unipathsdb_file, vec<tagged_rpint>, unipathsdb );

  // Generate adjacency graph and then HyperKmerPath
  cout << Date( ) << ": generating HyperKmerPath" << endl;
  digraph theDigraph;
  BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb,
			      unipaths, unipathsdb, theDigraph );

  HyperKmerPath theHKP;
  BuildUnipathAdjacencyHyperKmerPath( K, theDigraph, unipaths, theHKP );
  
  // Save HyperKmerPath.
  cout << Date( ) << ": saving HyperKmerPath" << endl;
  BinaryWriter::writeFile( hyperKmerPath_file, theHKP );

  // Build HyperFastavector.
  if ( HYPER_FASTAVECTOR ) {
    cout << Date( ) << ": generating HyperFastavector" << endl;
    vec<fastavector> edges;
    {
      vecbvec bases( unibases_file );
      edges.reserve( bases.size( ) );
      for (size_t ii=0; ii<bases.size( ); ii++)
	edges.push_back( fastavector( bases[ii] ) );
    }
    HyperFastavector theHFV( K, theHKP, edges );
    
    cout << Date( ) << ": saving  HyperFastavector" << endl;
    BinaryWriter::writeFile( hyperFastavector_file, theHFV );
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

