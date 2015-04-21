/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: AlignHyperToReference

   Take an assembly (a set of connected components of one HyperKmerPath), and
   for each <component>, find stretches of reference <captured> by that component.
   
   (original comment:)
   Apply FilterByReference() to the aligns for a given HyperKmerPath.  
   Output the "unwindings" of the graph components onto the finished sequence.

  @file
*/

#include "MainTools.h"

#include "Alignment.h"
#include "paths/KmerBaseBroker.h"
#include "paths/AlignHyperKmerPath.h"
#include "feudal/BinaryStream.h"

int main( int argc, char** argv )
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_Doc( HYPERKMERPATH, "complete pathname to the file containing the assembly as one HyperKmerPath with multiple components"  );

  // Either ALIGNS is specified...
  CommandArgument_String_OrDefault( ALIGNS, "" );

  // or TMPDIR is specified.
  CommandArgument_String_OrDefault( TMPDIR, "" ); // relative to cwd

  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String_OrDefault( GENOME, "genome" );
  CommandArgument_String_Doc( WORKDIR, "complete pathname to the run directory" );

  CommandArgument_Bool_OrDefault( DUMP_PATHS, False );

  CommandArgument_UnsignedInt_OrDefault( MIN_LENGTH, 1000 );
  CommandArgument_UnsignedInt_OrDefault( MIN_PERCENT, 10 );

  CommandArgument_Bool_OrDefault( REORDER, True );

  CommandArgument_Bool_OrDefault( ALIGN_DOM_FILTER, True );
  CommandArgument_Bool_OrDefault( EDGE_DOM_FILTER, True );
  CommandArgument_Bool_OrDefault( VERTEX_DOM_FILTER, False );

  CommandArgument_Bool_OrDefault( BRIEF, True );
  CommandArgument_Bool_OrDefault( DIPLOID, False );

  EndCommandArguments;

  if ( EDGE_DOM_FILTER && VERTEX_DOM_FILTER ) {
    cout << endl;
    cout << "Only one of EDGE_DOM_FILTER and VERTEX_DOM_FILTER may be True." << endl;

    TracebackThisProcess();
  }

  cout << Date() << ": loading graph" << endl;

  HyperKmerPath theGraph;
  BinaryReader::readFile( HYPERKMERPATH, &theGraph );
  theGraph.TestValid();

  unsigned int K = theGraph.K();

  ofstream out( "hkp.dot" );
  theGraph.PrintSummaryDOT0w( out, True, True, True );
  out.close();

  vec<look_align> aligns;
  vec< vec<int> > alignsIndex;

  KmerBaseBroker kbb( WORKDIR, K );

  if ( ALIGNS.empty() || ! IsRegularFile( ALIGNS ) ) {
    cout << Date() << ": generating aligns" << endl;

    String dataDir = PRE + "/" + DATA;
    AlignHyperKmerPath( theGraph, &kbb, dataDir + "/" + GENOME,
                        TMPDIR, aligns, alignsIndex );

    if ( ! ALIGNS.empty() )
      WriteLookAlignBinary( ALIGNS, aligns );
  }
  else {
    cout << Date() << ": loading aligns" << endl;

    LoadLookAlignBinary( ALIGNS, aligns );
  
    int numEdges = theGraph.Edges().size();
    alignsIndex.resize( numEdges );
    for ( int i = 0; i < aligns.isize(); ++i )
      alignsIndex[ aligns[i].query_id ].push_back( i );
  }

  cout << Date() << ": calling FilterByReference()" << endl;

  vec<TrustedPath> trustedPaths;
  FilterByReference( theGraph, K, aligns, alignsIndex, trustedPaths );
  cout << " got " << trustedPaths.size() << " trusted paths." << endl;

  if ( MIN_LENGTH > 0 ) {
    cout << Date() << ": calling FilterPathsByLength() on " << trustedPaths.size() << " trusted paths." << endl;
    FilterPathsByLength( trustedPaths, MIN_LENGTH, MIN_PERCENT );
  }

  if ( ALIGN_DOM_FILTER ) {
    cout << Date() << ": calling FilterPathsByAlignDominance() on " << trustedPaths.size() << " trusted paths." << endl;
    FilterPathsByAlignDominance( trustedPaths );
  }

  if ( EDGE_DOM_FILTER ) {
    cout << Date() << ": calling FilterPathsByEdgeDominance() on " << trustedPaths.size() << " trusted paths." << endl;
    FilterPathsByEdgeDominance( trustedPaths, theGraph.Edges().size() );
  }

  if ( VERTEX_DOM_FILTER ) {
    cout << Date() << ": calling FilterPathsByVertexDominance() on " << trustedPaths.size() << " trusted paths." << endl;
    FilterPathsByVertexDominance( trustedPaths );
  }

  if ( DUMP_PATHS ) {
    cout << Date() << ": dumping paths for " << trustedPaths.size() << " trusted paths." << endl;
    copy( trustedPaths.begin(), trustedPaths.end(),
          ostream_iterator<TrustedPath>( cout, "\n" ) );
  }

  cout << Date() << ": calling TrustedPathsToIndexedAligns() on " << trustedPaths.size() << " trusted paths." << endl;
  TrustedPathsToIndexedAligns( trustedPaths, theGraph.EdgeObjectCount( ),
                               aligns, alignsIndex );
  if ( REORDER ) {
    cout << Date() << ": calling ReorderToFollowReference()" << endl;
    ReorderToFollowReference( theGraph, aligns, alignsIndex );
  }

  vecbasevector genome( PRE + "/" + DATA + "/" + GENOME + ".fastb" );
  PrintAlignedHyperKmerPath( cout, theGraph, &kbb, genome, aligns, alignsIndex,
                             True, &trustedPaths, BRIEF, DIPLOID );    

  const int MIN_LEN = 3000;
  const int MIN_LEN_NO_REPORT = 1000;
  const int MIN_SEP = 2000;
  cout << "\n";
  ReportMisassemblies( cout, theGraph, aligns, alignsIndex, 
       MIN_LEN, MIN_SEP, MIN_LEN_NO_REPORT );
  cout << "\n";

}
