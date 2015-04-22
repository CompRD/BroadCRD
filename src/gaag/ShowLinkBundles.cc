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
#include "PairsManager.h"
#include "paths/Alignlet.h"
#include "paths/Sepdev.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/PairDistCorrection.h"

/**
 * ShowLinkBundles
 *
 * READS: basename of reads file; input: READS.pairs
 * CONTIGS: basename of fastb file with READS alinged to them; input:
 *    CONTIGS.fastb
 *    CONTIGS.qltoutlet
 *    CONTIGS.qltoutlet.index
 * MIN_LINKS: min number of links between vertices
 * REGAP: use better pair separate hueristics
 * OUTPUT: output file for bundle information (shown 0-based).
 * --bruce 11 Dec 2012
 */

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( CONTIGS );
  CommandArgument_String( READS );
  CommandArgument_Int_OrDefault( MIN_LINKS, 2 );
  CommandArgument_Bool_OrDefault( REGAP, True );
  CommandArgument_String_OrDefault( OUTPUT, "bundles.txt" );
  EndCommandArguments;

  // Load contigs and generate trivial supers (one scaffold per contig)
  cout << Date() << ": loading contigs" << endl;
  BaseVecVec contigs;
  contigs.ReadAll(CONTIGS + ".fastb");
  size_t ntigs = contigs.size();

  vec<superb> supers(contigs.size());
  for (int i=0; i<supers.isize(); i++) {
    supers[i].PlaceFirstTig(i, contigs[i].size());
  }

  digraphE<CLinkBundle> bundle_graph;

  // Load read alignment info
  PairsManager pairs;
  cout << Date() << ": loading reads pairing info" << endl;
  pairs.Read(READS + ".pairs");
    
  vec<alignlet> aligns;
  cout << Date() << ": loading alignments" << endl;
  BinaryReader::readFile( CONTIGS + ".qltoutlet", &aligns );
    
  cout << Date() << ": loading alignment index" << endl;
  vec<int> idx;
  BinaryReader::readFile( CONTIGS + ".qltoutlet.index" , &idx );
    
  if (REGAP) {
    cout << Date( ) << ": pair distance correction" << endl;
    // pair length correction 
    vec<int> seps;
    vec<int> sds;
    PairDistCorrectionFromIntDist(READS, pairs, aligns, idx, seps, sds, False);
    pairs.AddSeps(seps);
  }

  digraphE<sepdev> unused;
  cout << Date( ) << ": building scaffold graph" << endl;
  BuildScaffoldGraph( pairs, supers, aligns, idx, unused, &bundle_graph, NULL, NULL);

  cout << Date( ) << ": writing bundle information to " << OUTPUT << endl;
  ofstream fout;
  fout.open(OUTPUT.c_str());
  fout << "# c1 c2 gap/stdev weight score c1window c2window c1size c2size" << endl;

  vec<int> to_left;
  vec<int> to_right;
  bundle_graph.ToLeft( to_left );
  bundle_graph.ToRight( to_right );
  for (int edge_id=0; edge_id < bundle_graph.EdgeObjectCount( ); ++edge_id) {
    CLinkBundle bundle = bundle_graph.EdgeObject(edge_id);
    if (bundle.Weight() < MIN_LINKS) continue;

    int v1 = to_left[edge_id];
    int c1 = v1 / 2;
    int rc1 = v1 & 1;
    int v2 = to_right[edge_id];
    int c2 = v2 / 2;
    int rc2 = v2 & 1;

    fout << c1 << (rc1 ? "-" : "+") << " " 
	 << c2 << (rc2 ? "-" : "+") << " " 
	 << bundle.AsString() << " "
	 << contigs[c1].size() << " "
	 << contigs[c2].size() << endl;

  }
  // Done.
  cout << Date( ) << ": done" << endl;
}

