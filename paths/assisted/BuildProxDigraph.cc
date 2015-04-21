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
#include "lookup/LookAlign.h"
#include "random/Shuffle.h"
#include "paths/Alignlet.h"
#include "paths/assisted/LibLinksUtils.h"
#include "paths/assisted/Proxgraph.h"
#include "util/RunCommand.h"

#include "paths/Sepdev.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/MakeScaffoldsScoredBest.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/PairDistCorrection.h"

/**
 * BuildProxDigraph
 *
 * Generate a digraphE of CProx edges. Vertices are oriented contigs
 * 0[+], 0[-], 1[+], 1[-], ...
 *
 * RUN_DIR: where jump and long_jump reads are
 * JUMPS: head of jump reads
 * LONG_JUMPS: head of long jump reads
 * USE_REF: generate proximity by using aligns on ref
 * DISCARD_REPEATS: discard contigs aligning multiply on ref
 * MAX_CN: discard contigs with high copy number
 * MIN_GAP: min allowed gap size
 * MAX_GAP: max allowed gap size
 * MIN_LINKS: min number of links between vertices
 * VERBOSITY: from 0 (almost no logging), to 1 (max logging)
 */


/*
 * Translate unibase align triplets to alignlets for use in scaffolding code.
 * We randomize aligments to randomize placements of multiple aligns for given read.
 */

void UnibaseAlignsToContigAlignlets(const vecbvec &reads,
				    const vecbvec &unibases,
				    const vec<int> &origIds,
				    const vec< triple<int64_t,int64_t,int> > & aligns,
				    const int randomSeed,
				    vec<alignlet> &alignlets,
				    vec<int> & alignletIndex)
{
  // We only care about a subset of unibases corresponding to contigs;
  // flag them and create a reverse map back to contig id.
  size_t nUni = unibases.size();
  vec<bool> uniKeepers(nUni, false);
  vec<int> uniMap(nUni, -1);
  for (size_t i = 0; i < origIds.size(); ++i) {
    int u = origIds[i];
    uniKeepers[u] = true;
    uniMap[u] = i;
  }
  
  alignlets.reserve(aligns.size());
  alignletIndex.clear();
  alignletIndex.resize(reads.size(), -2);
  vec<int> randomize;
  Shuffle(aligns.size(), randomize, randomSeed);

  for (size_t i = 0; i < aligns.size(); i++) {
    int ai = randomize[i];
    int r = aligns[ai].first, u = aligns[ai].second, p = aligns[ai].third;
    // Filter by unibases of interest
    if (!uniKeepers[u]) continue;
    int contig = uniMap[u];
    int pos2 = ( p >= 0 ? p : -p - 1 );
    alignletIndex[r] = alignlets.size();
    alignlets.push(pos2, pos2 + reads[r].isize(),
		   contig, unibases[u].size(), p >= 0);
  }
}


int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( CONTIGS_HEAD );
  CommandArgument_String( HEAD_UNI );
  CommandArgument_String( RUN_DIR );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_String_OrDefault( JUMPS, "" );
  CommandArgument_String_OrDefault( LONG_JUMPS, "" );
  CommandArgument_Bool_OrDefault( USE_REF, True );
  CommandArgument_Bool_OrDefault( DISCARD_REPEATS, False );
  CommandArgument_Bool_OrDefault( REGAP, True );
  CommandArgument_Int_OrDefault( MAX_CN, 999999 );
  CommandArgument_Int_OrDefault( MIN_GAP, -500 );
  CommandArgument_Int_OrDefault( MAX_GAP, 2500 );
  CommandArgument_Int_OrDefault( MIN_LINKS, 4 );
  CommandArgument_Int_OrDefault( SEED, 42 );
  CommandArgument_Int_OrDefault( VERBOSITY, 0 );
  CommandArgument_Bool_OrDefault( USE_CLB, True );
  EndCommandArguments;

  String strK = ToString( K );

  // Debug output loggging
  ostream *log = VERBOSITY > 0 ? (&cout) : NULL;

  // Dir and file names.
  String contigsF = CONTIGS_HEAD + ".fastb";
  String contigOrigIdsF = CONTIGS_HEAD + ".orig.ids";

  String digraphF = CONTIGS_HEAD + ".CProx.digraphE";

  String uniOnRefAlignsF = CONTIGS_HEAD + ".on_reference.qlt";
  String uniCnsF = CONTIGS_HEAD + ".cns";
  
  String jumpHead = RUN_DIR + "/" + JUMPS;
  String jumpPairsF = jumpHead + ".pairs";
  String jumpAlignsF = jumpHead + ".aligns";
  String jumpReadsF = jumpHead + ".fastb";
  
  String unibasesF = HEAD_UNI + ".unibases.k" + strK;

  bool useJumps = ( JUMPS != "" );

  if ( ! ( USE_REF || useJumps ) ) {
    cout << "Fatal error. It is not possible to build a graph\n"
	 << "without jump reads and without a reference.\n" << endl;
    return 1;
  }
  

  vec<String> needed;
  needed.push_back(contigsF);
  needed.push_back(unibasesF);
  needed.push_back(uniCnsF);
  if ( USE_REF ) needed.push_back(uniOnRefAlignsF);
  if ( useJumps ) {
    needed.push_back( jumpPairsF );
    needed.push_back( jumpReadsF );
    needed.push_back( jumpAlignsF );
  }
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Raw ref-based and link-based info.
  CLibLinks raw_onref;
  vec<CLibLinks> raw_jumps;
  vec<CLibLinks> raw_Jumps;
  
  // Load contigs and generate trivial superb
  cout << Date() << ": loading contigs" << flush;
  BaseVecVec contigs;
  contigs.ReadAll(contigsF);
  size_t ntigs = contigs.size();
  cout << " (" << contigs.size() << ")" << endl;

  cout << Date() << ": loading unibases" << flush;
  BaseVecVec unibases;
  unibases.ReadAll(unibasesF);
  cout << " (" << unibases.size() << ")" << endl;

  // Load CN info.
  cout << Date( ) << ": loading contig->unibase map" << flush;
  vec<int> contigOrigIds;
  READX( contigOrigIdsF, contigOrigIds );
  ForceAssertEq( contigOrigIds.size( ), contigs.size( ) );

  /* only used if we plug in old scaffolding code
  vec<fastavector> fasta_contigs;
  for (size_t i = 0; i < ntigs; ++i) {
    fastavector fasta(contigs[i]);
    fasta_contigs.push_back(fasta);
  }
  */

  // Generate initial supers.
  vec<superb> supers(contigs.size());
  for (int i=0; i<(int)supers.size( ); i++) {
    supers[i].PlaceFirstTig(i, contigs[i].size());
  }

  // Add linking info from jumps and long jumps.
  PairsManager jpairs;
  PairsManager Jpairs;

  // Is repetitive tag (only used in some cases).
  vec<bool> isrep( ntigs, false );
  
  // Load CN info.
  cout << Date( ) << ": loading cns" << endl;
  vec<int> cns;
  READX( uniCnsF, cns );
  ForceAssertEq( cns.size( ), isrep.size( ) );

  // Add proximity info from aligns of contigs on reference.
  if ( USE_REF ) {
    cout << Date( ) << ": loading aligns of contigs on reference" << flush;
    vec<look_align> hits;
    LoadLookAligns( uniOnRefAlignsF, hits );
    cout << " (" << hits.size() << ")" << endl;

    if ( DISCARD_REPEATS ) {
      cout << Date( ) << ": tagging repetitive contigs... " << flush;
      vec<look_align> unique_hits;
      unique_hits.reserve( hits.size( ) );
      vec<int> nhits( ntigs, 0 );
      for (size_t ii=0; ii<hits.size( ); ii++)
	nhits[ hits[ii].query_id ] += 1;
      for (size_t ii=0; ii<hits.size( ); ii++)
	if ( nhits[ii] < 2 ) unique_hits.push_back( hits[ii] );
	else isrep[ hits[ii].query_id ] = true;
      int n_rep = 0;
      for (size_t ii=0; ii<isrep.size( ); ii++) {
	//if (cns[ii] > MAX_CN) isrep[ii] = true;
	if ( isrep[ii] ) n_rep++;
      }
      cout << n_rep << " tagged" << endl;
      swap( hits, unique_hits );
    }

    order_lookalign_TargetBeginEnd sorter;
    if ( ! is_sorted( hits.begin( ), hits.end( ) ) )
      sort( hits.begin( ), hits.end( ), sorter );
    
    cout << Date( ) << ": collecting separation on reference info" << endl;
    FillRefLinks( MIN_GAP, MAX_GAP, MAX_CN, hits, cns, raw_onref );
  }

  digraphE<CLinkBundle> clbGraph;

  if (useJumps) {
    cout << Date( ) << ": loading jump pairs" << endl;
    jpairs.Read(jumpPairsF);
    
    cout << Date( ) << ": loading jump reads" << flush;
    vecbvec jumpReads(jumpReadsF);
    cout << " (" << jumpReads.size() << ")" << endl;

    cout << Date( ) << ": loading jump aligns" << flush;
    vec< triple<int64_t,int64_t,int> > multiAligns;
    BinaryReader::readFile(jumpAlignsF, &multiAligns );
    cout << " (" << multiAligns.size() << ")" << endl;

    cout << Date() << ": translating aligns to alignlets" << flush;
    vec<alignlet> aligns;
    vec<int> idx;
    UnibaseAlignsToContigAlignlets(jumpReads, unibases, contigOrigIds,
				   multiAligns, SEED, aligns, idx);
    cout << " (" << aligns.size() << ")" << endl;
    

    if (REGAP) {
      cout << Date( ) << ": pair distance correction" << endl;
      // pair length correction 
      // correction
      vec<int> seps_empty;
      vec<int> seps;
      vec<int> sds;
      PairDistCorrectionFromIntDist(jumpHead, jpairs, aligns, idx, seps, sds, VERBOSITY>0);
      jpairs.AddSeps(seps);
    }

    if (USE_CLB) {
      digraphE<sepdev> unused;
      cout << Date( ) << ": building scaffold graph" << endl;
      BuildScaffoldGraph( jpairs, supers, aligns, idx, unused, &clbGraph);
    } else {
      cout << Date( ) << ": collecting linking info" << endl;
      FillLibLinks( MAX_CN, jpairs, isrep, cns, aligns, idx, raw_jumps );
    }
  }
  
  // Dump raw links.
  if ( VERBOSITY > 4 ) {
    cout << "\n";
    if ( useJumps )
      for (size_t ii=0; ii<raw_jumps.size( ); ii++)
	raw_jumps[ii].PrintLinks( cout );
    if ( useJumps )
      for (size_t ii=0; ii<raw_jumps.size( ); ii++)
	raw_Jumps[ii].PrintLinks( cout );
    if ( USE_REF )
      raw_onref.PrintLinks( cout );
  }
  
  // Generate a CProxInfo metadata container.
  CProxInfo info( MIN_GAP, MAX_GAP, MIN_LINKS );
  info.SetLibTags( ( useJumps ? &jpairs : 0) , ( useJumps ? &Jpairs : 0 ) );

  // Combine raw info into a set of 
  cout << Date( ) << ": building digraphE" << endl;
  vec<CLibLinks> *jraw = useJumps ? &raw_jumps : 0;
  vec<CLibLinks> *Jraw = useJumps ? &raw_Jumps : 0;
  CLibLinks *refraw = USE_REF ? &raw_onref : 0;
  proxgraph graphE;
  LibLinksToDigraphE( ntigs, jraw, Jraw, refraw, &info, graphE );

  if (USE_CLB) {
    cout << Date( ) << ": transfering link bundles to CProx graph" << endl;
    AddBundlesToDigraphE( clbGraph, &info, graphE );
  }
  
  // Save.
  cout << Date( ) << ": saving digraphE" << endl;
  BinaryWriter::writeFile( digraphF, graphE );
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

