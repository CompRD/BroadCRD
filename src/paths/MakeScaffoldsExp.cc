///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


// MakeScaffoldsExp.  Build scaffolds from assembly graph.  Run this
// after AlignPairsToHyperLG.  This is an experimental version of
// MakeScaffolds using a different core algorithm to evaluate joins.
// It's sufficiently different to keep on a separate code track for
// now, though it seems to generally produce better results. It will
// probably either replace MakeScaffolds or wither and die.  Code is
// still unclean with a few sections commented out and/or alternative
// versions of things to experiment with. I will refactor, clean,
// and/or abandon over time. --bruce 16oct09

/* MakeScaffolds.  Build scaffolds from assembly graph.
 * Run this after AlignPairsToFasta.
 *
 * Input:
 * -- Reads and alignments in sub_dir
 * -- Contigs at sub_dir/<SCAFFOLDS_IN>.contigs.fasta
 * -- Scaffolds at sub_dir/<SCAFFOLDS_IN>.superb (may be trivial/degenerate)
 *
 * Output:
 * -- Contigs at sub_dir/<SCAFFOLDS_OUT>.contigs.fasta
 * -- Assembly at sub_dir/<SCAFFOLDS_OUT>.assembly.fasta
 * -- Scaffolds at sub_dir/<SCAFFOLDS_OUT>.superb
 *
 *
 *
 ******************************************************************************/

#include "Equiv.h"
#include "Fastavector.h"
#include "Charvector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "feudal/FeudalTools.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "paths/HyperFastavector.h"
#include "paths/PseudoUnipatherLib.h"
#include "paths/SuckScaffolds.h"
#include "paths/BreakBadScaffolds.h"
#include "paths/ContigsManager.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/Sepdev.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/CleanseScaffoldGraph.h"
#include "paths/ScaffoldGraphIntegrity.h"
#include "paths/MergeScaffoldGraph.h"
#include <sstream>

/**
 * Bscore: bruce's attempt to measure of bundle quality. don't take seriously
 */
double Bscore(const CLinkBundle &b)
{
  //return (score_ * score_ - 1.0) * sqrt((double) weight_);
  // actually seems to be distributed around 0.85 --bruce
  double delta = b.score_ - 0.85;
  double sepk = (b.sep_ < 1000) ? 1.0 : (b.sep_/1000.0);
  return delta * delta * b.dev_ / sepk;
}



// Clean graph (remove bad scores, bad gaps, etc)

void CleanupGraph(digraphE<CLinkBundle> &G, const int mingap, const int maxgap, const double maxscore, const int minlinks, const int VERBOSITY)
{
  vec<int> to_delete;

  for (int e = 0; e < G.EdgeObjectCount(); ++e) {
    bool del = false;
    CLinkBundle cl = G.EdgeObject(e);
    if (cl.Sep() < mingap || cl.Sep() > maxgap) del = true;
    else if (minlinks > 0 && cl.Weight() < minlinks) del = true;
    else if (Bscore(cl) > maxscore) del = true;
    if (del && VERBOSITY>=3) cout << "Delete edge: " << G.EdgeObject(e).AsString() << endl;
    if (del) to_delete.push_back(e);
  }
  G.DeleteEdges(to_delete);
}



// Strategies for taking a step forward in the link graph.  Return -1 or next graph node.

// Simple step means there is only one link forward, and we're the only link to the target.
int stepforward_simple(const int start,
		       const digraphE<CLinkBundle> &links,
		       const vec< triple<int,int,int> > &path,
		       const vec<Bool> &used)
{
  vec<int> const &from = links.From(start);
  // make sure there's one link foward
  if (from.size() != 1) return -1;
  int to = from[0];

  if (used[to/2]) return -1;
  
  vec<int> const &to2 = links.To(to);
  for (u_int i = 0; i < to2.size(); ++i) {
    bool in_path = false;
    if (to2[i] == start) continue;
    for (u_int j = 0; j < path.size(); ++j) {
      if (to2[i] == path[j].first) {
	in_path = true;
	break;
      }
    }
    if (!in_path) return -1;
  }
  return to;
}


// Find a link forward for which we are the only thing linking to it,
// and through which we can get to all our other forward links.

// looks forward recursively through all branches collecting nodes we
// can see up through a specified depth
void
peekforward(const int start, const digraphE<CLinkBundle> &links, const int depth, vec<int> &nodes) {
  vec<int> const &from = links.From(start);
  for (u_int i = 0; i < from.size(); ++i) {
    int n = from[i];
    nodes.push_back(n);
    if (depth > 0) peekforward(n, links, depth-1, nodes);
  }
}

#define PEEKFORWARD_DEPTH 4

int stepforward_graph(const int start, const digraphE<CLinkBundle> &links, vec<Bool> &used)
{
  vec<int> const &from = links.From(start);
  u_int nlinks = from.size();

  // this only makes sense to do if there are two are more links
  if (nlinks < 2) return -1;

  for (u_int i = 0; i < nlinks; ++i) {
    int next = from[i];
    if (used[next/2]) continue;

    // don't consider if there are other unused links coming into it
    vec<int> const &to_next = links.To(next);
    bool OK = true;
    for (u_int j = 0; j < to_next.size(); ++j) {
      if (to_next[j] != start && !used[to_next[j]/2]) {
	OK = false;
	break;
      }
    }
    if (!OK) continue;

    // collect the nodes we can see above the horizon
    vec<int> nodes;
    peekforward(next, links, PEEKFORWARD_DEPTH, nodes);

    // for each other link neighbor, can we get there through "next"?
    OK = true;
    for (u_int j = 0; j < nlinks; ++j) {
      if (j == i) continue;
      int neighbor = from[j];
      if (find(nodes.begin(), nodes.end(), neighbor) == nodes.end()) {
	  OK = false;
	  break;
      }
    }
    // if still OK, we've seen all neighbors, we we have a winner
    if (OK) return next;
  }
  return -1;
}

// Look for a dominant link forward; that is, one which has
// dramatically more read link support than the others.
#define DOMINATION 10.0		// link wins if it has 5x more reads supporting it than others

int stepforward_dominant(const int start, const digraphE<CLinkBundle> &links, vec<Bool> &used)
{
  vec<int> const &from = links.From(start);
  u_int nlinks = from.size();

  // this only makes sense to do if there are two are more links
  if (nlinks < 2) return -1;

  int max_strength = 0, strongest = -1;
  vec<int> strengths(nlinks);

  // Record strengths & strongest link
  for (u_int i = 0; i < nlinks; ++i) {
    int to = from[i];
    CLinkBundle &s = links.EdgeObjectsBetween(start, to)[0];
    strengths[i] = s.Weight();
    if (s.Weight() > max_strength) {
      max_strength = s.Weight();
      strongest = i;
    }
  }

  int next = from[strongest];
  if (used[next/2]) return -1;

  // Compute ratio of strongest to next strongest
  float dominance = 1e6;
  for (u_int i = 0; i < nlinks; ++i) {
    if ((int)i == strongest) continue;
    float d = (float)max_strength / (float)strengths[i];
    dominance = min(dominance, d);
  }

  // if sufficiently dominant & not already used, it wins
  if (dominance >= DOMINATION) return next;

  return -1;
}

#define DEV_MULT  3

int stepforward_nearest(const int start, const digraphE<CLinkBundle> &links,
			const vec<superb> &scaffolds, vec<Bool> &used)
{
  vec<int> const &from = links.From(start);
  u_int nlinks = from.size();

  // this only makes sense to do if there are two are more links
  if (nlinks < 2) return -1;

  int nearest = -1, min_sep = 1000000000;
  int sep2 = 0, dev2 = 0;

#ifdef notdef
  // Record nearest neighbor
  for (u_int i = 0; i < nlinks; ++i) {
    int to = from[i];
    CLinkBundle b = links.EdgeObjectsBetween(start, to)[0];
    if (b.Sep() < min_sep) {
      min_sep = b.Sep();
      nearest = i;
    }
  }

  int s2 = from[nearest]/2;
  if (used[s2]) return -1;

  // Now check all other links to see if any could be closer
  sep2 = min_sep;
  const superb &S2 = scaffolds.at(s2);
  int s2_len = S2.SubSuperLength(0, S2.Ntigs() - 1);
  int s2_dev = S2.SubSuperLengthDev(0, S2.Ntigs() - 1);
  //PRINT4(s2, sep2, s2_len, s2_dev);
  for (u_int i = 0; i < nlinks; ++i) {
    if ((int)i == nearest) continue;
    int to = from[i];
    CLinkBundle b = links.EdgeObjectsBetween(start, to)[0];
    int s2x = to/2;
    int sep2x = b.Sep();
    int dev2x = b.Dev();
    const superb &S2x = scaffolds.at(s2x);
    int s2x_len = S2x.SubSuperLength(0, S2x.Ntigs() - 1);
    int s2x_dev = S2x.SubSuperLengthDev(0, S2x.Ntigs() - 1);
    //PRINT4(s2x, sep2x, s2x_len, s2x_dev);
    if (sep2x + s2x_len - DEV_MULT * (s2x_dev + dev2x) <= sep2 + s2_len)
      return -1;
  }
#else
  int min_dist = 1000000000;
  // Record nearest neighbor, using distance to far end
  for (u_int i = 0; i < nlinks; ++i) {
    int to = from[i];
    int s2 = to/2;
    const superb &S2 = scaffolds[s2];
    int s2_len = S2.SubSuperLength(0, S2.Ntigs() - 1);
    CLinkBundle b = links.EdgeObjectsBetween(start, to)[0];
    int dist = s2_len + b.Sep();
    
    //PRINT4(s2,s2_len,b.Sep(), dist);
    if (dist < min_dist) {
      min_dist = dist;
      sep2 = b.Sep();
      dev2 = b.Dev();
      nearest = i;
    }
  }

  int s2 = from[nearest]/2;
  if (used[s2]) return -1;

  for (u_int i = 0; i < nlinks; ++i) {
    if ((int)i == nearest) continue;
    int to = from[i];
    CLinkBundle b = links.EdgeObjectsBetween(start, to)[0];
    int s2x = to/2;
    int sep2x = b.Sep();
    int dev2x = b.Dev();
    // if our candidate coudn't reasonably be closest on the near end, punt
    //PRINT4(sep2,dev2,sep2x,dev2x);
    if (sep2x+DEV_MULT*dev2x/2 < sep2-DEV_MULT*dev2/2)
      return -1;
  }
#endif

  return from[nearest];
}


int simple = 0;

vec< triple<int,int,int> >
walkforward(const int start, const digraphE<CLinkBundle> &links, const vec<superb> &scaffolds, vec<Bool> &used) {
  vec< triple<int,int,int> > path;
  if (used[start/2]) return path;
  int v1 = start;
  while (v1 >= 0) {
    u_int n = links.From(v1).size();
    if (n == 0) break;

    int v2 = -1;

    // try several strategies
    // first, is there only one possibility forward
    v2 = stepforward_simple(v1, links, path, used);
    if (v2 >= 0) {
      ++simple;
    } else if (n > 1) {
      int v2g, v2d, v2n;
      
      // if there is more than one possible forward link, try these
      // methods to see if there is a unique answer or a majority
      v2g = stepforward_graph(v1, links, used);
      v2n = stepforward_nearest(v1, links, scaffolds, used);
      v2d = stepforward_dominant(v1, links, used);
      if (v2g >= 0 && ((v2g == v2n || v2g == v2d) || (v2n < 0 && v2d < 0)))
	v2 = v2g;
      else if (v2n >= 0 && ((v2n == v2d) || (v2g < 0 && v2d < 0)))
	v2 = v2n;
      else if (v2d >= 0 && (v2g < 0 && v2n < 0))
	v2 = v2d;
      //PRINT4(v1,v2g,v2n,v2d);
    }
    if (v2 >= 0) {
      CLinkBundle s = links.EdgeObjectsBetween(v1, v2)[0];
      if (path.size() == 0) path.push(start, 0, 0);
      path.push(v2, s.Sep(), s.Dev());
      used[v1/2] = used[v2/2] = True;
    }
    v1 = v2;
  }
  return path;
}

vec< triple<int,int,int> >
ScaffoldComponent(const digraphE<CLinkBundle> &links, const vec<int> &component, vec<Bool> &used)
{
  vec< triple<int,int,int> > path;
  if (component.size() < 2) return path;
  PRINT2(component.size(), component[0]);
  for (u_int c = 1; c < component.size(); ++c) {
    if (component[c] == (component[c-1] ^ 1)) {
      PRINT2(component[c], component[c-1]);
      return path;
    }
  }
  vec<int> sources;
  links.SubgraphSources(component, sources);
  if (sources.size() > 1) PRINT(sources.size());
  // shave bogus sources
  for (u_int i = 0; i < sources.size(); ++i) {
    int s = sources[i];
    vec<int> nexts = links.From(s);
  }
  return path;
}


int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String(SUBDIR);

  // Input from <PRE>/<DATA>/<RUN>:
  //   <READS>.{fastb,pairs}

  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads" );

  // Input from <PRE>/<DATA>/<RUN>/ASSEMBLIES/<SUBDIR>:
  //   <ALIGNS>.qltoutlet{,.index},
  //   <SCAFFOLDS_IN>.{contigs.fasta,superb}
     
  CommandArgument_String_OrDefault( SCAFFOLDS_IN, "initial_scaffolds" );

  // Output saved in <PRE>/<DATA>/<RUN>/<SUBDIR>/<ASSEMBLIES>:
  //   <SCAFFOLDS_OUT>.{contigs.fasta,assembly.fasta,superb}
     
  CommandArgument_String_OrDefault( SCAFFOLDS_OUT, "scaffolds" );

  // If true, run SuckScaffold to absorb small scaffolds into the
  // gaps of larger ones.
  CommandArgument_Bool_OrDefault( SUCK_SCAFFOLDS, True );

  // These are passed to BuildScaffoldGraph
  CommandArgument_Double_OrDefault( MAX_SCORE, 1.65 );
  CommandArgument_Int_OrDefault( MIN_LINKS, 4 );
  CommandArgument_Double_OrDefault( MIN_SCORE, 0.5);
  CommandArgument_Int_OrDefault( RATIO_MIN_LINKS, 2 );
  CommandArgument_Double_OrDefault( RATIO_TO_SECOND, 3.0);
  CommandArgument_Int_OrDefault( MAX_OVERLAP, 350 );
  CommandArgument_Int_OrDefault( LOW_WEIGHT, 6 );
  CommandArgument_Bool_OrDefault(MERGE, False);
  CommandArgument_Bool_OrDefault(CLEANSE, True);
  CommandArgument_Bool_OrDefault(CLEANSE_ISOLATE_ONLY, False);
  CommandArgument_Int_OrDefault(VERBOSITY, 0);
  CommandArgument_Bool_OrDefault(WRITE, True);
  CommandArgument_Int_OrDefault( MIN_ALLOWED_GAP, -800 );
  CommandArgument_Int_OrDefault( MAX_ALLOWED_GAP, 1000000 );
  CommandArgument_Double_OrDefault( MAX_BSCORE, 1000000 );
  CommandArgument_Int_OrDefault(MAX_ITER, 1000000);
  CommandArgument_Bool_OrDefault(JOIN, True);
  CommandArgument_String_OrDefault( SCAFFOLD_GRAPH_OUT, "");
  CommandArgument_Bool_OrDefault( DOT, "False");
  CommandArgument_Bool_OrDefault( DOT_PER_ITER, "False");
  CommandArgument_String_OrDefault( DEBUG_CONTIGS, "{}" );

  EndCommandArguments;

  // for future use
  vec<int> debug_contigs;
  ParseIntSet(DEBUG_CONTIGS, debug_contigs, False);

  // Define heuristic constants.
  static const int min_allowed_gap = MIN_ALLOWED_GAP;
  static const int max_allowed_gap = MAX_ALLOWED_GAP;

  // Dir and file names.

  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String reads_file = run_dir + "/" + READS + ".fastb";
  String pairs_file = run_dir + "/" + READS + ".pairs";
  String contigs_in_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
  String scaffolds_in_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file = sub_dir + "/" + ALIGNS + ".qltoutlet.index";
     
  String out_suffix = SCAFFOLDS_OUT + ( JOIN ? "" : ".unjoined" );
  String contigs_out_file = sub_dir + "/" + out_suffix + ".contigs.fasta";
  String assembly_out_file = sub_dir + "/" + out_suffix + ".assembly.fasta";
  String scaffolds_out_file = sub_dir + "/" + out_suffix + ".superb";
     
  cout << Date( ) << ": loading scaffolds" << endl;
  vec<superb> scaffolds;
  ReadSuperbs( scaffolds_in_file, scaffolds );
  vec<fastavector> scaffold_contigs;
  LoadFromFastaFile( contigs_in_file, scaffold_contigs );
  size_t ntigs = scaffold_contigs.size();
     
  // Check that the files are consistent.
  size_t n_tigs_in_scaffolds = 0;
  for ( size_t i = 0; i < scaffolds.size(); i++ ) {
    n_tigs_in_scaffolds += scaffolds[i].Ntigs();
    for ( int j = 0; j < scaffolds[i].Ntigs(); j++ ) {
      int tig = scaffolds[i].Tig(j);
      ForceAssertEq( (int)scaffold_contigs[tig].size(), scaffolds[i].Len(j) );
    }
  }
  ForceAssertGe( ntigs, n_tigs_in_scaffolds );
     
  // Load alignments and convert to lightweight form.
  // aligns0_index works as this: say aligns0_index[id] = ii, then:
  //  ii = -2 :   no look_aligns for this read
  //  ii = -1 :   more than one look_align for this read (id)
  //  ii >= 0 :   aligns0[ii] is the only look_align for this read

  cout << Date( ) << ": loading aligns" << flush;
  vec<alignlet> aligns0;
  BinaryReader::readFile( aligns_file, &aligns0 );
  vec<int> aligns0_index;
  BinaryReader::readFile( index_file, &aligns0_index );
  longlong nreads = MastervecFileObjectCount( reads_file );
  cout << " (" << aligns0.size( ) << " aligns found)" << endl;

  // Load read pairs.

  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );

  vec<alignlet> dummy(0);
  vec< vec<int> > dummy_index(0);
  ContigsManager contigs_info( scaffold_contigs, aligns0, aligns0_index, dummy, dummy_index );

  // Track orientations of contigs.

  vec<Bool> rctig;

  cout << Date( ) << ": " << scaffolds.isize() << " initial scaffolds, " << ntigs << " initial contigs" << endl;

  // Now begin the scaffold joining operation.  This is repeated until there
  // no more joins are made.

  if ( VERBOSITY >= 5 )
    {    cout << "\ninitial scaffolds:\n\n";
      for ( int i = 0; i < scaffolds.isize( ); i++ )
	scaffolds[i].Print( cout, "scaffold_" + ToString(i) );    }

  int iter = 0;

  int joins;

  do {
    joins = 0;

    // Generate scaffold link graph.
	  
    cout << Date( ) << ": generating scaffold link graph" << endl;

    digraphE<sepdev> G;
    digraphE<CLinkBundle> Gbundles;
    ostream *log = ( VERBOSITY >= 5 ) ? &cout : 0;
    
    if (iter == 0 && MERGE) {
      // Merge overlapping edges (log within).
      MergeScaffoldGraph( scaffolds, scaffold_contigs, aligns0, aligns0_index, cout, pairs, VERBOSITY>0 );
      ForceAssert( ScaffoldGraphIntegrity( Gbundles ) );
    }
    BuildScaffoldGraph( pairs, scaffolds, aligns0,
			aligns0_index, G, &Gbundles, 0, log,
			MAX_SCORE, MIN_LINKS, MIN_SCORE, RATIO_MIN_LINKS,
			RATIO_TO_SECOND, MAX_OVERLAP, LOW_WEIGHT);
    if (CLEANSE) {
      // Note that MIN_LINKS and MAX_OVERLAP to Cleanse are used differrently
      // by CleanseScaffoldGraph, and should probably be separate parameters.
      // Must ISOLATE_ONLY on 2nd and subsequent iterations
      int cuts = CleanseScaffoldGraph( Gbundles,
				       scaffolds,
				       contigs_info,
				       cout,
				       VERBOSITY>0, 
				       CLEANSE_ISOLATE_ONLY || (iter>0),
				       MIN_LINKS,
				       MAX_OVERLAP);
      if (cuts)
	cout << Date( ) << ": Cleansing resulted in " << cuts << " cuts" << endl;
      ForceAssert( ScaffoldGraphIntegrity( Gbundles ) );
    }

    if ( SCAFFOLD_GRAPH_OUT != "" ) {
      String out_file = SCAFFOLD_GRAPH_OUT + ".before." + ToString( iter );
      BinaryWriter::writeFile( out_file, Gbundles );
    }


    int nscaff = scaffolds.size( );

    if ( SCAFFOLD_GRAPH_OUT != "" ) {
      String out_file = SCAFFOLD_GRAPH_OUT + ".after." + ToString( iter );
      BinaryWriter::writeFile( out_file, Gbundles );
    }

    // number of contigs can change as a result of CleanseScaffolds splitting some contigs;
    // this should only happen on the first iteration, or something is wrong.
    ntigs = scaffold_contigs.size();
    if (iter == 0) {
      rctig.resize(ntigs);
      for (u_int i = 0; i < ntigs; ++i) rctig[i] = false;
    } else ForceAssertEq(ntigs, rctig.size());

    cout << Date( ) << ": " << nscaff << " scaffolds, " << scaffold_contigs.isize() << " contigs" << endl;
        
    CleanupGraph( Gbundles, MIN_ALLOWED_GAP, MAX_ALLOWED_GAP, MAX_BSCORE, 0, VERBOSITY);

    if (VERBOSITY>=3) {
      for (int i = 0; i < nscaff; ++i) {
	cout << "Scaffold " << i << " links:" << endl;
	int v1 = 2*i;
	for (int j = 0; j < Gbundles.From(v1).isize(); ++j) {
	  int v2 = Gbundles.From(v1)[j];
	  vec<CLinkBundle> vb = Gbundles.EdgeObjectsBetween(v1,v2);
	  cout << "  fw>" << v2/2 << ((v2&1)?"rc":"fw") << " " << vb.size();
	  for (int k = 0; k < vb.isize(); ++k) cout << " " << vb[k].AsString() << " " << ToString(Bscore(vb[k]),2);
	  cout << endl;
	}
	++v1;
	for (int j = 0; j < Gbundles.From(v1).isize(); ++j) {
	  int v2 = Gbundles.From(v1)[j];
	  vec<CLinkBundle> vb = Gbundles.EdgeObjectsBetween(v1,v2);
	  cout << "  rc>" << v2/2 << ((v2&1)?"rc":"fw") << " " << vb.size();
	  for (int k = 0; k < vb.isize(); ++k) cout << " " << vb[k].AsString() << " " + ToString(Bscore(vb[k]),2);
	  cout << endl;
	}
      }
    }
    if (DOT || DOT_PER_ITER) {
      vec< vec<String> > labels;
      for (int v1 = 0; v1 < Gbundles.N(); ++v1) {
	int s1 = v1 / 2;
	String s1dir = (v1 & 1) ? "rc" : "fw";
	vec<String> v1labels;
	vec<int> from = Gbundles.From(v1);
	for (int i = 0; i < from.isize(); ++i) {
	  ostringstream label;
	  int v2 = from[i];
	  int s2 = v2 / 2;
	  String s2dir = (v2 & 1) ? "rc" : "fw";
	  CLinkBundle b = Gbundles.EdgeObjectsBetween(v1,v2)[0];
	  //label << s1 << s1dir << ">" << s2 << s2dir << "(" << b.AsString() << ")";
	  label << s1 << s1dir << ">" << s2 << s2dir << "(" << b.Weight() << "," << ToString(b.Score(),2) << "," << b.Sep() << ")";
	  v1labels.push_back(label.str());
	}
	labels.push_back(v1labels);
      }
      ostringstream dotfile;
      dotfile << SCAFFOLDS_OUT;
      if (DOT_PER_ITER) dotfile << "." << iter;
      dotfile << ".dot";
      Ofstream(dot, dotfile.str());
      Gbundles.DOT(dot, labels);
      dot.close();
    }

    if ( SCAFFOLD_GRAPH_OUT != "" ) {
      String out_file = SCAFFOLD_GRAPH_OUT + ".after." + ToString( iter );
      BinaryWriter::writeFile( out_file, Gbundles );
    }

    // Components, sources, etc
    vec< vec<int> > comp;
    Gbundles.Components(comp);
    vec<int> sources;
    Gbundles.Sources(sources);
    cout << Date( ) << ": " << Gbundles.N() << " nodes, " << Gbundles.EdgeObjectCount() << " edges, " << comp.size() << " components, " << sources.isize() << " sources" << endl;

    if ( !JOIN ) break;


    // "Lines" are lists of scaffolds to be joined together.
    // The triples are target,sep,dev (sep and dev are 0 for the
    // first element.  Should be made into something more
    // elegant, like CLinkBundles or its own class.
	  
    vec< vec< triple<int,int,int> > > lines;

    if (VERBOSITY>=1) cout << Date( ) << ": Walking" << endl;
    // Keep track of which scaffolds have been included (note this
    // is per-scaffold, not per node
    vec<Bool> used(nscaff, False);

#ifdef notdef
    // try to scaffold component by component. not implemented
    for (u_int c = 0; c < comp.size(); ++c) {
      vec< triple<int,int,int> > path;
      path = ScaffoldComponent(Gbundles, comp[c], used);
      if (path.isize() > 1) lines.push(path);
    }

#endif

    // First loop over the sources in the graph, going as far as we
    // comfortably can.  
    for (u_int s = 0; s < sources.size(); ++s) {
      vec< triple<int,int,int> > path = walkforward(sources[s], Gbundles, scaffolds, used);
      if (path.isize() > 1) lines.push(path);
    }

    // Now try everything else which has links but not used
    for (int n = 0; n < Gbundles.N(); ++n) {
      if (!used[n/2] && Gbundles.From(n).isize() > 0) {
	vec< triple<int,int,int> > path = walkforward(n, Gbundles, scaffolds, used);
	if (path.isize() > 1) lines.push(path);
      }
    }

    // Make sure we include scaffolds which weren't touched
    for (int i = 0; i < nscaff; ++i) {
      vec< triple<int,int,int> > line;
      if (!used[i]) {
	line.push(2*i, 0, 0);
	lines.push_back(line);
      }
    }
    // Merge the scaffolds in the lines to form new scaffolds.  This is a
    // purely formal operation.

    if (VERBOSITY >= 1)
      cout << Date( ) << ": merging scaffolds" << endl;

    vec<superb> new_scaffolds;
    vec<Bool> reversed( ntigs, False );
    if ( VERBOSITY >= 2 )
      {    cout << "\nlines:\n";
	for ( int i = 0; i < lines.isize( ); i++ )
	  {  cout << i << ": ";
	    for ( int j = 0; j < lines[i].isize( ); j++ )
	      {    if ( j > 0 ) cout << " --> ";
		int s = lines[i][j].first;
		cout << s/2
		     << ( s % 2 == 0 ? "fw" : "rc" );
	      }
	    cout << endl;    }

      }
    for ( int i = 0; i < lines.isize( ); i++ )
      {
	joins += lines[i].size( ) - 1;
	superb x;
	int ntigs = 0;
	for ( int j = 0; j < lines[i].isize( ); j++ )
	  ntigs += scaffolds[ lines[i][j].first/2 ].Ntigs( );
	x.SetNtigs(ntigs);
	int count = 0;
	for ( int j = 0; j < lines[i].isize( ); j++ )
	  {    int y = lines[i][j].first;
	    int s = y/2;
	    superb S = scaffolds[s];


	    // In the rc case, reverse S, and also track the change.

	    if ( y % 2 == 1 )
	      {    S.Reverse( );
		for ( int k = 0; k < S.Ntigs( ); k++ )
		  {    rctig[ S.Tig(k) ] = !rctig[ S.Tig(k) ];
		    reversed[ S.Tig(k) ] = True;
		  }    }

	    // Add S to the scaffold.

	    for ( int k = 0; k < S.Ntigs( ); k++ )
	      {
		x.SetTig( count, S.Tig(k) );
		x.SetLen( count, S.Len(k) );
		if ( k < S.Ntigs( ) - 1 )
		  {    x.SetGap( count, S.Gap(k) );
		    x.SetDev( count, S.Dev(k) );
		    count++;    }    }

	    // Add the between scaffold gap.
	    
	    if ( j < lines[i].isize( ) - 1 )
	      {    x.SetGap( count, lines[i][j+1].second );
		x.SetDev( count, lines[i][j+1].third );    }
	    count++;    }
		    	       
	new_scaffolds.push_back(x);
      }
    
    // Reverse the alignments lying over contigs that have just been reversed.
    for ( int i = 0; i < aligns0.isize( ); i++ )
      if ( reversed[ aligns0[i].TargetId( ) ] ) aligns0[i].Reverse( );


    scaffolds = new_scaffolds;
    cout << Date( ) << ": made " << joins << " joins on iter " << iter << endl;

    // Absorb small scaffolds into larger ones.
    if ( SUCK_SCAFFOLDS ) {
      ofstream devnull ( "/dev/null" );
      ostream *plog = ( VERBOSITY > 0 ? &cout : &devnull );
      int n_sucked = SuckScaffolds( pairs, aligns0_index, aligns0,
				    scaffolds, rctig, plog );
      cout << Date( ) << ": " << n_sucked << " scaffolds embedded" << endl;
    }

    vec<int> scaffold_sizes(scaffolds.size());
    for ( size_t i = 0; i < scaffolds.size(); i++ )
      scaffold_sizes[i] = scaffolds[i].ReducedLength( );
    Sort(scaffold_sizes);
    scaffold_sizes.EraseValue( 0 );
    cout << Date( ) << ": N50 scaffold size (without gaps) = " << N50(scaffold_sizes) << endl;

    // Check to see if done.
    ++iter;
    if (iter >= MAX_ITER) break;
    
  } while (joins > 0);

  // Flip rc contigs.
     
  cout << Date() << ": flipping rc contigs in the fastavector" << endl;
  for ( size_t i = 0; i < ntigs; i++ )
    if ( rctig[i] ) {
      scaffold_contigs[i].ReverseComplement();
    }
     
  // Report some stats.

  vec<int> contig_sizes, scaffold_sizes;
  for ( size_t i = 0; i < scaffolds.size(); i++ )
    {    for ( int j = 0; j < scaffolds[i].Ntigs( ); j++ )
	contig_sizes.push_back(  scaffolds[i].Len(j) );
      scaffold_sizes.push_back( scaffolds[i].ReducedLength( ) );    }
  Sort(contig_sizes), Sort(scaffold_sizes);
  contig_sizes.EraseValue( 0 );
  scaffold_sizes.EraseValue( 0 );
     
  if ( contig_sizes.empty() )
    cout << "There are no contigs.  Something is very wrong." << endl;
  else
    cout << "N50 contig size = " << N50(contig_sizes) 
	 << " (" << contig_sizes.size() << " contigs)" << endl;
  if ( scaffold_sizes.empty() )
    cout << "There are no scaffolds.  Something is very wrong." << endl;
  else
    cout << "N50 scaffold size (without gaps) = " << N50(scaffold_sizes) 
	 << " (" << scaffold_sizes.size() << " scaffolds)" << endl;
  cout << "total bases in contigs = " << BigSum(scaffold_sizes) << endl;

  // Save results.

  if ( !WRITE ) exit(0);
     
  cout << "\n" << Date( ) << ": writing output files" << endl;
     
  WriteSuperbs( scaffolds_out_file , scaffolds );
  Ofstream( contig_fasta, contigs_out_file );
  for ( size_t i = 0; i < ntigs; i++ )
    scaffold_contigs[i].Print( contig_fasta, "contig_" + ToString(i) );
  contig_fasta.close();
     
  WriteScaffoldedFasta( assembly_out_file, scaffold_contigs, scaffolds );

  cout << Date( ) << ": done" << endl;
}
