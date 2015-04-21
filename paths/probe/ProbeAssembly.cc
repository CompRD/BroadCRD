///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ProbeAssembly.  Given a range of contigs on a scaffold, gather information
// about them.  [EXPERIMENTAL, IN PROGRESS.]
//
// Note that the way this computes the unipath graph is buggy.
//
// Requirements:
// 1. You need to have "dot" in your path.  It is part of the Graphviz package.

// MakeDepend: dependency DisplayLocs
// MakeDepend: dependency LocalizeReadsLG
// MakeDepend: dependency MergeNeighborhoods1
// MakeDepend: dependency MergeNeighborhoods2
// MakeDepend: dependency MergeNeighborhoods3
// MakeDepend: dependency DumpHyper
// MakeDepend: library OMP
#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "paths/GetNexts.h"
#include "paths/UnibaseUtils.h"
#include "util/SearchFastb2Core.h"

// Sort aligns by: id1, start, id2
struct order_aligns_local
  : public binary_function< const triple<int64_t,int64_t,int> &,
			    const triple<int64_t,int64_t,int> &,
			    bool >
{
public:
  bool operator( ) ( const triple<int64_t,int64_t,int> &left,
		     const triple<int64_t,int64_t,int> &right ) const {
    if ( left.first != right.first )
      return ( left.first < right.first );
    if ( left.third != right.third )
      return ( left.third < right.third );
    return ( left.second < right.second );
  }
};

// Map local unibases onto global unibases.
void MapLocalToGlobal ( const int K,
			const String &fasta_file,
			const String &map_file,
			const String &unibases_file ) 
{
  // Stream for map, and file names.
  ofstream map_out( map_file.c_str( ) );
  String kmers_file = fasta_file + ".kmers";

  // Load fasta and names, save fastb of all K-mers.
  vec<String> names;
  vec<fastavector> fastas;
  LoadFromFastaFile( fasta_file, fastas, names );
       
  int nkmers = 0;
  vec<int> first_kmer( fastas.size( ), -1 );
  for (int ii=0; ii<fastas.isize( ); ii++) {
    first_kmer[ii] = nkmers;
    nkmers += fastas[ii].size( ) - K + 1;
  }
       
  vecbvec kmers;
  kmers.reserve( nkmers );
  for (size_t ii=0; ii<fastas.size( ); ii++) {
    bvec fastb = fastas[ii].ToBasevector( );
    for (int jj=0; jj<fastb.isize( ) - K + 1; jj++)
      kmers.push_back( bvec( fastb, jj, K ) );
  }
       
  kmers.WriteAll( kmers_file );
       
  // Align and sort aligns (remove rc's).
  vec< triple<int64_t,int64_t,int> > aligns;
  {
    SearchFastb2( kmers_file, unibases_file, K, &aligns, 0, -1, 0.9, False );
	 
    vec< triple<int64_t,int64_t,int> > fwaligns;
    fwaligns.reserve( aligns.size( ) / 2 );
    for (size_t ii=0; ii<aligns.size( ); ii++)
      if ( aligns[ii].third > -1 )
	fwaligns.push_back( aligns[ii] );
	 
    order_aligns_local sorter;
    sort( fwaligns.begin( ), fwaligns.end( ), sorter );
	 
    swap( fwaligns, aligns );
  }
       
  // Clean up.
  Remove( kmers_file );

  // Map kmer pos to align.
  vec<int> to_align( kmers.size( ), -1 );
  for (int ii=0; ii<aligns.isize( ); ii++)
    to_align[ aligns[ii].first ] = ii;
       
  // Loop over all fastas.
  for (int fasta_id=0; fasta_id<fastas.isize( ); fasta_id++) {
    const fastavector &fasta = fastas[fasta_id];
    bool is_last = fasta_id == fastas.isize( ) - 1;
    int kbeg = first_kmer[fasta_id];
    int kend = is_last ? kmers.size( ) : first_kmer[fasta_id+1];
	 
    // Loop over all kmers in this fasta.
    vec<int> uids;
    bool missing_kmers = false;
    for (int kid=kbeg; kid<kend; kid++) {
      int align_id = to_align[kid];
      if ( align_id < 0 ) {
	missing_kmers = true;
	continue;
      }
      int uid = aligns[align_id].second;
      if ( uids.size( ) < 1 || uids.back( ) != uid )
	uids.push_back( uid );
    }
	 
    for (int ii=0; ii<uids.isize( ); ii++)
      map_out << BaseAlpha(fasta_id) << "\t" << uids[ii] << "\n";

    if ( missing_kmers )
      cout << " WARNING! There are missing kmers in fasta_" << fasta_id
	   << " from local assembly " << fasta_file << endl;
  }
  
  // Close map stream.
  map_out.close( );
}

void DotToPng( const String& fn )
{    SystemSucceed( "cat " + fn + " | dot -Tpng > " + fn + ".png" );    }

int main(int argc, char *argv[])
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(PROBEDIR, "probe");
     CommandArgument_String_OrDefault(RESULTSDIR, "");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_Int(K);
     CommandArgument_Int(TIG1);
     CommandArgument_Int(TIG2);
     CommandArgument_Bool_OrDefault(LOCS, True);
     CommandArgument_Bool_OrDefault(PNG, True);
     CommandArgument_Int_OrDefault(ADDS, 6);
     // LocalizeReadsLG options
     CommandArgument_Bool_OrDefault(ASSEMBLE, True);
     CommandArgument_Bool_OrDefault(BRIDGE_BETWEEN_CLOUD_UNIPATHS, False);
     // MergeNeighborhoods options
     CommandArgument_Bool_OrDefault(MERGE_NHOODS, False);
     CommandArgument_Int_OrDefault(MN_MIN_ALIGN_LENGTH, 500);
     CommandArgument_Int_OrDefault(MN_MIN_OVERLAP, 10000);
     CommandArgument_Int_OrDefault(MN_MIN_PROPER_OVERLAP, 2000);
     CommandArgument_Int_OrDefault(MN_MIN_PROPER_OVERLAP_FINAL, 5000);
     EndCommandArguments;

     // Define directories.

     PROBEDIR = SUBDIR + "/" + PROBEDIR;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String probe_dir = run_dir + "/ASSEMBLIES/" + PROBEDIR;
     String results_dir = (RESULTSDIR != "" ? RESULTSDIR : probe_dir + "/results");
     Mkdir777(probe_dir);
     Mkdir777(results_dir);
     cout << "probe directory: " << probe_dir << endl;
     cout << "results directory: " << results_dir << endl << endl;

     // Load scaffolds, check TIG1 and TIG2, output summary.

     int ntigs
          = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
     cout << Date( ) << ": loading scaffolds" << endl;
     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    int n = scaffolds[i].Ntigs( );
          for ( int j = 0; j < n; j++ )
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;    }    }
     int s1 = to_super[TIG1], s2 = to_super[TIG2];
     const superb& S = scaffolds[s1];
     int p1 = to_super_pos[TIG1], p2 = to_super_pos[TIG2];
     if ( !( s1 == s2 && p1 <= p2 ) )
     {    cout << "TIG1 and TIG2 don't make sense" << endl;
          exit(1);    }
     Ofstream( sout, results_dir + "/summary" );
     for ( int j = p1; j <= p2; j++ )
     {    int m = S.Tig(j);
          sout << m << " (l = " << S.Len(j) << ")\n";
          if ( j < p2 )
               sout << " -- (" << S.Gap(j) << " +/- " << S.Dev(j) << ") --> ";    }

     // Load unibases and generate ancillary structures.

     cout << Date( ) << ": loading unibases" << endl;
     vecbasevector unibases( run_dir + "/all_reads.unibases.k" + ToString(K) );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc );
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );

     // Load seeds.

     vec<int> seeds;
     fast_ifstream in( sub_dir + "/seeds.ids" );
     String line;
     getline( in, line );
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          seeds.push_back( line.Int( ) );    }

     // Go through the contigs.

     vec<int> seeds_used, all_unis;
     Ofstream( aout, results_dir + "/seed_map" );
     for ( int pp = p1; pp <= p2; pp++ )
     {    int tig = scaffolds[s1].Tig(pp);
          aout << endl;
          PRINT3_TO( aout, s1, p2, tig );
	  cout << Date( ) << ": dealing with contig " << pp
	       << " (" << pp - p1 + 1 << " of " << p2 - p1 + 1 << ")" << endl;

          // Run DisplayLocs on the contig.
	  
          if (LOCS)
          {    SystemSucceed( "DisplayLocs" + ARGC(PRE) + ARGC(DATA) + ARGC(RUN)
                    + ARGC(SUBDIR) + ARGC(ASSEMBLY) + ARG(TIGS, tig) + ARG(NH, True)
                    + " > " + results_dir + "/tig." + ToString(tig) + ".locs" );    }

          // Align the contig to the unibases.  

          vecbasevector kmers, atig;
          atig.ReadOne( sub_dir + "/" + ASSEMBLY + ".contigs.fastb", tig );
          for ( int j = 0; j <= atig[0].isize( ) - K; j++ )
          {    basevector b;
               b.SetToSubOf( atig[0], j, K );
               kmers.push_back(b);    }
          kmers.WriteAll( results_dir + "/kmers" );
          vec< triple<int64_t,int64_t,int> > ALIGNS;
          SearchFastb2( results_dir + "/kmers", 
               run_dir + "/all_reads.unibases.k" + ToString(K), K, &ALIGNS, 0,
               -1, 0.90, False );
          Remove( results_dir + "/kmers" );
          vec<int> unis;
          for ( int i = 0; i < ALIGNS.isize( ); i++ )
               if ( ALIGNS[i].third >= 0 ) unis.push_back( ALIGNS[i].second );
          UniqueSort(unis);
          for ( int j = 0; j < unis.isize( ); j++ )
               aout << "hits unipath " << unis[j] << "\n";
          all_unis.append(unis);

          // See if any of the unibases are seeds.
     
          aout << "\n";
          for ( size_t i = 0; i < unis.size( ); i++ )
          {    int u = unis[i];
               int p = Position( seeds, u ), prc = Position( seeds, to_rc[u] );
               if ( p >= 0 ) 
               {    aout << "unipath " << u << " is seed " << p << endl;
                    seeds_used.push_back(p);    }
               if ( prc >= 0 )
               {    aout << "reverse complement of unipath " << u << " (" 
                         << to_rc[u] << ") is seed " << prc << endl;    
	       seeds_used.push_back(prc);    }    }    }

     // Go through passes in which we add successively more adjacent unibases.

     UniqueSort(all_unis);
     vec<int> all_unis_orig(all_unis);
     cout << Date( ) << ": starting addition passes" << endl;
     for ( int opass = 0; opass <= ADDS; opass++ )
     {    all_unis = all_unis_orig;
     
          // Add adjacent unibases, multiple times.

          for ( int pass = 1; pass <= opass; pass++ )
          {    int nu = all_unis.size( );
               for ( int i = 0; i < nu; i++ )
               {    int u = all_unis[i];
                    for ( int j = 0; j < nexts[u].isize( ); j++ )
                    {    int v = nexts[u][j];
                         all_unis.push_back(v);    }
                    for ( int j = 0; j < nexts[ to_rc[u] ].isize( ); j++ )
                    {    int v = to_rc[ nexts[ to_rc[u] ][j] ];
                         all_unis.push_back(v);    }    }
               UniqueSort(all_unis);    }
     
          // Dump these unibases.

          Ofstream( uout, 
               results_dir + "/selected_unibases." + ToString(opass) + ".fasta" );
          for ( size_t i = 0; i < all_unis.size( ); i++ )
          {    int u = all_unis[i];
               unibases[u].Print( uout, "unibase_" + ToString(u) );    }

          // Generate graph that shows the relationship between the unibases that
          // have been found.  The way we define vertices is terrible.  IF TWO
          // EDGES GO UP TO A VERTEX BUT NO EDGES BEYOND IT ARE PRESENT IN THE
          // SET, THEN THE VERTEX WILL BE SHOWN AS TWO VERTICES!

          UniqueSort(all_unis);
          vec<int> edges(all_unis);
          vec< pair<int,int> > verts0;
          for ( int i = 0; i < edges.isize( ); i++ )
          {    verts0.push( i, 0 );
               verts0.push( i, 1 );    }
          equiv_rel e( verts0.size( ) );
          for ( int i1 = 0; i1 < edges.isize( ); i1++ )
          {    int u = edges[i1];
               for ( int j = 0; j < nexts[u].isize( ); j++ )
               {    int v = nexts[u][j];
                    int i2 = Position( edges, v );
                    if ( i2 >= 0 ) e.Join( 2*i1 + 1, 2*i2 );    }    }
          vec<int> reps;
          e.OrbitRepsAlt(reps);
          int N = reps.size( );
          vec< vec<int> > from(N), to(N), from_edge_obj(N), to_edge_obj(N);
          for ( int i = 0; i < edges.isize( ); i++ )
          {    int u = edges[i];
               int x = BinPosition( reps, e.ClassId( 2*i ) );
               int y = BinPosition( reps, e.ClassId( 2*i + 1 ) );
               from[x].push_back(y), to[y].push_back(x);    
               from_edge_obj[x].push_back(i), to_edge_obj[y].push_back(i);    }
          vec<Bool> complete( N, True );
          for ( int i1 = 0; i1 < edges.isize( ); i1++ )
          {    int u = edges[i1];
               for ( int j = 0; j < nexts[u].isize( ); j++ )
               {    int v = nexts[u][j];
                    int i2 = Position( edges, v );
                    if ( i2 < 0 ) 
                    {    complete[ BinPosition( reps, e.ClassId( 2*i1 + 1 ) ) ] 
                              = False;
                         aout << "see " << u << " --> " << v << "\n";    }    }
               for ( int j = 0; j < nexts[ to_rc[u] ].isize( ); j++ )
               {    int v = to_rc[ nexts[ to_rc[u] ][j] ];
                    int i2 = Position( edges, v );
                    if ( i2 < 0 ) 
                    {    complete[ BinPosition( reps, e.ClassId( 2*i1 ) ) ] = False;
                         aout << "see " << v << " --> " << u << "\n";    }    }    }
          for ( int i = 0; i < N; i++ )
          {    SortSync( from[i], from_edge_obj[i] );
               SortSync( to[i], to_edge_obj[i] );    }
          digraphE<int> G( from, to, edges, to_edge_obj, from_edge_obj );
          vec<double> lengths( edges.size( ) );
          vec<String> edge_labels( edges.size( ) );
          for ( size_t i = 0; i < edges.size( ); i++ )
          {    edge_labels[i] = ToString( edges[i] );
               if ( BinMember( all_unis_orig, edges[i] ) )
                    edge_labels[i] += "*";    }
          for ( size_t i = 0; i < edges.size( ); i++ )
               lengths[i] = unibases[ edges[i] ].size( );
          {    Ofstream( dout, results_dir + "/main.dot.prelim" );

               G.PrettyDOT( dout, lengths, 
                    digraphE<int>::edge_label_info( digraphE<int>::edge_label_info::EXPLICIT, False, True,
                         &edge_labels ) );    }
               
          {    Ofstream( dout2, results_dir + "/main." + ToString(opass) + ".dot" );
               fast_ifstream din( results_dir + "/main.dot.prelim" );
               for ( int j = 0; j < 7; j++ )
               {    getline( din, line );
                    dout2 << line << "\n";    }
               for ( int i = 0; i < N; i++ )
               {    if ( complete[i] ) continue;
                    dout2 << i 
                         << " [width=0.1,height=0.1,shape=circle,fillcolor=white,"
                         << "label=\"\"];\n";    }
               while(1)
               {    getline( din, line );
                    if ( din.fail( ) ) break;
                    dout2 << line << "\n";    }    }
          if (PNG) DotToPng( results_dir + "/main." + ToString(opass) + ".dot" );
          Remove( results_dir + "/main.dot.prelim" );    }
     
     // Rerun local assemblies.

     if ( !ASSEMBLE ) return 0;
     cout << Date( ) << ": running local assemblies" << endl;
     if ( seeds_used.empty( ) ) return 0;
     UniqueSort(seeds_used);
     {    Ofstream( sout, results_dir + "/seeds_used" );
          for ( size_t i = 0; i < seeds_used.size( ); i++ )
               sout << seeds_used[i] << "\n";    }

     Cp2( results_dir + "/seeds_used", probe_dir + "/seeds_used");
     SymlinkForce( "../unipath_link_graph.cloud.k96", 
          probe_dir + "/unipath_link_graph.cloud.k96");
     SymlinkForce( "../seeds.ids", probe_dir + "/seeds.ids");

     int nthreads = omp_get_max_threads( );
     SystemSucceed( "LocalizeReadsLG" + ARGC(PRE) + ARGC(DATA) + ARGC(RUN)
		    + ARGC(K) + ARG(SUBDIR,PROBEDIR) 
		    + " READS=all_reads " + ARG(NUM_THREADS, nthreads)
          + ARG(USE_TRUTH, False) + ARG(LOCAL_DOT, True) + ARG(LOCAL_FASTA, True)
          + ARG(DUMP_LOCAL_UNIBASES, True) + ARG(LOCAL_LOG, True)
          + " OVERRIDE_TIMEOUTS=True MAX_COPY_NUMBER_OTHER=2 "
          + ARG(SEED_IDS, "@" + results_dir + "/seeds_used")
          + ARGC(BRIDGE_BETWEEN_CLOUD_UNIPATHS)
          + " > " + probe_dir + "/LocalizeReadsLG.log" );
     int xid = 0;
     while(1)
     {    String xdir = probe_dir + "/seed/" + ToString(xid++);
          if ( !IsDirectory(xdir) ) break;
          for ( size_t i = 0; i < seeds_used.size( ); i++ )
          {    int s = seeds_used[i];
               String head = xdir + "/" + ToString(s);
               String rhead = results_dir + "/" + ToString(s);
               if ( IsRegularFile( head + ".log" ) )
               {    Cp2( head + ".log", rhead + ".log" );
                    Cp2( head + ".hbv.fasta", rhead + ".hbv.fasta" );
                    Cp2( head + ".hbv.dot", rhead + ".hbv.dot" );
                    if (PNG) DotToPng( rhead + ".hbv.dot" );
                    Cp2( head + ".local_unibases.fasta", 
                         rhead + ".local_unibases.fasta" );
                    Cp2( head + ".local_unibases.dot", 
                         rhead + ".local_unibases.dot" );
                    // if (PNG) DotToPng( rhead + ".local_unibases.dot" ); // slow!
                         }    }    }

     // Merge local assemblies.

     if (MERGE_NHOODS) {
       cout << Date( ) << ": merging local assemblies" << endl;
       SystemSucceed( "MergeNeighborhoods1" 
            + ARGC(PRE) + ARGC(DATA) + ARGC(RUN) + ARGC(K) + ARG(SUBDIR, PROBEDIR) 
            + ARG(NUM_THREADS, nthreads) + " SEEDS=@" + probe_dir + "/seeds_used"
            + " > " + probe_dir + "/MergeNeighborhoods1.log" );
       SystemSucceed( "MergeNeighborhoods2" 
            + ARGC(PRE) + ARGC(DATA) + ARGC(RUN) + ARGC(K) + ARG(SUBDIR, PROBEDIR) 
            + ARG(NUM_THREADS, nthreads) + ARG(MIN_OVERLAP, MN_MIN_OVERLAP)
            + ARG(MIN_PROPER_OVERLAP, MN_MIN_PROPER_OVERLAP)
		      + " > " + probe_dir + "/MergeNeighborhoods2.log" );
       SystemSucceed( "MergeNeighborhoods3" 
            + ARGC(PRE) + ARGC(DATA) + ARGC(RUN) + ARGC(K) + ARG(SUBDIR, PROBEDIR)
            + " > " + probe_dir + "/MergeNeighborhoods3.log" );
       SystemSucceed( "DumpHyper" 
            + ARGC(PRE) + ARGC(DATA) + ARGC(RUN) + ARG(SUBDIR, PROBEDIR)
            + ARG(WRUN, "") + " > " + probe_dir + "/MergeNeighborhoods3.log" );
       Cp2( probe_dir + "/hyper.fasta", results_dir + "/hyper.fasta");
       Cp2( probe_dir + "/hyper.dot", results_dir + "/hyper.dot");       
     }
     
     // Map local unibases to global unibases.
     
     if ( MERGE_NHOODS ) {
       cout << Date( ) << ": mapping merged neighboorhoods" << endl;
       
       String fasta_file = results_dir + "/hyper.fasta";
       String map_file = results_dir + "/hyper.map";
       String unibases_file = run_dir + "/all_reads.unibases.k" + ToString( K );
       
       MapLocalToGlobal( K, fasta_file, map_file, unibases_file );
     }
       for (size_t select=0; select<seeds_used.size( ); select++) {
	 int seed_id = seeds_used[select];
	 cout << Date( ) << ": mapping neighboorhood " << seed_id << endl;
	 
	 String str_seed = ToString( seed_id );
	 String fasta_file = results_dir + "/" + str_seed + ".hbv.fasta";
	 String map_file = results_dir + "/" + str_seed + ".hbv.map";
	 String unibases_file = run_dir + "/all_reads.unibases.k" + ToString( K );
	 
	 MapLocalToGlobal( K, fasta_file, map_file, unibases_file );
       }

     cout << Date( ) << ": done" << endl;    }
