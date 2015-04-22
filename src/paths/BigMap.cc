///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BigMap.  Build a map of the genome.
//
// Note a problem with the handling of circular scaffolds: we don't distinguish
// between the case where there's a gap at the end, and the case where there's
// contiguous circular sequence.  Also we don't save any knowledge about 
// circularity.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmerRecord.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/AlignReadsOnUnibases.h"
#include "paths/AssemblyCleanupTools.h"
#include "fastg/FastgTools.h"
#include "paths/BigMapDot.h"
#include "paths/BigMapTools.h"
#include "paths/GetHomes.h"
#include "paths/GetNexts.h"
#include "paths/LongReadTools.h"
#include "paths/PdfEntry.h"
#include "paths/ProcessGap.h"
#include "paths/RemodelGapTools.h"
#include "paths/UnibaseUtils.h"
#include "paths/UnipathScaffold.h"
#include "paths/Uniseq.h"
#include "paths/long/CreateGenome.h"

class hit_thingee {

     public:

     hit_thingee( ) { }
     hit_thingee( const int u2, const int pos2, const int id2,
          const int dist_from_end1, const Bool unique1 ) 
          : u2(u2), pos2(pos2), id2(id2),
          dist_from_end1(dist_from_end1), unique1(unique1) { }

     int u2;
     int pos2;
     int id2;
     int dist_from_end1;
     Bool unique1;

};

void FilterByJumps( snark& S, const int K2, const vecbasevector& jbases,
     const PairsManager& jpairs, 
     const vec< vec< triple<int,int,Bool> > >& placements_by_read,
     const vec< vec< triple<int64_t,int,Bool> > >& placements_by_unipath,
     const int VERBOSITY, const int FILTER_VERBOSITY,
     const int FILTER_VERBOSE_U1, const Bool FILTER_SHOW_BASES_AFTER ){    
  while(1){    
    cout << Date( ) << ": filtering walks using jumps2" << endl;
    vec<int> to_left, to_right;
    S.G( ).ToLeft(to_left), S.G( ).ToRight(to_right);
    int delta = 0;
    for ( int edge_id = 0; edge_id < S.EdgeN( ); edge_id++ ){    
      if ( S.Edge(edge_id).Closed( ) ){    
	int closures_before = S.Edge(edge_id).ClosureCount( );
	int verbosity = FILTER_VERBOSITY;
	if ( FILTER_VERBOSE_U1 >= 0 ){    
	  if ( S.Vert( to_left[edge_id] ).U( ).back( ) != FILTER_VERBOSE_U1 ){    
	    verbosity = 0;    }    }
	FilterWalksUsingJumps2( edge_id, S, to_left, to_right, jbases, jpairs, 
				placements_by_read, placements_by_unipath, verbosity );    
                    int closures_after = S.Edge(edge_id).ClosureCount( );
                    delta += closures_before - closures_after;    
      }    
    }
    S.SwallowSimpleGaps( );
    S.BringOutTheDead( );
    cout << Date( ) << ": " << delta << " closures eliminated" << endl;
    if ( delta == 0 ) break;    
  }
  for ( int x = 0; x < S.VertN( ); x++ ){    
    for ( int j = 0; j < S.G( ).From(x).isize( ); j++ ){    
      int y = S.G( ).From(x)[j];
      const gapster& g = S.G( ).EdgeObjectByIndexFrom( x, j );
      if ( g.Open( ) ) continue;
      if ( VERBOSITY >= 1 ){    
	cout << "\nafter filtering, walks from " << S.Vert(x).U( ).back( ) 
	     << " to " << S.Vert(y).U( ).front( ) << "\n";
	for ( int l = 0; l < g.ClosureCount( ); l++ ){    
	  cout << "[" << l << "] ";
	  g.Closure(l).Print( cout, K2 );
	  cout << "\n";
	  if (FILTER_SHOW_BASES_AFTER)
	    g.Closure(l).Bases( ).Print( cout, l );      
	}    
      }    
    }    
  }    
}

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, "long_direct");
     CommandArgument_Int_OrDefault_Doc(K1, 96, "little K");
     CommandArgument_Int_OrDefault_Doc(K2, 640, "big K");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_Doc(HEAD1, "head for K1 unibases");
     CommandArgument_String_Doc(HEAD2, "head for K2 unibases");
     CommandArgument_String_OrDefault(JUMPS_IN, "jump_reads_ec");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault(WRITE_SNARK, True);
     CommandArgument_Int_OrDefault(FORCE_FIRST, -1);
     CommandArgument_Double_OrDefault(MAX_DEV_DIFF, 3.0);
     CommandArgument_Int_OrDefault(MIN_LINKS_FOR_UNBRIDGEABLE, 0);
     CommandArgument_Bool_OrDefault(DELETE_COMPLEMENTS, True);
     CommandArgument_Bool_OrDefault_Doc(REMOVE_DUPLICATES, False, 
        "Remove duplicate (or almost duplicate) scaffolds");
     
     // logging options

     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     CommandArgument_Bool_OrDefault(VALIDATE, False);
     CommandArgument_Bool_OrDefault(AUTO_VALIDATE, False);
     CommandArgument_Bool_OrDefault(SHOW_TRANSLATION, False);
     CommandArgument_Bool_OrDefault(SHOW_LINKS, False);
     CommandArgument_Bool_OrDefault(SHOW_LINKS2, False);
     CommandArgument_Bool_OrDefault(SHOW_INITIAL_LINKS, False);
     CommandArgument_Bool_OrDefault(SHOW_INITIAL_LINKS2, False);
     CommandArgument_Bool_OrDefault(PRINT_GRAPH1, False);
     CommandArgument_Bool_OrDefault(PRINT_GRAPH2, False);
     CommandArgument_Bool_OrDefault(SCAFFOLDING_VERBOSE, False);
     CommandArgument_Bool_OrDefault(RAISE_VERBOSE, False);
     CommandArgument_Bool_OrDefault(DELETE_VERBOSE, False);
     CommandArgument_String_OrDefault(DOT, "");
     CommandArgument_Bool_OrDefault(DOT_SIMPLE, False);
     CommandArgument_Bool_OrDefault(CIRCO, True);
     CommandArgument_String_OrDefault(DOT_LEGENDS,
        "{Legend,Summary,Discrepancies}");
     CommandArgument_Bool_OrDefault(AUTO_DOT, False);
     CommandArgument_Int_OrDefault(FILTER_VERBOSITY, 0);
     CommandArgument_Int_OrDefault(FILTER_VERBOSE_U1, -1);
     CommandArgument_Bool_OrDefault(FILTER_SHOW_BASES_AFTER, False);
     CommandArgument_Bool_OrDefault(SHOW_PARTNERS_IN_GAP, False);
     CommandArgument_Bool_OrDefault(DUMP, False);
     CommandArgument_Bool_OrDefault(DUMP_VERBOSE, False);
     CommandArgument_Bool_OrDefault(DUMP_PAIRS, False);
     CommandArgument_Bool_OrDefault(ANALYZE_EDGES, False);
     CommandArgument_Bool_OrDefault(UNWIND, False);
     CommandArgument_Int_OrDefault(MAX_WALKS_TO_PRINT, 20);
     CommandArgument_Int_OrDefault(MIN_KMERS1, 800);

     EndCommandArguments;

     // Thread control, etc.

     double clock = WallClockTime( );
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Heuristics.

     const int min_kmers2 = 100;

     // Define directories.

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     Mkdir777( run_dir + "/tmp" );

     // Parse funny arguments.

     if (AUTO_DOT)
     {    if ( DATA.Contains( "/babies/" ) 
               && DATA.After( "/babies/" ).Contains( "/" ) )
          {    DOT = DATA.Between( "/babies/", "/" ) + ".dot";    }    }
     if ( AUTO_VALIDATE && IsRegularFile( data_dir + "/genome.fastb" ) )
          VALIDATE = True;

     // Load unibases and generate ancillary data.

     cout << Date( ) << ": loading unibases" << endl;
     String KS1 = ToString(K1), KS2 = ToString(K2);
     vecbasevector unibases1( run_dir + "/" + HEAD1 + ".unibases.k" + KS1 );
     String unibases2_fn = run_dir + "/" + HEAD2 + ".unibases.k" + KS2;
     vecbasevector unibases2(unibases2_fn);
     int nuni1 = unibases1.size( ), nuni2 = unibases2.size( );
     vec<int> to_rc1, to_rc2;
     UnibaseInvolution( unibases1, to_rc1 );
     UnibaseInvolution( unibases2, to_rc2 );
     vec< vec<int> > nexts1;
     GetNexts( K1, unibases1, nexts1 );

     // Unwind unibases2.

     const int L = 96;
     if (UNWIND)
     {    cout << "\nunwinding of unibases2 onto unibases1:\n";
          vec< triple<kmer<L>,int,int> > kmers_plusx;
          MakeKmerLookup0( unibases1, kmers_plusx );
          vec< kmer<L> > kmersx( kmers_plusx.size( ) );
          for ( size_t i = 0; i < kmersx.size( ); i++ )
               kmersx[i] = kmers_plusx[i].first;
          for ( int u = 0; u < nuni2; u++ )
          {    cout << u << " -->";
               int last_u1 = -1;
               for ( int p = 0; p <= unibases2[u].isize( ) - K1; p++ )
               {    kmer<L> x;
                    x.SetToSubOf( unibases2[u], p );
                    int m = BinPosition( kmersx, x );
                    ForceAssertGe( m, 0 );
                    int u1 = kmers_plusx[m].second;
                    if ( u1 != last_u1 ) cout << " " << u1;
                    last_u1 = u1;    }
               cout << "\n";    }
          cout << "\n";    }

     // Define alignments between unibases2.  There are two sets.  The first set
     // set (nexts2) consists of K2-1 base alignments.  The second set (nexts2x) 
     // consists of alignments of size between K1-1 and K2-1.

     vec< vec<int> > nexts2;
     GetNexts( K2, unibases2, nexts2 );
     const int K0 = 80;
     ForceAssertLe( K0, K1 );
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases2.size( ); i++ )
     {    const basevector& u = unibases2[i];
          starts.push_back( starts.back( ) + u.isize( ) - K0 + 1 );    }
     vec< triple<kmer<K0>,int,int> > kmers_plus( starts.back( ) );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( size_t i = 0; i < unibases2.size( ); i++ )
     {    const basevector& u = unibases2[i];
          for ( int j = 0; j <= u.isize( ) - K0; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf( u, j );
               kmers_plus[r].second = i, kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);
     vec< kmer<K0> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     cout << Date( ) << ": mapping from 2 to 2" << endl;
     vec< vec< pair<int,int> > > nexts2x( unibases2.size( ) );
     #pragma omp parallel for
     for ( int u1 = 0; u1 < (int) unibases2.size( ); u1++ )
     {    kmer<K0> x;
          x.SetToSubOf( unibases2[u1], unibases2[u1].isize( ) - K0 );
          int64_t low = LowerBound( kmers, x ), high = UpperBound( kmers, x );
          for ( int64_t z = low; z < high; z++ )
          {    int u2 = kmers_plus[z].second, pos2 = kmers_plus[z].third;
               if ( u2 == u1 && pos2 == unibases2[u1].isize( ) - K0 ) continue;
               int overlap = pos2 + K0;
               if ( overlap < K1 || overlap > unibases2[u1].isize( ) ) continue;
               Bool match = True;
               for ( int j = 0; j < overlap - K0; j++ )
               {    if ( unibases2[u2][j] 
                         != unibases2[u1][ unibases2[u1].isize( ) - overlap + j ] )
                    {    match = False;
                         break;    }    }
               if (match) nexts2x[u1].push( u2, overlap );    }    }
     vec< vec< pair<int,int> > > nexts2y(nexts2x);
     for ( int i = 0; i < nexts2y.isize( ); i++ )
     {    vec<Bool> to_remove( nexts2y[i].size( ), False );
          for ( int j = 0; j < nexts2y[i].isize( ); j++ )
          {    if ( nexts2y[i][j].second < K2 - 1 )
                    to_remove[j] = True;    }
          EraseIf( nexts2y[i], to_remove );    }

     // Clean nexts2x.

     int total_nextsx = 0, total_nextsy = 0;
     for ( int u = 0; u < nuni2; u++ )
     {    total_nextsx += nexts2x[u].size( );
          total_nextsy += nexts2y[u].size( );    }
     int killed = 0;
     for ( int u1 = 0; u1 < nuni2; u1++ )
     {    for ( int j1 = 0; j1 < nexts2x[u1].isize( ); j1++ )
          {    int u3 = nexts2x[u1][j1].first;
               int over13 = nexts2x[u1][j1].second;
               int start3 = unibases2[u1].isize( ) - over13;
               for ( int j2 = 0; j2 < nexts2x[u1].isize( ); j2++ )
               {    int u2 = nexts2x[u1][j2].first;
                    int over12 = nexts2x[u1][j2].second;
                    int start2 = unibases2[u1].isize( ) - over12;
                    for ( int j3 = 0; j3 < nexts2x[u2].isize( ); j3++ )
                    {    if ( nexts2x[u2][j3].first != u3 ) continue;
                         int over23 = nexts2x[u2][j3].second;
                         if ( start2 + unibases2[u2].isize( ) - over23 != start3 )
                              continue;
                         nexts2x[u1].erase( nexts2x[u1].begin( ) + j1 );
                         killed++;
                         int u1rc = to_rc2[u1], u3rc = to_rc2[u3];
                         if ( u1rc != u3 || u3rc != u1 )
                         {    for ( int k = 0; k < nexts2x[u3rc].isize( ); k++ )
                              {    if ( nexts2x[u3rc][k].first != u1rc ) continue;
                                   if ( nexts2x[u3rc][k].second != over13 ) continue;
                                   nexts2x[u3rc].erase( nexts2x[u3rc].begin( ) + k );
                                   killed++;
                                   break;    }    }    }    }    }    }

     // Locate unibases on genome.

     const int LG = 12;
     VecIntPairVec Glocs;
     vecbasevector genome;
     vec< vec<placementy> > Ulocs1, Ulocs2;
     if ( AUTO_VALIDATE && IsRegularFile( data_dir + "/genome.fastb" ) )
          VALIDATE = True;
     if (VALIDATE)
     {    cout << Date( ) << ": locating on genome" << endl;
          genome.ReadAll( data_dir + "/genome.fastb" );
          CreateGlocs( genome, LG, Glocs );
          Ulocs1.resize(nuni1);
          for ( int u = 0; u < nuni1; u++ )
          {    Ulocs1[u] = FindGenomicPlacementsY( u, unibases1[u], LG, genome, 
                    Glocs );    }
          Ulocs2.resize(nuni2);

          // Find perfect placements of unibases2.  Then use QueryLookupTable
          // to find imperfect placements.

          vecbasevector genome2( genome.size( ) );
          for ( size_t g = 0; g < genome.size( ); g++ )
               genome2[g] = Cat( genome[g], genome[g] );
          VecIntPairVec Glocs2;
          CreateGlocs( genome2, LG, Glocs2 );
          int count = 0;
          {    Ofstream( out, run_dir + "/tmp/BigMap.ids" );
               for ( int u = 0; u < nuni2; u++ )
               {    Ulocs2[u] = FindGenomicPlacementsY( 
                         u, unibases2[u], LG, genome2, Glocs2 );    
                    if ( Ulocs2[u].empty( ) ) 
                    {    count++;
                         out << u << "\n";    }    
                    vec<placementy> newlocs;
                    for ( int j = 0; j < Ulocs2[u].isize( ); j++ )
                    {    placementy p = Ulocs2[u][j];
                         if ( p.pos >= genome[p.g].isize( ) ) continue;
                         if ( p.Pos > genome[p.g].isize( ) )
                              p.Pos -= genome[p.g].isize( );  
                         newlocs.push_back(p);    }
                    Ulocs2[u] = newlocs;    }    }
          if ( count > 0 )
          {    // Technically unsound:
               if ( !IsRegularFile( data_dir + "/genome2.lookup" ) )
               {    genome2.WriteAll( data_dir + "/genome2.fastb" );
                    SystemSucceed( "MakeLookupTable SOURCE=" + data_dir 
                         + "/genome2.fastb OUT_HEAD=" + data_dir + "/genome2 "
                         + "LO=True" );    }
               fast_pipe_ifstream qlt( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS="
                    + run_dir + "/" + HEAD2 + ".unibases.k" + KS2 + " L=" + data_dir 
                    + "/genome2.lookup SEQS_IS_FASTB=True PARSEABLE=True "
                    + "SEQS_TO_PROCESS=@" + run_dir + "/tmp/BigMap.ids" );
               String line;
               look_align la;
               while(1)
               {    getline( qlt, line );
                    if ( qlt.fail( ) ) break;
                    if ( line.Contains( "QUERY", 0 ) )
                    {    la.ReadParseable(line);
                         int g = la.target_id;
                         int pos2 = la.pos2( ), Pos2 = la.Pos2( );
                         if ( pos2 >= genome[g].isize( ) ) continue;
                         if ( Pos2 > genome[g].isize( ) ) Pos2 -= genome[g].isize( );
                         Ulocs2[la.query_id].push( la.query_id, g, pos2, Pos2, 
                              la.Fw1( ), la.mutations, la.indels, la.pos1( ) 
                              + la.query_length - la.Pos1( ) );    }    }    }
          Remove( run_dir + "/tmp/BigMap.ids" );    }

     // Get copy number.

     VecPdfEntryVec CN1( ( run_dir + "/" + HEAD1 + ".unipaths.predicted_count.k"
          + KS1 ).c_str( ) );
     vec<int> predicted_CN1( nuni1, -1 );
     for ( int i = 0; i < nuni1; i++ )
          GetMostLikelyValue( predicted_CN1[i], CN1[i] );
     vec<double> CN_raw;
     BinaryReader::readFile( run_dir + "/" + HEAD1 + ".unipaths.cn_raw.k" + KS1,
                             &CN_raw );

     // Define copy number for each position in each unibase2.

     vec< vec<double> > cn2;
     String cnpos2_fn = run_dir + "/cnpos2";
     if ( IsRegularFile(cnpos2_fn) && IsOlder( unibases2_fn, cnpos2_fn ) )
          BinaryReader::readFile( cnpos2_fn.c_str( ), &cn2 );
     else
     {    cout << Date( ) << ": compute position copy number" << endl;
          ForceAssertEq( L, K1 );
          vec< triple<kmer<L>,int,int> > kmers_plusx;
          MakeKmerLookup( unibases1, kmers_plusx );
          vec< kmer<L> > kmersx( kmers_plusx.size( ) );
          for ( size_t i = 0; i < kmersx.size( ); i++ )
               kmersx[i] = kmers_plusx[i].first;
          cn2.resize(nuni2);
          for ( int u = 0; u < nuni2; u++ )
          {    cn2[u].resize( unibases2[u].size( ), -1.0 );
               for ( int p = 0; p <= unibases2[u].isize( ) - K1; p++ )
               {    kmer<L> x;
                    x.SetToSubOf( unibases2[u], p );
                    int m = BinPosition( kmersx, x );
                    if ( m < 0 ) continue;
                    cn2[u][p] = CN_raw[ kmers_plusx[m].second ];    }    }
          BinaryWriter::writeFile( cnpos2_fn.c_str( ), cn2 );    }

     // Load jumps.

     PairsManager jpairs;
     cout << Date( ) << ": loading jump data" << endl;
     vecbasevector jbases( run_dir + "/" + JUMPS_IN + ".fastb" );
     jpairs.Read( run_dir + "/" + JUMPS_IN + ".pairs" );
     jpairs.makeCache( );

     // Generate jbases_sorted and jbases_sorted_id.

     vec<basevector> jbases_sorted;
     vec<int64_t> jbases_sorted_id;
     int min_read = 1000000000;
     for ( size_t id = 0; id < jbases.size( ); id++ )
          min_read = Min( min_read, jbases[id].isize( ) );
     DPRINT(min_read);
     cout << Date( ) << ": sorting jbases" << endl;
     jbases_sorted.resize( jbases.size( ) );
     jbases_sorted_id.resize( jbases.size( ) );
     for ( size_t id = 0; id < jbases.size( ); id++ )
     {    jbases_sorted[id] = jbases[id];
          if ( jbases[id].isize( ) > min_read )
          {    jbases_sorted[id].SetToSubOf( 
                    jbases[id], jbases[id].isize( ) - min_read, min_read );    }
          jbases_sorted_id[id] = id;    }
     ParallelSortSync( jbases_sorted, jbases_sorted_id );

     // Compute placements_by_read.  

     vec< vec< triple<int,int,Bool> > > placements_by_read( jbases.size( ) );
     cout << Date( ) << ": looking up unibases" << endl;
     vec< triple<int64_t,int,int> > jaligns0;
     #pragma omp parallel for
     for ( size_t m = 0; m < unibases2.size( ); m++ )
     {    vec< triple<int64_t,int,int> > jaligns0m;
          for ( int p = 0; p <= unibases2[m].isize( ) - min_read; p++ )
          {    basevector b( unibases2[m], p, min_read );
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) b.ReverseComplement( );
                    int64_t low = LowerBound( jbases_sorted, b );
                    int64_t high = UpperBound( jbases_sorted, b );
                    for ( int64_t l = low; l < high; l++ )
                    {    int64_t id = jbases_sorted_id[l];
                         int mx = m;
                         if ( pass == 2 ) mx = -m-1;
                         jaligns0m.push( id, mx, p );    }    }    }
          #pragma omp critical
          {    jaligns0.append(jaligns0m);    }    }
     cout << Date( ) << ": sorting aligns" << endl;
     ParallelSort(jaligns0);
     vec< vec< triple<int64_t,int,Bool> > > 
          placements_by_unipath( unibases2.size( ) );
     cout << Date( ) << ": generating placements_by_read" << endl;
     for ( size_t j = 0; j < jaligns0.size( ); j++ )
     {    int u = jaligns0[j].second;
          Bool fw = True;
          if ( u < 0 )
          {    u = -u-1;
               fw = False;    }
          placements_by_read[ jaligns0[j].first ].push( u, jaligns0[j].third, fw );
          placements_by_unipath[u].push( jaligns0[j].first, jaligns0[j].third, fw );
               }

     // Generate an alternate version of jump alignments in which we use the filt
     // reads and truncate to get perfect aligments.  There has to be a better way
     // of doing this.  The reason we use the filt reads instead of the ec reads is
     // in repeat regions, some ec reads appear to be missing or conceivably 
     // miscorrected.

     cout << Date( ) << ": generating truncated filtered jump alignments" << endl;
     vecbasevector jbases_t;
     PairsManager jpairs_t;
     jpairs_t.Read( run_dir + "/jump_reads_filt.pairs" );
     jpairs_t.makeCache( );
     vec< vec< triple<int64_t,int,Bool> > > placements_by_unipath_t;
     vec< vec< triple<int,int,Bool> > > placements_by_read_t;
     String jump_t_fn = run_dir + "/jump_t";
     if ( IsRegularFile(jump_t_fn) && IsOlder( unibases2_fn, jump_t_fn ) )
     {    BinaryReader reader(jump_t_fn);
          reader.read( &jbases_t );
          reader.read( &placements_by_unipath_t );
          reader.read( &placements_by_read_t );    }
     else
     {    jbases_t.ReadAll( run_dir + "/jump_reads_filt.fastb" );
          int64_t nfjumps = jbases_t.size( );
          const int delta = 10;
          const int min_read_size = 40;
          int max_read = 0;
          for ( int64_t id = 0; id < nfjumps; id++ )
               max_read = Max( max_read, jbases_t[id].isize( ) );
          int start_size = (max_read/delta) * delta;
          placements_by_unipath_t.resize(nuni2);
          placements_by_read_t.resize(nfjumps);
          for ( int rs = start_size; rs >= min_read_size; rs -= delta )
          {    for ( int64_t id = 0; id < nfjumps; id++ )
               {    if ( placements_by_read_t[id].nonempty( ) ) continue;
                    basevector& b = jbases_t[id];
                    if ( b.isize( ) >= rs ) 
                         b.SetToSubOf( b, b.isize( ) - rs, rs );    }
               vec<basevector> jbases_sorted_t(nfjumps);
               vec<int64_t> jbases_sorted_id_t(nfjumps);
               for ( int64_t id = 0; id < nfjumps; id++ )
               {    if ( placements_by_read_t[id].empty( ) )
                    {    jbases_sorted_t[id] = jbases_t[id];
                         jbases_sorted_id_t[id] = id;    }    }
               ParallelSortSync( jbases_sorted_t, jbases_sorted_id_t );
               vec< triple<int64_t,int,int> > jaligns0_t;
               #pragma omp parallel for
               for ( int m = 0; m < nuni2; m++ )
               {    vec< triple<int64_t,int,int> > jaligns0m_t;
                    for ( int p = 0; p <= unibases2[m].isize( ) - rs; p++ )
                    {    basevector b( unibases2[m], p, rs );
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    if ( pass == 2 ) b.ReverseComplement( );
                              int64_t low = LowerBound( jbases_sorted_t, b );
                              int64_t high = UpperBound( jbases_sorted_t, b );
                              for ( int64_t l = low; l < high; l++ )
                              {    int64_t id = jbases_sorted_id_t[l];
                                   int mx = m;
                                   if ( pass == 2 ) mx = -m-1;
                                   jaligns0m_t.push( id, mx, p );    }    }    }
                    #pragma omp critical
                    {    jaligns0_t.append(jaligns0m_t);    }    }
               ParallelSort(jaligns0_t);
               for ( int64_t j = 0; j < jaligns0_t.isize( ); j++ )
               {    int u = jaligns0_t[j].second;
                    Bool fw = True;
                    if ( u < 0 )
                    {    u = -u-1;
                         fw = False;    }
                    int64_t id = jaligns0_t[j].first;
                    placements_by_read_t[id].push( u, jaligns0_t[j].third, fw );
                    placements_by_unipath_t[u].push( 
                         id, jaligns0_t[j].third, fw );    }    }
          BinaryWriter writer(jump_t_fn);
          writer.write(jbases_t);
          writer.write(placements_by_unipath_t);
          writer.write(placements_by_read_t);    }
     if (DUMP_PAIRS)
     {    cout << "\nDUMP OF PAIRS\n\n";
          int64_t total_bases = 0, total_aligned_bases = 0;
          for ( int64_t pid = 0; pid < (int64_t) jpairs_t.nPairs( ); pid++ )
          {    int64_t id1 = jpairs_t.ID1(pid), id2 = jpairs_t.ID2(pid);
               total_bases += jbases_t[id1].isize( ) + jbases_t[id2].isize( );
               if ( placements_by_read_t[id1].nonempty( )
                    && placements_by_read_t[id2].nonempty( ) )
               {    total_aligned_bases 
                         += jbases_t[id1].isize( ) + jbases_t[id2].isize( );     }
               for ( int j1 = 0; j1 < placements_by_read_t[id1].isize( ); j1++ )
               for ( int j2 = 0; j2 < placements_by_read_t[id2].isize( ); j2++ )
               {    if ( !placements_by_read_t[id1][j1].third ) continue;
                    if ( placements_by_read_t[id2][j2].third ) continue;
                    int u1 = placements_by_read_t[id1][j1].first;
                    int start1 = placements_by_read_t[id1][j1].second;
                    int stop1 = start1 + jbases_t[id1].isize( );
                    int u2 = placements_by_read_t[id2][j2].first;
                    int start2 = placements_by_read_t[id2][j2].second;
                    int stop2 = start2 + jbases_t[id2].isize( );
                    cout << u1 << "." << start1 << "-" << stop1 << "    " << u2 
                         << "." << start2 << "-" << stop2 << "\n";    }    }
          cout << "\n";
          PRINT2( total_bases, total_aligned_bases );
          Scram();     }

     // Get quick and dirty upper bounds on jump distribution.  Bad.

     const int min_size = 10000;
     const int max_sep = 20000;
     vec<int> dists;
     for ( int64_t id1 = 0; id1 < (int64_t) jbases.size( ); id1++ )
     {    int64_t id2 = jpairs.getPartnerID(id1);
          for ( int j1 = 0; j1 < placements_by_read[id1].isize( ); j1++ )
          for ( int j2 = 0; j2 < placements_by_read[id2].isize( ); j2++ )
          {    int u = placements_by_read[id1][j1].first;
               if ( unibases2[u].isize( ) < min_size ) continue;
               if ( placements_by_read[id2][j2].first != u ) continue;
               int pos1 = placements_by_read[id1][j1].second;
               int pos2 = placements_by_read[id2][j2].second;
               Bool fw1 = placements_by_read[id1][j1].third;
               Bool fw2 = placements_by_read[id2][j2].third;
               int len1 = jbases[id1].size( ), len2 = jbases[id2].size( );
               if ( !fw1 || fw2 ) continue;
               int Pos2 = pos2 + len2;
               int dist = Pos2 - pos1;
               if ( dist >= 0 && dist <= max_sep ) dists.push_back(dist);    }    }
     ReverseSort(dists);
     int distsq1 = 2000, distsq2 = 5000, distsq3 = 10000;
     size_t distsSz = dists.size();
     if ( distsSz )
     {    distsq1 = dists[ 3 * distsSz / 4 ];
          distsq2 = dists[ distsSz / 2 ];
          distsq3 = dists[ distsSz / 4 ];    }
     DPRINT4( distsq1, distsq2, distsq3, distsSz );
     Sort(dists);
     
     // Dump.

     if (DUMP)
     {
     const int max_far = 10000;
     for ( int u1 = 0; u1 < (int) unibases2.size( ); u1++ )
     {    cout << "\nunique lookaheads from " << u1 << "\n";
          vec<int> ahead;
          for ( int j = 0; j < placements_by_unipath[u1].isize( ); j++ )
          {    int64_t id1 = placements_by_unipath[u1][j].first;
               if ( placements_by_read[id1].size( ) != 2 ) continue;
               int pos1 = placements_by_unipath[u1][j].second;
               int dist1 = unibases2[u1].isize( ) - pos1 - jbases[id1].isize( );
               if ( dist1 > max_far ) continue;
               Bool fw1 = placements_by_unipath[u1][j].third;
               if ( !fw1 ) continue;
               int64_t id2 = jpairs.getPartnerID(id1);
               if ( placements_by_read[id2].size( ) != 2 ) continue;
               int i2 = 0;
               if ( placements_by_read[id2][i2].third ) i2 = 1;
               int u2 = placements_by_read[id2][i2].first;
               if ( u2 == u1 ) continue;
               int pos2 = placements_by_read[id2][i2].second;
               int dist2 = pos2;
               if ( dist2 > max_far ) continue;
               Bool fw2 = placements_by_read[id2][i2].third;
               if (fw2) continue;
               if (DUMP_VERBOSE) PRINT6( u1, u2, dist1, dist2, id1, id2 );
               ahead.push_back(u2);    }
          Sort(ahead);
          vec< pair<int,int> > aheadx;
          for ( int i = 0; i < ahead.isize( ); i++ )
          {    int j = ahead.NextDiff(i);
               aheadx.push( j-i, ahead[i] );   
               i = j - 1;    }
          ReverseSort(aheadx);
          for ( int i = 0; i < aheadx.isize( ); i++ )
          {    int u = aheadx[i].second;
               cout << u << " (mult=" << aheadx[i].first 
                    << ",len=" << unibases2[u].isize( ) - K2 + 1 << ")\n";    }    }
     return 0;
     }
               
     // Map large unibases1 to unibases2.

     vec< vec< pair<int,int> > > to_big;
     if ( K2 == 640 ) ToBig<640>( unibases1, unibases2, to_big );
     else
     {    cout << "Not implemented for K2 = " << K2 << endl;
          exit(1);    }
     if (SHOW_TRANSLATION)
     {    cout << "\nMap from unibases1 to unibases2\n\n";
          for ( int u1 = 0; u1 < nuni1; u1++ )
          {    int nkmers1 = unibases1[u1].isize( ) - K1 + 1;
               if ( to_big[u1].size( ) > 1 ) cout << u1 << " MULTIMAPS!\n";
               for ( int j = 0; j < to_big[u1].isize( ); j++ )
               {    int u2 = to_big[u1][j].first, pos2 = to_big[u1][j].second;
                    int cn = predicted_CN1[u1];
                    prob_t maxp = 0;
                    for ( uint z = 0; z < CN1[u1].size( ); z++ )
                         if ( CN1[u1][z].first == cn ) maxp = CN1[u1][z].second;
                    String prob = ToString( 100.0 * maxp, 1 );
                    cout << u1 << "[l=" << nkmers1 << ",rc=" << to_rc1[u1] << ",cn="
                         << cn << "(" << prob << "%)" << ",cn_raw=" 
                         << setiosflags(ios::fixed) << setprecision(2) << CN_raw[u1]
                         << resetiosflags(ios::fixed) << "] --> " << u2 << "[l=" 
                         << unibases2[u2].isize( ) - K2 + 1 << ",rc=" << to_rc2[u2] 
                         << "] at " << pos2 << "\n";    }    }
          cout << "\n";    }

     // Define bridgeable gaps.  This uses a hardcoded upper limit on the
     // library, which is bad.

     cout << Date( ) << ": defining bridgeable gaps" << endl;
     vec< vec<int> > bridge_to(nuni1), bridge_to_dist(nuni1);
     const int max_dist = 15000;
     #pragma omp parallel for
     for ( int u = 0; u < nuni1; u++ )
     {    vec< triple<int,int,Bool> > x;
          x.push( u, 0, False );
          while(1)
          {    Bool progress = False;
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( x[j].third ) continue;
                    int v = x[j].first, d = x[j].second;
                    if ( d > max_dist ) continue;
                    for ( int i = 0; i < nexts1[v].isize( ); i++ )
                    {    int vn = nexts1[v][i];
                         int dn = d + unibases1[vn].isize( ) - (K1-1);
                         Bool found = False;
                         for ( int r = 0; r < x.isize( ); r++ )
                         {    if ( x[r].first == vn )
                              {    found = True;
                                   if ( x[r].second > dn )
                                   {    x[r].second = dn;
                                        progress = True;
                                        x[r].third = False;    }
                                   break;    }    }
                         if ( !found ) 
                         {    x.push( vn, dn, dn >= max_dist );    
                              progress = True;    }    }    }
               if ( !progress ) break;    }
          for ( int j = 0; j < x.isize( ); j++ )
          {    bridge_to[u].push_back( x[j].first );
               bridge_to_dist[u].push_back( x[j].second 
                    - ( unibases1[ x[j].first ].isize( ) - (K1-1) ) );    }
          SortSync( bridge_to[u], bridge_to_dist[u] );    }

     // Load predicted gaps and use them to make a digraph.

     cout << Date( ) << ": loading predicted gaps" << endl;
     digraphE<linklet> G;
     {    vec< vec<int> > from(nuni1), to(nuni1);
          vec< vec<int> > from_edge_obj(nuni1), to_edge_obj(nuni1);
          vec<linklet> edges;
          String line;
          fast_ifstream in( 
               run_dir + "/" + HEAD1 + ".unibases.k" + KS1 + ".predicted_gaps.txt" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               int u1, u2, sep, dev, nlinks;
               istringstream iline( line.c_str( ) );
               iline >> u1 >> u2 >> sep >> dev >> nlinks;
               int nkmers1 = unibases1[u1].isize( ) - K1 + 1;
               int nkmers2 = unibases1[u2].isize( ) - K1 + 1;
               if ( nkmers1 < MIN_KMERS1 || nkmers2 < MIN_KMERS1 ) continue;
               if ( u2 == to_rc1[u1] ) continue; 
               if ( nlinks < MIN_LINKS_FOR_UNBRIDGEABLE )
               {    int p = BinPosition( bridge_to[u1], u2 );
                    const int max_devs = 6;
                    if ( p < 0 || bridge_to_dist[u1][p] > sep + max_devs * dev )
                    {    continue;    }    }
               from[u1].push_back(u2), to[u2].push_back(u1);
               from_edge_obj[u1].push_back( edges.size( ) );
               to_edge_obj[u2].push_back( edges.size( ) );
               edges.push( sep, dev, nlinks, 0 );    }
          for ( int u = 0; u < nuni1; u++ )
          {    SortSync( from[u], from_edge_obj[u] );
               SortSync( to[u], to_edge_obj[u] );    }
          G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

     // Parameters for identifying unipaths that don't have copy number one.

     const int max_link_ratio = 5;
     const int dev_mult = 3;

     // Print initial links.

     if (SHOW_INITIAL_LINKS)
     {    cout << "\nInitial unibase1 links\n\n";
          for ( int u1 = 0; u1 < nuni1; u1++ )
          for ( int j = 0; j < G.From(u1).isize( ); j++ )
          {    int u2 = G.From(u1)[j];
               int nkmers1 = unibases1[u1].isize( ) - K1 + 1;
               int nkmers2 = unibases1[u2].isize( ) - K1 + 1;
               const linklet& l = G.EdgeObjectByIndexFrom( u1, j );
               cout << "u1 = " << u1 << "[l=" << nkmers1
                    /* << ",cn=" << predicted_CN1[u1] */
                    << "], u2 = " << u2
                    << "[l=" << nkmers2
                    /* << ",cn=" << predicted_CN1[u2] */
                    << "], ";
               int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
               cout << "sep = " << sep << " +/- " << dev << ", nlinks = "
                    << nlinks << "\n";    }
          cout << "\n";    }

     // Derive initial links2.

     digraphE<linklet> G2I;
     RaiseLinks( K2, unibases1, unibases2, to_rc2, to_big, G, G2I );

     // Compute raw copy number for unibases2.

     vec<double> raw2(nuni2, -1.0);
     {    vec< vec<int> > from_big(nuni2);
          for ( int u = 0; u < nuni1; u++ )
          {    for ( int j = 0; j < to_big[u].isize( ); j++ )
                    from_big[ to_big[u][j].first ].push_back(u);    }
          for ( int u = 0; u < nuni2; u++ )
          {    int total_len = 0;
               double weighted_raw = 0.0;
               for ( int j = 0; j < from_big[u].isize( ); j++ )
               {    int v = from_big[u][j];
                    double raw = CN_raw[v];
                    int len = unibases1[v].isize( ) - K1 + 1;
                    total_len += len;
                    weighted_raw += raw * double(len);    }
               if ( from_big[u].size( ) > 0 )
               raw2[u] = weighted_raw / double(total_len);    }    }

     // Remove weak links.

     cout << Date( ) << ": removing weak links" << endl;
     RemoveWeakLinks2( K2, unibases2, to_rc2, raw2, G2I, nexts2x, max_link_ratio,
          dev_mult, DELETE_VERBOSE, MAX_DEV_DIFF );

     // Print links.

     cout << Date( ) << ": printing links" << endl;
     if (SHOW_INITIAL_LINKS2)
     {    cout << "\nInitial unibase2 links\n" << endl;
          for ( int u1 = 0; u1 < nuni2; u1++ )
          for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
          {    int u2 = G2I.From(u1)[j];
               int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
               int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
               const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
               cout << "u1 = " << u1 << "[l=" << nkmers1 << "], u2 = " << u2
                    << "[l=" << nkmers2 << "], ";
               int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
               cout << "sep = " << sep << " +/- " << dev << ", nlinks = " 
                    << nlinks << "\n";    }
          cout << "\n";    }

     // Print unibase graphs.

     if (PRINT_GRAPH1)
     {    cout << "\nunibases1 graph:\n\n";
          for ( int u = 0; u < nuni1; u++ )
          {    cout << u << "[l=" << unibases1[u].isize( ) - K1 + 1
                    << ",cn=" << CN_raw[u] << ",rc=" << to_rc1[u] << "] -->";
               for ( int j = 0; j < nexts1[u].isize( ); j++ )
                    cout << " " << nexts1[u][j];
               cout << "\n";    }
          cout << "\n";    }
     if (PRINT_GRAPH2)
     {    cout << "\nunibases2 graph:\n" << endl;
          for ( int u = 0; u < nuni2; u++ )
          {    cout << u << "[l=" << unibases2[u].isize( ) - K2 + 1
                    << ",rc=" << to_rc2[u] << "] -->";
               for ( int j = 0; j < nexts2[u].isize( ); j++ )
                    cout << " " << nexts2[u][j];
               cout << "\n";    }
          cout << "\n";    }

     // Identify certain unipaths that likely repeats and which are 
     // 'linked around'.  These will then be excluded from subsequent building.
     // (Started, turned off.)

     /*
     for ( int u = 0; u < nuni2; u++ )
     {    int n1 = G2I.To(u).size( ), n2 = G2I.From(u).size( );
          equiv_rel e1(n1), e2(n2);
          for ( int xa = 0; xa < n1; xa++ )
          {    int v = G2I.To(u)[xa];
               for ( int j = 0; j < G2I.From(v).isize( ); j++ )
               {    int xb = Position( G2I.To(u), G2I.From(v)[j] );
                    if ( xb < 0 ) continue;
                    e1.Join( xa, xb );    }    }
          for ( int xa = 0; xa < n2; xa++ )
          {    int v = G2I.From(u)[xa];
               for ( int j = 0; j < G2I.From(v).isize( ); j++ )
               {    int xb = Position( G2I.From(u), G2I.From(v)[j] );
                    if ( xb < 0 ) continue;
                    e2.Join( xa, xb );    }    }
     */

     vec< pair<int,int> > uni2_len;
     for ( int u = 0; u < nuni2; u++ )
          uni2_len.push( unibases2[u].size( ), u );
     ReverseSort(uni2_len);

     vec<String> reports(nuni2);
     vec< vec<int> > rights(nuni2), lefts(nuni2);
     vec< vec<int> > unexts(nuni2);
     cout << Date( ) << ": walking right" << endl;
     for ( int ju = 0; ju < nuni2; ju++ )
     {    int u = uni2_len[ju].second;
          if ( unibases2[u].isize( ) - K2 + 1 < min_kmers2 )
          {    continue;    }
          ostringstream out;
          out << "\nwalking right from " << u 
               << "[" << unibases2[u].isize( ) - K2 + 1 << ",cn=" 
               << setiosflags(ios::fixed) << setprecision(2) << raw2[u]
               << resetiosflags(ios::fixed) << "]";
          if (VALIDATE)
          {    out << " (";
               for ( int j = 0; j < Ulocs2[u].isize( ); j++ )
                    Print( out, Ulocs2[u][j], Ulocs2[u].size( ) );
               out << " )";    }
          out << "\n";
          Advance( u, unibases2, to_rc2, raw2, K2, G2I, 
               max_link_ratio, out, unexts[u] );
          if ( unexts[u].solo( ) )
          {    int unext = unexts[u][0];
               rights[u].push_back(unext);
               lefts[unext].push_back(u);    
               Bool sym = True;
               if (sym)
               {    lefts[ to_rc2[u] ].push_back( to_rc2[unext] );
                    rights[ to_rc2[unext] ].push_back( to_rc2[u] );    }    }
          reports[u] = out.str( );    }
     if ( VERBOSITY >= 1 )
     {    cout << Date( ) << ": print reports" << endl;
          for ( int u = 0; u < nuni2; u++ )
               if ( reports[u].size( ) > 0 ) cout << reports[u];    }

     // Now, very judiciously, bring back some of the cases where Advance returned
     // more than one solution.

     /*
     cout << "\nrehabilitating:";
     vec<int> use2;
     for ( int u = 0; u < nuni2; u++ )
     {    if ( unexts[u].size( ) <= 1 ) continue;
          if ( lefts[u].empty( ) && rights[u].empty( ) ) continue;
          Bool indirect = True;
          for ( int j = 0; j < unexts[u].isize( ); j++ )
          {    int v = unexts[u][j];
               Bool hit = False;
               for ( int l = 0; l < G2I.To(u).isize( ); l++ )
               {    int w = G2I.To(u)[l];
                    for ( int m = 0; m < G2I.From(w).isize( ); m++ )
                    {    int x = G2I.From(w)[m];
                         if ( x == v ) hit = True;    }    }
               if ( !hit ) indirect = False;    }
          if (indirect) continue;
          cout << " " << u;
          use2.push_back(u);    }
     cout << "\n";
     for ( int i = 0; i < use2.isize( ); i++ )
     {    int u = use2[i];
          for ( int j = 0; j < unexts[u].isize( ); j++ )
          {    int unext = unexts[u][j];
               rights[u].push_back(unext);
               lefts[unext].push_back(u);    
               Bool sym = True;
               if (sym)
               {    lefts[ to_rc2[u] ].push_back( to_rc2[unext] );
                    rights[ to_rc2[unext] ].push_back( to_rc2[u] );    }    }    }
     */

     // Define vertices to use.  Don't use vertices in reverse complement
     // components.

     vec<int> verts;
     equiv_rel e(nuni2);
     for ( int u1 = 0; u1 < nuni2; u1++ )
     {    for ( int j = 0; j < rights[u1].isize( ); j++ )
               e.Join( u1, rights[u1][j] );    }
     for ( int u = 0; u < nuni2; u++ )
     {    if ( e.ClassId( to_rc2[u] ) < e.ClassId(u) ) continue;
          if ( rights[u].nonempty( ) || lefts[u].nonempty( )
               || unibases2[u].isize( ) >= 2000 )
          {    verts.push_back(u);    }    }

     // Make graph.

     cout << Date( ) << ": making graph" << endl;
     int nv = verts.size( );
     vec< vec<int> > from(nv), to(nv);
     for ( int x1 = 0; x1 < nv; x1++ )
     {    int u1 = verts[x1];
          for ( int j = 0; j < rights[u1].isize( ); j++ )
          {    int u2 = rights[u1][j];
               int x2 = Position( verts, u2 );
               if ( x2 < 0 ) continue;
               from[x1].push_back(x2), to[x2].push_back(x1);    }    }
     for ( int x = 0; x < nv; x++ )
     {    UniqueSort(from[x]), UniqueSort(to[x]);    }
     digraph H( from, to );

     // Linearize squares   b ----> d
     //                     ^       ^    ====>  a --> b --> c --> d
     //                     |       |
     //                     a ----> c
     // where possible.

     for ( int a = 0; a < H.N( ); a++ )
     {    if ( H.From(a).size( ) != 2 ) continue;
          for ( int pass = 1; pass <= 2; pass++ )
          {    int b = H.From(a)[0], c = H.From(a)[1];
               if ( pass == 2 ) swap( b, c );
               if ( !H.From(b).solo( ) || !H.To(b).solo( ) ) continue;
               if ( !H.From(c).solo( ) || !H.To(c).solo( ) ) continue;
               if ( H.From(b)[0] != H.From(c)[0] ) continue;
               int d = H.From(b)[0];
               if ( H.To(d).size( ) != 2 ) continue;
               if ( a == b || a == c || a == d || b == c || b == d || c == d ) 
                    continue;
               if ( !Member( G2I.From( verts[b] ), verts[c] ) ) continue;
               H.FromMutable(b)[0] = c;
               H.ToMutable(c)[0] = b;
               H.ToMutable(d).resize(1);
               H.ToMutable(d)[0] = c;
               H.FromMutable(a).resize(1);
               H.FromMutable(a)[0] = b;    }    }

     // Linearize triangles   b        
     //                       ^          ====>  a --> b --> c
     //                       |        
     //                       a ----> c
     // where possible.

     for ( int count = 0; count < 2; count++ )
     {
     for ( int xpass = 1; xpass <= 2; xpass++ )
     {    H.Reverse( );
          G2I.Reverse( );
          for ( int a = 0; a < H.N( ); a++ )
          {    if ( H.From(a).size( ) != 2 ) continue;
               for ( int pass = 1; pass <= 2; pass++ )
               {    int b = H.From(a)[0], c = H.From(a)[1];
                    if ( pass == 2 ) swap( b, c );
                    if ( !H.To(b).solo( ) ) continue;
                    if ( !Member( G2I.From( verts[b] ), verts[c] ) ) continue;
                    int pa = Position( H.To(c), a ), pb = Position( H.To(c), b );
                    H.FromMutable(a).resize(1);
                    H.FromMutable(a)[0] = b;
                    if ( pb < 0 )
                    {    H.FromMutable(b).push_back(c);
                         H.ToMutable(c)[pa] = b;
                         Sort( H.FromMutable(b) );
                         Sort( H.ToMutable(c) );    }
                    else 
                    {    H.ToMutable(c).erase( 
                              H.ToMutable(c).begin( ) + pa );    }    
                    break;    }    }    }
     }

     // Build a snark.

     vec<uniseq> seq;
     for ( int i = 0; i < verts.isize( ); i++ )
     {    vec<int> u(1), overlap;
          u[0] = verts[i];
          seq.push( u, overlap );    }
     vec< vec<int> > from_edge_obj( verts.size( ) ), to_edge_obj( verts.size( ) );
     vec<gapster> edges;
     for ( int x1 = 0; x1 < verts.isize( ); x1++ )
     {    for ( int j = 0; j < H.From(x1).isize( ); j++ )
          {    int x2 = H.From(x1)[j];
               int u1 = verts[x1], u2 = verts[x2];
               int lx = BinPosition( G2I.From(u1), u2 );
               ForceAssertGe( lx, 0 );
               const linklet& l = G2I.EdgeObjectByIndexFrom( u1, lx );
               from_edge_obj[x1].push( edges.size( ) );
               edges.push( l.sep, l.dev );    }    }
     for ( int x2 = 0; x2 < verts.isize( ); x2++ )
     {    for ( int j = 0; j < H.To(x2).isize( ); j++ )
          {    int x1 = H.To(x2)[j];
               int k = BinPosition( H.From(x1), x2 );
               to_edge_obj[x2].push( from_edge_obj[x1][k] );    }    }
     digraphE<gapster> gapsters( 
          H.From( ), H.To( ), edges, to_edge_obj, from_edge_obj );
     snark S( gapsters, seq );
     S.SetUnibases(unibases2);
     S.SetToRc(to_rc2);
     uniseq dummy;
     dummy.SetUnibases(unibases2);

     // Handle simple inverted repeats.  This might better be done later, after 
     // the graph is simplified.

     S.HandleSimpleInvertedRepeats( );

     // Look for closures.

     cout << Date( ) << ": looking for unique paths" << endl;
     for ( int x1 = 0; x1 < S.VertN( ); x1++ )
     {    int u1 = S.Vert(x1).U(0);
          for ( int j = 0; j < S.From(x1).isize( ); j++ )
          {    int x2 = S.From(x1)[j];
               int u2 = S.Vert(x2).U(0);
               vec< vec< pair<int,int> > > walks1;

               vec<int> use;
               const int flank = 3000;
               if ( S.Vert(x1).Len( ) >= flank && S.Vert(x2).Len( ) >= flank )
               {    vec< triple<int,int,int> > between;
                    FindPartnersInGap( x1, x2, flank, S, jbases, jpairs,
                         placements_by_read, placements_by_unipath, between );
                    if (SHOW_PARTNERS_IN_GAP)
                    {    cout << "\npartners in gap:\n";
                         for ( int j = 0; j < between.isize( ); j++ )
                         {    int u = between[j].first;
                              cout << u << "." << between[j].second
                                   << "-" << between[j].third << " (of " 
                                   << unibases2[u].isize( ) << ")\n";    }    }
                    for ( int j = 0; j < between.isize( ); j++ )
                         use.push_back( between[j].first );
                    UniqueSort(use);    }

               int bad, level = 1;
               int sep = S.G( ).EdgeObjectByIndexFrom( x1, j ).Sep( );
               int dev = S.G( ).EdgeObjectByIndexFrom( x1, j ).Dev( );
               GetWalks( u1, u2, sep, dev, unibases2, K2, to_rc2, nexts2y, use,
                    walks1, bad, MAX_DEV_DIFF );

               if ( bad == 0 && walks1.empty( ) )
               {    level = 2;
                    GetWalks( u1, u2, sep, dev, unibases2, K2, to_rc2, nexts2x, use,
                         walks1, bad, MAX_DEV_DIFF );    }

               if ( bad == 0 && walks1.nonempty( ) )
               {    vec<uniseq> closures;
                    for ( int l = 0; l < walks1.isize( ); l++ )
                    {    vec<int> u, over;
                         for ( int m = 0; m < walks1[l].isize( ); m++ )
                         {    u.push_back( walks1[l][m].first );
                              if ( m > 0 ) 
                                   over.push_back( walks1[l][m].second );    }
                         closures.push( u, over );    }
                    S.CloseGap( x1, x2, closures );    }

               if ( VERBOSITY >= 1 )
               {    cout << "\ninitial walks from " << u1 << " to " << u2 
                         << ", level = " << level << ", bad = " << bad 
                         << ", nwalks = " << walks1.size( ) << "\n";
                    for ( int j = 0; j < Min( MAX_WALKS_TO_PRINT, walks1.isize( ) ); 
                         j++ )
                    {    cout << "[" << j << "] ";
                         const vec< pair<int,int> >& w = walks1[j];
                         vec<int> u, over;
                         for ( int m = 0; m < w.isize( ); m++ )
                         {    u.push_back( w[m].first );
                              if ( m > 0 ) over.push_back( w[m].second );    }
                         uniseq s( u, over );
                         s.Print( cout, K2 );
                              cout << "\n";    }
                    if ( walks1.isize( ) > MAX_WALKS_TO_PRINT ) 
                         cout << "...\n";    }    }    }
     S.SwallowSimpleGaps( );
     S.BringOutTheDead( );

     // Prefilter, looking for the absence of links.

     cout << Date( ) << ": prefilter, looking for absence of links" << endl;
     for ( int e = 0; e < S.EdgeN( ); e++ )
     {    gapster& g = S.Gmutable( ).EdgeObjectMutable(e);
          vec<Bool> to_remove( g.ClosureCount( ), False );
          #pragma omp parallel for
          for ( int j = 0; j < g.ClosureCount( ); j++ )
          {    const uniseq& u = g.Closure(j);
               for ( int j1 = 0; j1 < u.N( ); j1++ )
               {    int u1 = u.U(j1); 
                    int n1 = unibases2[u1].isize( ) - (K2-1);
                    const int min_kmers = 1000;
                    if ( n1 < min_kmers ) continue;
                    for ( int j2 = j1 + 1; j2 < u.N( ); j2++ )
                    {    int u2 = u.U(j2), dist = -K2+1;
                         if ( u2 == u1 || u2 == to_rc2[u1] ) continue;
                         int n2 = unibases2[u2].isize( ) - (K2-1);
                         if ( n2 < min_kmers ) continue;
                         for ( int l = j1 + 1; l < j2; l++ )
                              dist += unibases2[ u.U(l) ].isize( ) - (K2-1);
                         if ( dist <= 0 ) continue;
                         if ( dist > 10000 ) break;
                         vec<int> links( 2, 0 );
                         vec<double> predicted( 2, 0.0 );
                         // These two passes don't really make sense, not
                         // sure what I was thinking.
                         for ( int pass = 0; pass < 2; pass++ )
                         {    int v1 = to_rc2[u2], v2 = to_rc2[u1];
                              u1 = v1, u2 = v2; 
                              for ( int l = 0; 
                                   l < placements_by_unipath[u1].isize( ); l++ )
                              {    if ( !placements_by_unipath[u1][l].third ) 
                                        continue;
                                   int64_t id1 = placements_by_unipath[u1][l].first;
                                   int64_t id2 = jpairs.getPartnerID(id1);
                                   if ( placements_by_read[id2].empty( ) ) continue;
                                   int pos1 = placements_by_unipath[u1][l].second;
                                   int low = unibases2[u1].isize( ) 
                                        - pos1 - jbases[id1].isize( ) + dist;
                                   int high = low + unibases2[u2].isize( );
                                   int64_t start = LowerBound( dists, low );
                                   int64_t stop = UpperBound( dists, high );
                                   predicted[pass] += double(stop-start) 
                                        / double( dists.size( ) );
                                   for ( int m = 0; 
                                        m < placements_by_read[id2].isize( ); m++ )
                                   {    if ( placements_by_read[id2][m].third ) 
                                             continue;
                                        if ( placements_by_read[id2][m].first != u2 )
                                             continue;
                                        links[pass]++;    }    }    }
                         const double min_pred = 20.0;
                         const double div = 10.0;
                         if ( predicted[0] < min_pred || predicted[1] < min_pred ) 
                              continue;
                         if ( links[0] >= predicted[0]/div || links[1] 
                              >= predicted[1]/div )
                         {    continue;    }
                         to_remove[j] = True;
                         if ( VERBOSITY >= 1 )
                         {    cout << "\nnot enough links, screening out\n";
                              int pred1 = int(round(predicted[0]));
                              int pred2 = int(round(predicted[1]));
                              int links1 = links[0], links2 = links[1];
                              PRINT9( u1, u2, n1, n2, dist, links1, pred1, 
                                   links2, pred2 );    }    }    }    }
          g.RemoveSomeClosures(to_remove);    }

     // Attempt to remove fw/rc duplication.  Note that filtering walks using jumps
     // can introduce asymmetry (not sure how), so either the removal needs to be
     // here or we need to figure out how the asymmetry gets in.

     if (DELETE_COMPLEMENTS) S.DeleteComplements( );

     // Filter walks using jumps.

     cout << Date( ) << ": filter walks using jumps" << endl;
     FilterByJumps( S, K2, jbases, jpairs, placements_by_read, 
          placements_by_unipath, VERBOSITY, FILTER_VERBOSITY, FILTER_VERBOSE_U1, 
          FILTER_SHOW_BASES_AFTER );
     S.RemoveSubsumedStuff( );
     S.SwallowSimpleGaps( );
     S.BringOutTheDead( );

     // Look again at walks.  Find jump pair placements on each.  We ignore 
     // placements for which both reads are placed on an end unipath. 

     for ( int xv = 0; xv < S.VertN( ); xv++ )
     for ( int j = 0; j < S.G( ).From(xv).isize( ); j++ )
     {    int y = S.G( ).From(xv)[j];
          gapster& g = S.Gmutable( ).EdgeObjectByIndexFromMutable( xv, j );
          int e = S.G( ).EdgeObjectIndexByIndexFrom( xv, j );
          if ( g.Open( ) ) continue;
          int v1 = g.Closure(0).U( ).front( ), v2 = g.Closure(0).U( ).back( );
          vec<int> uncov_count( g.ClosureCount( ) );
          for ( int i = 0; i < g.ClosureCount( ); i++ )
          {    const uniseq& x = g.Closure(i);
               int xkmers = 0;
               for ( int l = 1; l < x.N( ) - 1; l++ )
                    xkmers += unibases2[ x.U(l) ].isize( ) - (K2-1);
               vec< pair<int,int> > phits;
               vec< pair<int64_t,int64_t> > pids;
               for ( int m1 = 0; m1 < x.N( ); m1++ )
               {    int u1 = x.U(m1);
                    for ( int l1 = 0; 
                         l1 < placements_by_unipath_t[u1].isize( ); l1++ )
                    {    if ( !placements_by_unipath_t[u1][l1].third ) continue;
                         int64_t id1 = placements_by_unipath_t[u1][l1].first;
                         int pos1 = placements_by_unipath_t[u1][l1].second;
                         int64_t id2 = jpairs_t.getPartnerID(id1);
                         for ( int l2 = 0; l2 < placements_by_read_t[id2].isize( );
                              l2++ )
                         {    if ( placements_by_read_t[id2][l2].third ) continue;
                              int u2 = placements_by_read_t[id2][l2].first;
                              int pos2 = placements_by_read_t[id2][l2].second;
                              int Pos2 = pos2 + jbases_t[id2].isize( );
                              for ( int m2 = m1; m2 < x.N( ); m2++ )
                              {    if ( x.U(m2) != u2 ) continue;
                                   if ( m1 == m2 && ( m1 == 0 || m1 == x.N( ) - 1 ) )
                                        continue;
                                   int pos1_alt = pos1, Pos2_alt = Pos2;
                                   for ( int r = 0; r < m1; r++ )
                                   {    pos1_alt += 
                                             unibases2[ x.U(r) ].isize( )-(K2-1);   }
                                   for ( int r = 0; r < m2; r++ )
                                   {    Pos2_alt += 
                                             unibases2[ x.U(r) ].isize( )-(K2-1);   }
                                   int len = Pos2_alt - pos1_alt;
                                   pos1_alt -= unibases2[v1].isize( );
                                   Pos2_alt -= unibases2[v1].isize( );
                                   if ( len < distsq1 || len > distsq3 ) continue;
                                   // int64_t pid = jpairs_t.getPairID(id1); // XXXX
                                   // PRINT(pid); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                                   phits.push( pos1_alt, Pos2_alt );    
                                   pids.push( id1, id2 );    }    }    }    }
               UniqueSortSync( phits, pids );
               if ( VERBOSITY >= 1 )
               {    cout << "\nlooking for uncovered bases on edge " << e << endl;
                    PRINT5( xv, j, i, v1, v2 );    }
               int n1 = unibases2[v1].size( ), n2 = unibases2[v2].size( );
               int total = xkmers + n1 + n2 - (K2-1);
               vec<ho_interval> rcov, runcov;
               rcov.push(0, n1), rcov.push(total-n2, total);
               const int autocov = 4000;
               rcov.push( 0, Min( autocov, total ) );
               rcov.push( Max( 0, total - autocov ), total );
               for ( int m = 0; m < phits.isize( ); m++ )
               {    pair<int,int> h = phits[m];
                    int64_t pid = h.first;
                    int64_t id1 = pids[m].first, id2 = pids[m].second;
                    int rl1 = jbases_t[id1].size( ), rl2 = jbases_t[id2].size( );
                    if ( n1 + h.first >= 0 && n1 + h.first + rl1 <= total )
                         rcov.push( n1 + h.first, n1 + h.first + rl1 );
                    if ( n1 + h.second - rl2 >= 0 && n1 + h.second <= total )
                         rcov.push( n1 + h.second - rl2, n1 + h.second );     }
               Uncovered( total, rcov, runcov );
               uncov_count[i] = Sum(runcov);
               if ( VERBOSITY >= 2 )
               {    cout << "uncovered intervals:\n";
                    for ( int m = 0; m < runcov.isize( ); m++ )
                         cout << runcov[m] << "\n";
                    cout << "uncovered bases = " << uncov_count[i] << endl;    }    }
          int min_uncov = uncov_count.size() > 0 ? Min(uncov_count) : 0;
          const int uncov_mult = 2.0;
          vec<Bool> to_delete( g.ClosureCount( ), False );
          if ( VERBOSITY >= 1 ) cout << "\ncoverage summary for edge " << e << endl;
          for ( int i = 0; i < g.ClosureCount( ); i++ )
          {    Bool bad = ( uncov_count[i] > double(min_uncov) * uncov_mult );
               if ( VERBOSITY >= 1 ) 
                    PRINT4( i, uncov_count[i], min_uncov, int(bad) );
               if (bad) to_delete[i] = True;    }
          g.RemoveSomeClosures(to_delete);     }
     S.RemoveSubsumedStuff( );
     S.SwallowSimpleGaps( );
     S.BringOutTheDead( );

     // Look out from terminal vertices to see what we can see.  First find vertices
     // that are part of big stuff.

     const int min_mass = 10000;
     vec<Bool> big1( S.VertN( ), False );
     vec<Bool> big2( S.Unibases( ).size( ), False );
     vec< vec<int> > comp;
     S.G( ).Components(comp);
     for ( int c = 0; c < comp.isize( ); c++ )
     {    int mass = 0;
          for ( int j = 0; j < comp[c].isize( ); j++ )
               mass += S.Vert( comp[c][j] ).Len( );
          if ( mass >= min_mass )
          {    for ( int j = 0; j < comp[c].isize( ); j++ )
               {    int x = comp[c][j];
                    big1[x] = True;
                    big2[ S.Vert(x).U(0) ] = True;    }    }    }

     // Now go through the vertices.

     for ( int x = 0; x < S.VertN( ); x++ )
     {    if ( S.From(x).nonempty( ) ) continue;

          // Vertex must be connected to a lot.

          if ( !big1[x] ) continue;

          // Establish the dominant copy number of x.

          int max_size = 0, best = -1;
          for ( int j = 0; j < S.Vert(x).N( ); j++ )
          {    int u = S.Vert(x).U(j);
               if ( S.Unibase(u).isize( ) > max_size )
               {    max_size = S.Unibase(u).size( );
                    best = j;    }    }
          double cndom = raw2[ S.Vert(x).U(best) ];

          // Go through jump placements on x.

          double max_big = 1.5;
          const int max_dist_from_end = 10000;
          cout << "\n" << Date() << ": looking forward from " << S.Vert(x).U( ).back( ) << endl;
          vec<hit_thingee> hits;
          vec<Bool> hit( unibases2.size( ), False );
	  #pragma omp parallel for
          for ( int j = 0; j < S.Vert(x).N( ); j++ )
	  {    int u1 = S.Vert(x).U(j);
               if ( raw2[u1] < 0 || raw2[u1] > max_big * cndom ) continue;
               for ( int l = 0; l < placements_by_unipath[u1].isize( ); l++ )
               {    int64_t id1 = placements_by_unipath[u1][l].first;
                    int pos1 = placements_by_unipath[u1][l].second;
                    Bool unique1 
                         = cn2[u1][pos1] >= 0 && cn2[u1][pos1] <= max_big * cndom;
                    Bool fw1 = placements_by_unipath[u1][l].third;
                    if ( !fw1 ) continue;
                    int dist_from_end1 = unibases2[u1].isize( ) - pos1;
                    for ( int k = j + 1; k < S.Vert(x).N( ); k++ )
                    {    dist_from_end1 += S.Unibase( S.Vert(x).U(k) ).isize( )
                              - ( k < S.Vert(x).Over( ).isize( ) 
                                   ? S.Vert(x).Over(k) : 0 );    }
                    if ( dist_from_end1 > max_dist_from_end ) continue;
                    int64_t id2 = jpairs.getPartnerID(id1);
                    for ( int m = 0; m < placements_by_read[id2].isize( ); m++ )
                    {    int u2 = placements_by_read[id2][m].first;
                         int pos2 = placements_by_read[id2][m].second;
                         Bool fw2 = placements_by_read[id2][m].third;
                         if (fw2) continue;
			 #pragma omp critical
			 { hits.push( u2, pos2, id2, dist_from_end1, unique1 );
			   hit[u2] = True; }
                         /*
                         cout << "- see " << u2 << "." << pos2 
                              << ", id2 = " << id2 << "\n";    
                         */
                              }    }    }

          // Find paths extending u.

          int u = S.Vert(x).U( ).back( );
          vec< vec<int> > paths;
          vec<int> up;
          up.push_back(u);
          paths.push_back(up);
          const int max_paths = 10000;
          const int max_kmers = 10000;
          Bool more = True;
          while(1)
          {    if ( paths.isize( ) >= max_paths ) break;
               if ( !more ) break;
               more = False;
               vec< vec<int> > jpaths;
               for ( int j = 0; j < paths.isize( ); j++ )
               {    const vec<int>& p = paths[j];
                    int v = p.back( );
                    Bool extended = False;
                    int kmers = 0;
                    for ( int m = 1; m < p.isize( ); m++ )
                         kmers += S.Unibase( p[m] ).isize( ) - K2 + 1;
                    if ( kmers <= max_kmers )
                    {    for ( int l = 0; l < nexts2y[v].isize( ); l++ )
                         {    int w = nexts2y[v][l].first;
     
                              // We require that w is hit by a jump partner.
                              // This may be too stringent.
     
                              if ( !hit[w] ) continue;
     
                              // We do not allow a vertex to appear more than a
                              // an arbitrarily fixed number of times.
     
                              const int max_mult = 5;
                              int wcount = 0;
                              for ( int m = 0; m < p.isize( ); m++ )
                                   if ( p[m] == w ) wcount++;
                              if ( wcount == max_mult ) continue;
     
                              vec<int> pp(p);
                              pp.push_back(w);
                              extended = True;
                              jpaths.push_back(pp);    
                              more = True;    }    }
                    if ( !extended ) jpaths.push_back(p);    }
               paths = jpaths;    }
          if ( paths.isize( ) >= max_paths ) continue;

          // Define coverage and truncate paths that have too large a coverage hole.

          for ( int i = 0; i < paths.isize( ); i++ )
          {    vec<int>& p = paths[i];
               vec<int> h;
               int pos = 0;
               for ( int k = 0; k < p.isize( ); k++ )
               {    int v = p[k];
                    for ( int j = 0; j < hits.isize( ); j++ )
                    {    if ( hits[j].u2 != v ) continue;
                         h.push_back( pos + hits[j].pos2 );    }
                    pos += unibases2[v].isize( ) - K2 + 1;    }
               UniqueSort(h);
               const int max_hole = 500;
               for ( int j = 0; j < h.isize( ) - 1; j++ )
               {    if ( h[j] < unibases2[ p[0] ].isize( ) - K2 + 1 ) continue;
                    if ( h[j+1] - h[j] > max_hole )
                    {    int pos = 0;
                         for ( int k = 0; k < p.isize( ); k++ )
                         {    int v = p[k];
                              pos += unibases2[v].isize( ) - K2 + 1;
                              if ( pos > h[j] )
                              {    p.resize(k+2);
                                   /*
                                   cout << "would resize path " << i 
                                        << " to " << k+2 << "\n";
                                   PRINT3( j, h[j], h[j+1] );
                                   */
                                   break;    }    }
                         break;    }    }    }

          // Find acceptable read placements.  Note temporary hardcoded distance
          // limit.

          const int max_dist = 10000;
          const int unique_weight = 20;
          vec< vec< pair<int64_t,int> > > hids( paths.size( ) );
          for ( int i = 0; i < paths.isize( ); i++ )
          {    const vec<int>& p = paths[i];
               int ppos = 0;
               for ( int k = 0; k < p.isize( ); k++ )
               {    int v = p[k];
                    for ( int j = 0; j < hits.isize( ); j++ )
                    {    if ( hits[j].u2 != v ) continue;
                         int64_t id2 = hits[j].id2;
                         int pos2 = ppos + hits[j].pos2;
                         int dist = hits[j].dist_from_end1 + pos2 
                              + jbases[id2].isize( ) - S.Unibase( p[0] ).isize( );
                         if ( dist < 0 || dist > max_dist ) continue;
                         hids[i].push( id2, 
                              hits[j].unique1 ? unique_weight : 1 );    }
                    ppos += unibases2[v].isize( ) - K2 + 1;    }  
               UniqueSort( hids[i] );    }

          // Delete inferior paths.

          vec<Bool> to_delete( paths.size( ), False );
          for ( int j1 = 0; j1 < paths.isize( ); j1++ )
          {    if ( to_delete[j1] ) continue;
               for ( int j2 = 0; j2 < paths.isize( ); j2++ )
               {    if ( to_delete[j2] ) continue;
                    int better1 = 0, better2 = 0;
                    for ( int l = 0; l < hids[j1].isize( ); l++ )
                    {    if ( !BinMember( hids[j2], hids[j1][l] ) ) 
                              better1 += hids[j1][l].second;    }
                    for ( int l = 0; l < hids[j2].isize( ); l++ )
                    {    if ( !BinMember( hids[j1], hids[j2][l] ) ) 
                              better2 += hids[j2][l].second;    }
                    const int min_mult = 4;
                    const int min_count = 2;
                    if ( better1 >= min_count && better1 >= min_mult * better2 )
                         to_delete[j2] = True;    }    }
          EraseIf( paths, to_delete );
          EraseIf( hids, to_delete );
     
          // Report coverage along each path.
	  
	  if ( VERBOSITY >= 1 ){ 
	    for ( int i = 0; i < paths.isize( ); i++ )
	      {    vec<int> p = paths[i];
		vec<int64_t> h;
		int pos = -S.Unibase( p[0] ).isize( );
		for ( int k = 0; k < p.isize( ); k++ )
		  {    int v = p[k];
                    for ( int j = 0; j < hits.isize( ); j++ )
		      {    if ( hits[j].u2 != v ) continue;
			h.push_back( pos + hits[j].pos2 );    }
                    pos += unibases2[v].isize( ) - K2 + 1;    }
		UniqueSort(h);
		cout << "coverage of path [" << i << "] =";
		for ( int j = 0; j < p.isize( ); j++ )
		  {    cout << " " << p[j];
                    if ( big2[ p[j] ] ) cout << "*";    }
		Sort(p);
		Bool doubled = False;
		for ( int j1 = 0; j1 < p.isize( ); j1++ )
		  {    int j2 = p.NextDiff(j1);
                    if ( j2 - j1 == 2 ) doubled = True;
                    j1 = j2 - 1;    }
		if (doubled) cout << " (at repeat)";
		cout << "; len = " << pos << "\n";
		const int dist = 20000;
		const int window = 1000;
		vec<int> cov( dist / window, 0.0 );
		for ( int j = 0; j < h.isize( ); j++ )
		  {    int a = h[j];
                    if ( a < 0 || a > dist ) continue;
                    cov[a/window]++;    }
		cout << "coverage: " << hids[i].size( ) << " =";
		for ( int j = 0; j < cov.isize( ); j++ )
		  cout << " " << cov[j];
		cout << "\n";    }
	  }

          // Compare paths.

          if ( VERBOSITY >= 1 ) cout << "\ncomparison of paths:\n";
          for ( int j1 = 0; j1 < paths.isize( ); j1++ )
          for ( int j2 = j1 + 1; j2 < paths.isize( ); j2++ )
          {    int better1 = 0, better2 = 0;
               for ( int l = 0; l < hids[j1].isize( ); l++ )
               {    if ( !BinMember( hids[j2], hids[j1][l] ) ) 
                         better1 += hids[j1][l].second;    }
               for ( int l = 0; l < hids[j2].isize( ); l++ )
               {    if ( !BinMember( hids[j1], hids[j2][l] ) ) 
		          better2 += hids[j2][l].second;    }
	       if ( VERBOSITY >= 1 ) PRINT4( j1, j2, better1, better2 );    }    }

     // Choose alternatives based on presence in ucores.  First we have to translate
     // to unibases2.

     if ( IsRegularFile( run_dir + "/extended.long.ucores" ) )
     {    vec< vec<int> > ucores, ucores2;
          BinaryReader::readFile( run_dir + "/extended.long.ucores", &ucores );
          for ( int id = 0; id < ucores.isize( ); id++ )
          {    const vec<int>& x = ucores[id];
               if ( x.empty( ) ) continue;
               if ( !to_big[ x[0] ].solo( ) ) continue;
               basevector b1 = unibases1[ x[0] ];
               for ( int j = 1; j < x.isize( ); j++ )
               {    b1.resize( b1.isize( ) - (K1-1) );
                    b1 = Cat( b1, unibases1[ x[j] ] );    }
               int u2 = to_big[ x[0] ][0].first, p2 = to_big[ x[0] ][0].second;
               vec<int> y;
               y.push_back(u2);
               Bool fail = False;
               if ( p2 < 0 ) continue; // This really shouldn't happen!
               for ( int p1 = 0; p1 < b1.isize( ); p1++ )
               {    if ( b1[p1] != unibases2[u2][p2] )
                    {    fail = True;
                         break;    }
                    if ( p1 == b1.isize( ) - 1 ) break;
                    p2++;
                    if ( p2 < unibases2[u2].isize( ) ) continue;
                    p2 = 0;
                    Bool found_next = False;
                    for ( int j = 0; j < nexts2[u2].isize( ); j++ )
                    {    int v2 = nexts2[u2][j];
                         if ( unibases2[v2][K2-1] == b1[p1+1] )
                         {    found_next = True;
                              u2 = v2;
                              p2 = K2-1;
                              y.push_back(u2);
                              break;    }    }
                    if ( !found_next )
                    {    fail = True;
                         break;    }    }
               if ( !fail ) ucores2.push_back(y);    }
          if( VERBOSITY >= 1 ) PRINT2( ucores.size( ), ucores2.size( ) );
          for ( int e = 0; e < S.EdgeN( ); e++ )
          {    gapster& g = S.Gmutable( ).EdgeObjectMutable(e);
               if ( g.ClosureCount( ) <= 1 ) continue;
               vec<int> count( g.ClosureCount( ), 0 );
               for ( int i = 0; i < g.ClosureCount( ); i++ )
               {    for ( int j = 0; j < ucores2.isize( ); j++ )
                    {    if ( ucores2[j].Contains( g.Closure(i).U( ) ) ) 
                              count[i]++;    }    }
               vec<Bool> to_delete( g.ClosureCount( ), False );
               vec<int> ids( g.ClosureCount( ), vec<int>::IDENTITY );
               ReverseSortSync( count, ids );
               int min_count = 5;
               if ( VERBOSITY >= 1 )
               {    cout << "\nlong read votes for edge " << e 
                         << " (" << g.ClosureCount( ) << " closures)" << endl;
                    for ( int j = 0; j < count.isize( ); j++ )
                    {    if ( count[j] == 0 ) break;
                         cout << count[j] << " votes for closure "
                              << ids[j]+1 << endl;    }    }
               if ( count[0] >= min_count )
               {    for ( int j = 1; j < count.isize( ); j++ )
                    {    if ( count[j] > 0 ) continue;
                         to_delete[ ids[j] ] = True;    }    }
               g.RemoveSomeClosures(to_delete);     }
          S.RemoveSubsumedStuff( );
          S.SwallowSimpleGaps( );
          S.BringOutTheDead( );    }

     // Test for singleton plasmids.

     for ( int i = 0; i < S.VertN( ); i++ )
     {    if ( S.From(i).nonempty( ) || S.To(i).nonempty( ) ) continue;
          const uniseq& x = S.Vert(i);
          int u1 = x.U( ).front( ), u2 = x.U( ).back( );
          if ( !nexts2[u2].solo( ) ) continue;
          if ( nexts2[u2][0] != u1 ) continue;
          vec<int> u, overlap;
          u.push_back( u2, u1 );
          overlap.push_back( K2-1 );
          vec<uniseq> c;
          c.push( u, overlap );
          S.Gmutable( ).AddEdge( i, i, gapster(c) );    }

     /*

     // Attempt to pull apart repeats.

     cout << Date( ) << ": pulling apart" << endl;
     for ( int i = 0; i < S.VertN( ); i++ )
     {    int n = S.From(i).size( );
          if ( n <= 1 || n != S.To(i).isize( ) ) continue;
          Bool fail = False;
          for ( int j = 0; j < n; j++ )
          {    if ( S.G( ).EdgeObjectByIndexFrom( i, j ).ClosureCount( ) == 0  )
                    fail = True;
               if ( S.G( ).EdgeObjectByIndexTo( i, j ).ClosureCount( ) == 0  )
                    fail = True;    }
          if (fail) continue;
          if ( VERBOSITY >= 1 )
          {    cout << "\n";
               PRINT2( i, n );
               cout << "x =";
               for ( int j = 0; j < S.Vert(i).N( ); j++ )
                    cout << " " << S.Vert(i).U(j);
               cout << "\n";    }
          // vec<int> vall = S.From(i);
          // vall.append( S.To(i) );
          // vall.push_back(i);
          // UniqueSort(vall);
          // if ( vall.isize( ) != 2*n+1 ) continue;
          vec< triple<int64_t,int,int> > links;
          for ( int j1 = 0; j1 < S.To(i).isize( ); j1++ )
          {    int v1 = S.To(i)[j1];
               for ( int k1 = 0; k1 < S.Vert(v1).N( ); k1++ )
               {    int u1 = S.Vert(v1).U(k1);
                    for ( int r1 = 0; r1 < placements_by_unipath_t[u1].isize( ); 
                         r1++ )
                    {    int64_t id1 = placements_by_unipath_t[u1][r1].first;
                         int pos1 = placements_by_unipath_t[u1][r1].second;
                         if ( !placements_by_unipath_t[u1][r1].third ) continue;
                         int dist1 = unibases2[u1].isize( ) 
                              - pos1 - jbases_t[id1].isize( );
                         for ( int l = k1+1; l < S.Vert(v1).N( ); l++ )
                         {    dist1 += unibases2[ S.Vert(v1).U(l) ].isize( ) 
                                   - K2 + 1;    }
                         int m = 1000000000;
                         const gapster& g1 = S.G( ).EdgeObjectByIndexTo( i, j1 );
                         for ( int m1 = 0; m1 < g1.ClosureCount( ); m1++ )
                         {    const uniseq& x2 = g1.Closure(m1);
                              int d = 0;
                              for ( int l = 1; l < x2.N( ); l++ )
                              {    int w2 = x2.U(l);
                                   d += unibases2[w2].isize( ) - K2 + 1;    }
                              m = Min( m, d );    }
                         dist1 += m;
                         for ( int l = 1; l < S.Vert(i).N( ); l++ )
                              dist1 += unibases2[ S.Vert(i).U(l) ].isize( ) - K2 + 1;
                         const int max_dist = 15000; // ****************************
                         if ( dist1 > max_dist ) continue;
                         int64_t id2 = jpairs_t.getPartnerID(id1);
                         int64_t pid = jpairs_t.getPairID(id1);
                         for ( int s2 = 0; s2 < placements_by_read_t[id2].isize( );
                              s2++ )
                         {    int u2 = placements_by_read_t[id2][s2].first;
                              int pos2 = placements_by_read_t[id2][s2].second;
                              if ( placements_by_read_t[id2][s2].third ) continue;
                              for ( int j2 = 0; j2 < S.From(i).isize( ); j2++ )
                              {    int v2 = S.From(i)[j2];
                                   const gapster& g2
                                        = S.G( ).EdgeObjectByIndexFrom( i, j2 );
                                   for ( int k2 = 0; k2 < S.Vert(v2).N( ); k2++ )
                                   {    if ( S.Vert(v2).U(k2) != u2 ) continue;
                                        int dist = dist1 + pos2;
                                        for ( int l = 1; l < S.Vert(i).N( ); l++ )
                                        {    int v2 = S.Vert(i).U(l);
                                             dist += unibases2[v2].isize( )
                                                  - K2 + 1;    }
                                        for ( int l = 1; l < k2; l++ )
                                        {    int w2 = S.Vert(v2).U(l);
                                             dist += unibases2[w2].isize( )
                                                  - K2 + 1;    }
                                        int m = 1000000000;
                                        for ( int m2 = 0; m2 < g2.ClosureCount( ); 
                                             m2++ )
                                        {    const uniseq& x2 = g2.Closure(m2);
                                             int d = 0;
                                             for ( int l = 1; l < x2.N( ); l++ )
                                             {    int w2 = x2.U(l);
                                                  d += unibases2[w2].isize( )
                                                       - K2 + 1;    }
                                             m = Min( m, d );    }
                                        dist += m;
                                        if ( dist > max_dist ) continue;
                                        int m2 = -1;
                                        if ( VERBOSITY >= 1 ) 
                                        {    PRINT9( j1, j2, m2, id1, id2, pid,
                                                  u1, u2, dist );    }    
                                        links.push( pid, j1, j2 );    }
                                   for ( int m2 = 0; m2 < g2.ClosureCount( ); m2++ )
                                   {    const uniseq& x2 = g2.Closure(m2);
                                        for ( int k2 = 0; k2 < x2.N( ); k2++ )
                                        {    if ( x2.U(k2) != u2 ) continue;
                                             int dist = dist1 + pos2;
                                             for ( int l = 1; l < S.Vert(i).N( ); 
                                                  l++ )
                                             {    int v2 = S.Vert(i).U(l);
                                                  dist += unibases2[v2].isize( )
                                                       - K2 + 1;    }
                                             for ( int l = 1; l < k2; l++ )
                                             {    int w2 = x2.U(l);
                                                  dist += unibases2[w2].isize( )
                                                       - K2 + 1;    }
                                             if ( dist > max_dist ) continue;
                                             if ( VERBOSITY >= 1 ) 
                                             {    PRINT9( j1, j2, m2, id1, id2, pid,
                                                       u1, u2, dist );    }    
                                             links.push( pid, j1, j2 );    
                                                  }    }
                                                       }    }    }    }    }
          UniqueSort(links);
          cout << "\nsummary:\n";
          vec< pair<int,int> > totals;
          for ( int j = 0; j < links.isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < links.isize( ); k++ )
                    if ( links[k].first != links[j].first ) break;
               if ( k - j == 1 )
               {    int64_t pid = links[j].first;
                    int j1 = links[j].second, j2 = links[j].third;
                    PRINT3( pid, j1, j2 );
                    totals.push( j1, j2 );    }
               j = k - 1;    }    
          Sort(totals);
          cout << "\ntotals:\n";
          for ( int j = 0; j < totals.isize( ); j++ )
          {    int k = totals.NextDiff(j);
               int n = k - j;
               int j1 = totals[j].first, j2 = totals[j].second;
               PRINT3( n, j1, j2 );
               j = k - 1;    }    }
     if ( VERBOSITY >= 1 ) cout << "\n";
     */

     // Use GetHomes to eliminate some ambiguities.

     if ( VERBOSITY >= 1 ) cout << "\n";
     cout << Date( ) << ": calling GetHomes" << endl;
     GetHomes( run_dir, K2, unibases2, S, ( VERBOSITY > 0 ) );
     cout << Date() << ": GetHomes finished" << endl;

     // Remove funny links.  We look for cases where there should be a lot more
     // read pairs than there actually are.  I don't think that the math is quite
     // right.

     S.SetToRc(to_rc2); // should not be necessary, working around unknown bug
     for ( int v1 = 0; v1 < S.VertN( ); v1++ )
     {    for ( int j = 0; j < S.From(v1).isize( ); j++ )
          {    int v2 = S.From(v1)[j];
               const gapster& g = S.G( ).EdgeObjectByIndexFrom( v1, j );
               Bool missing = False;
               for ( int l1 = 0; l1 < S.Vert(v1).N( ); l1++ )
               for ( int l2 = 0; l2 < S.Vert(v2).N( ); l2++ )
               {    int u1 = S.Vert(v1).U(l1), u2 = S.Vert(v2).U(l2);
                    const int min_kmers = 1000;
                    if ( unibases2[u1].isize( ) - (K2-1) < min_kmers ) continue;
                    if ( unibases2[u2].isize( ) - (K2-1) < min_kmers ) continue;
                    const double dev_mult = 3.0;
                    int dist = g.MaxLen(dev_mult);
                    if ( g.Closed( ) )
                    {    dist -= unibases2[ g.Closure(0).U( ).front( ) ].isize( );
                         dist -= unibases2[ g.Closure(0).U( ).back( ) ].isize( );   }
                    for ( int m = l1+1; m < S.Vert(v1).N( ); m++ )
                         dist += unibases2[ S.Vert(v1).U(m) ].isize( ) - K2 + 1;
                    for ( int m = 0; m < l2; m++ )
                         dist += unibases2[ S.Vert(v2).U(m) ].isize( ) - K2 + 1;
                    if ( dist > 20000 ) continue; // *******************************
                    double predicted = 0.0, actual = 0.0;
                    for ( int l = 0; l < placements_by_unipath[u1].isize( ); l++ )
                    {    if ( !placements_by_unipath[u1][l].third ) continue;
                         int64_t id1 = placements_by_unipath[u1][l].first;
                         int64_t id2 = jpairs.getPartnerID(id1);
                         int pos1 = placements_by_unipath[u1][l].second;
                         if ( pos1 >= unibases2[u1].isize( ) - (K2-1) ) continue;
                         int low = unibases2[u1].isize( ) + (K2-1)
                              - pos1 - jbases[id1].isize( ) + dist;
                         int high = low + unibases2[u2].isize( );
                         int64_t start = LowerBound( dists, low );
                         int64_t stop = UpperBound( dists, high );
                         if ( start == (int64_t) dists.size( ) ) continue;
                         if ( placements_by_read[id2].empty( ) ) continue;
                         predicted += double(stop-start) / double( dists.size( ) );
                         for ( int m = 0; m < placements_by_read[id2].isize( ); 
                              m++ )
                         {    if ( placements_by_read[id2][m].third ) continue;
                              if ( placements_by_read[id2][m].first != u2 ) continue;
                              int pos2 = placements_by_read[id2][m].second;
                              if ( pos2 + jbases[id2].isize( ) <= K2-1 ) continue;
                              int d = dist + unibases2[u1].isize( ) 
                                   - pos1 - jbases[id1].isize( ) + pos2;
                              double add =
                                   ( double( dists.isize( ) ) - LowerBound(dists,d) )
                                   / double( dists.isize( ) - start );
                              // if ( add > 0 ) PRINT4( pos1, pos2, d, add ); // XXX
                              actual += add;    }    }
                    if ( predicted == 0 ) continue;
                    predicted *= raw2[u2] / raw2[u1];
                    const double min_pred = 20.0;
                    const double div = 10.0;
                    int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
                    int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
                    if ( predicted >= min_pred && actual < predicted/div )
                    {    missing = True;
                         if ( VERBOSITY >= 1 )
                         {    cout << "\nOUCH: missing expected links from\n";
                              S.Vert(v1).Print( cout, K2 );
                              cout << "\nto ";
                              S.Vert(v2).Print( cout, K2 );
                              cout << "\n";    }
                              cout << "dist = " << dist << ", u1 = " << u1 << "[l=" 
                                   << nkmers1 << "], u2 = " << u2 << "[l=" << nkmers2
                                   << "], predicted = " << int(round(predicted)) 
                                   << ", actual = " 
                                   << int(round(actual)) << "\n";    }    }    
               if (missing) S.Gmutable( ).DeleteEdgeFrom( v1, j );    }    }
     S.RemoveSubsumedStuff( );
     S.SwallowSimpleGaps( );
     S.BringOutTheDead( );

     // Enlarge vertices, preparatory to removing fw/rc duplication.  Basically 
     // this tries to suck as much as possible out of the edges and put it into the 
     // vertices instead, which has a canonicalizing effect, thereby making the
     // fw/rc duplication removal better able to recognize duplication.

     for ( int v = 0; v < S.VertN( ); v++ )
     {    uniseq& x = S.VertMutable(v);
          if ( S.From(v).solo( ) )
          {    gapster& g = S.Gmutable( ).EdgeObjectByIndexFromMutable( v, 0 );
               if ( g.ClosureCount( ) == 1 )
               {    uniseq& c = g.ClosureMutable(0);
                    int n = c.N( );
                    for ( int k = 0; k < n - 2; k++ )
                    {    x.UMutable( ).push_back( c.U(1) );
                         c.UMutable( ).erase( c.UMutable( ).begin( ) );
                         x.OverMutable( ).push_back( c.Over(0) );
                         c.OverMutable( ).erase( 
                              c.OverMutable( ).begin( ) );    }    }    }
          if ( S.To(v).solo( ) )
          {    gapster& g = S.Gmutable( ).EdgeObjectByIndexToMutable( v, 0 );
               if ( g.ClosureCount( ) == 1 )
               {    uniseq& c = g.ClosureMutable(0);
                    int n = c.N( );
                    for ( int k = 0; k < n - 2; k++ )
                    {    x.UMutable( ).push_front( c.U( c.U( ).isize( ) - 2 ) );
                         c.UMutable( ).pop_back( );
                         x.OverMutable( ).push_front( c.Over( ).back( ) );
                         c.OverMutable( ).pop_back( );    }    }    }    }

     // Again attempt to remove fw/rc duplication.  May not be necessary.

     if (DELETE_COMPLEMENTS) S.DeleteComplements( );

     // Assess links.  Recalibrate statistical model and exclude alternatives that
     // now appear outside its bounds.

     cout << "\nLINK ASSESSMENT:\n";
     vec<double> global_offbys;
     for ( int u1 = 0; u1 < nuni2; u1++ )
     for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
     {    int u2 = G2I.From(u1)[j];
          int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
          int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
          const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
          int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
          vec<int> actual;
          vec<double> offby;
          for ( int v = 0; v < S.VertN( ); v++ )
          {    const vec<int>& u = S.Vert(v).U( );
               for ( int k1 = 0; k1 < u.isize( ); k1++ )
               {    if ( u[k1] != u1 ) continue;
                    for ( int k2 = k1+1; k2 < u.isize( ); k2++ )
                    {    if ( u[k2] != u2 ) continue;
                         int x = -K2+1;
                         for ( int l = k1+1; l < k2; l++ )
                              x += unibases2[ u[l] ].isize( ) - K2 + 1;
                         actual.push_back(x);    
                         offby.push_back( Abs( double(sep-x)/double(dev) ) );   
                              }    }    }
          SortSync( offby, actual );
          if ( offby.nonempty( ) ) global_offbys.push_back( offby[0] );
          if ( VERBOSITY >= 1 )
          {    cout << "\nu1 = " << u1 << "[l=" << nkmers1 << "], u2 = " << u2
                    << "[l=" << nkmers2 << "], ";
               cout << "sep = " << sep << " +/- " << dev << ", nlinks = " << nlinks;
               cout << ", offby =";
               for ( int l = 0; l < offby.isize( ); l++ )
                    cout << " " << offby[l];
               cout << "\n";    }    }
     Sort(global_offbys);
     if ( global_offbys.empty( ) ) cout << "\nno offbys\n";
     else
     {    const double gtail = 0.9;
          int gstart = int( round( double( global_offbys.size( ) ) * gtail ) );
          const int min_off = 3;
          gstart = Min( gstart, global_offbys.isize( ) - min_off );
          gstart = Max( gstart, 0 );
          cout << "\ntotal offbys = " << global_offbys.size( ) << endl;
          cout << "worst offbys = ...";
          for ( int l = gstart; l < global_offbys.isize( ); l++ )
               cout << " " << global_offbys[l];
          cout << "\n";
          double gtail_mult = 1.5;
          double max_dev = gtail_mult * global_offbys[gstart];
          PRINT(max_dev);
          for ( int v1 = 0; v1 < S.VertN( ); v1++ )
          {    for ( int j = 0; j < S.From(v1).isize( ); j++ )
               {    int v2 = S.From(v1)[j];
                    gapster& g = S.Gmutable( ).EdgeObjectByIndexFromMutable( v1, j );
                    if ( g.ClosureCount( ) <= 1 ) continue;
                    vec<Bool> to_remove( g.ClosureCount( ), False );
                    int u1 = g.Closure(0).U( ).front( );
                    int u2 = g.Closure(0).U( ).back( );
                    int sep, dev;
                    Bool found = False;
                    for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
                    {    if ( u2 == G2I.From(u1)[j] )
                         {    const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
                              sep = l.sep, dev = l.dev;
                              found = True;
                              break;    }    }
                    if ( !found ) continue;
                    for ( int j = 0; j < g.ClosureCount( ); j++ )
                    {    int actual = g.Closure(j).ClosureLen( );
                         if ( Abs( double(actual-sep)/double(dev) ) > max_dev )
                              to_remove[j] = True;    }
                    if ( Sum(to_remove) < g.ClosureCount( ) )
                         g.RemoveSomeClosures(to_remove);     }    }
          S.RemoveSubsumedStuff( );
          S.SwallowSimpleGaps( );
          S.BringOutTheDead( );    }

     // Again use GetHomes to eliminate some ambiguities.

     if ( VERBOSITY >= 1 ) cout << "\n";
     cout << Date( ) << ": calling GetHomes" << endl;
     S.SetToRc(to_rc2); // should not be necessary, working around unknown bug
     GetHomes( run_dir, K2, unibases2, S, ( VERBOSITY > 0 ) );

     // Again filter walks using jumps.

     cout << Date( ) << ": filter walks using jumps" << endl;
     FilterByJumps( S, K2, jbases, jpairs, placements_by_read, 
          placements_by_unipath, VERBOSITY, FILTER_VERBOSITY, FILTER_VERBOSE_U1, 
          FILTER_SHOW_BASES_AFTER );
     S.RemoveSubsumedStuff( );
     S.SwallowSimpleGaps( );
     S.BringOutTheDead( );

     // ============================================================================
     // assembly is complete at this point
     // ============================================================================

     // Print assembly.

     if ( VERBOSITY >= 1 )
     {    cout << "\nFINAL ASSEMBLY\n";
          for ( int i = 0; i < S.VertN( ); i++ )
          {    cout << "\nvertex " << i << endl;
               S.Vert(i).Print( cout, K2 );
               cout << "\n";    }
          for ( int i = 0; i < S.EdgeN( ); i++ )
          {    const gapster& g = S.Edge(i);
               if ( g.ClosureCount( ) == 0 ) continue;
               cout << "\nedge from " << g.Closure(0).U( ).front( ) << " to "
                    << g.Closure(0).U( ).back( ) 
                    << endl;
               for ( int j = 0; j < g.ClosureCount( ); j++ )
               {    cout << "[" << j+1 << "] ";
                    g.Closure(j).Print( cout, K2 );
                    cout << " (" << g.Closure(j).ClosureLen( )
                         << " kmers)\n";    }    }    }

     // Analyze edges.

     if (ANALYZE_EDGES)
     {    cout << "\nANALYSIS OF EDGES UP TO SIX-FOLD\n";
          for ( int i = 0; i < S.EdgeN( ); i++ )
          {    const gapster& g= S.Edge(i);
               int n = g.ClosureCount( );
               const int max_n = 6;
               if ( n >= 2 && n <= max_n )
               {    int u1 = g.Closure(0).U( ).front( );
                    int u2 = g.Closure(0).U( ).back( );
                    cout << "\nedge from " << u1 << " to " << u2 << "\n";
                    for ( int j1 = 0; j1 < n; j1++ )
                    for ( int j2 = j1+1; j2 < n; j2++ )
                    {    cout << "\n[" << j1+1 << "] vs [" << j2+1 << "]\n";
                         uniseq s1 = g.Closure(j1), s2 = g.Closure(j2);
                         s1.TrimEnds(1,1), s2.TrimEnds(1,1);
                         basevector U1 = s1.Bases( ), U2 = s2.Bases( );
                         align a;
                         const int bandwidth = 1000;
                         int errors;
                         SmithWatFreeSym( U1, U2, a, True, True, 1, 1 );
                         // SmithWatBandedA( U1, U2, 0, bandwidth, a, errors );
                         PrintVisualAlignment( 
                              True, cout, U1, U2, a );    }    }    }    }

     // Print graph report.

     if ( DOT != "" ) 
     {    if (DOT_SIMPLE)
          {    String dot_file = DOT.Before( ".dot" ) + ".simple.dot";
               Ofstream( dout, dot_file );
               vec< vec<String> > edge_labels( S.VertN( ) );
               for ( int v = 0; v < S.VertN( ); v++ )
               {    for ( int j = 0; j < S.From(v).isize( ); j++ )
                    {    const gapster& g = S.G( ).EdgeObjectByIndexFrom( v, j );
                         edge_labels[v].push_back( 
                              ToString( g.ClosureCount( ) ) );    }    }
               S.G( ).DOT_el( dout, ( CIRCO ? "circo" : "" ), edge_labels );    }
          vec<String> legends_to_show;
          ParseStringSet( DOT_LEGENDS, legends_to_show );
          BigMapDot( DOT, S, K2, VALIDATE, genome, Ulocs2, legends_to_show,
               CIRCO );    }

     // Compute and write scaffolds
     
     if ( WRITE || WRITE_SNARK )
          BinaryWriter::writeFile( run_dir + "/" + SCAFFOLDS_OUT + ".snark", S );
     if ( !WRITE ) 
     {    cout << Date( ) << ": BigMap done!" << endl;
          return 0;    }
     cout << "\n" << Date() << ": computing scaffolds" << endl; 
     VecEFasta econtigs;
     vec<superb> scaffolds;
     digraphE<sepdev> SGorig;
     S.ComputeScaffolds( scaffolds, econtigs, SGorig );
     
     Assembly A( scaffolds, econtigs );
     A.remove_small_contigs( 1, 1 ); // remove empty contigs;
     A.remove_small_scaffolds( 1 ); // remove empty scaffolds;
     if ( REMOVE_DUPLICATES )
       A.dedup(False); // remove reverse complement and duplicate scaffolds
     //A.dedup2(); // remove likely duplicate scaffolds
     A.check_integrity();

     vec<String> scaffMap = A.scaffMap();
     scaffMap.Print(cout); cout << endl;
     vec<int> rverts; // remaining scaffold vertices
     for ( size_t si = 0; si < scaffMap.size(); si++ )
       rverts.push_back( atoi(scaffMap[si].c_str()) );
     rverts.Print(cout); cout << endl;
     PRINT( SGorig.N() );
     digraphE<sepdev> SG( SGorig.Subgraph(rverts) );
     PRINT( SG.N() );

     cout << "preparing to print" << endl;
     vec<String> v_labels(SG.N()), colors(SG.N());
     for ( int vi = 0; vi < SG.N(); vi++ ){
       int si = atoi( scaffMap[vi].c_str() );
       PRINT(si);
       v_labels[vi] = "sid=" + ToString(si) + "\\nLen=" + ToString( scaffolds.at(si).FullLength() );
     }
     vec< vec<String> > legends, e_labels(SG.N());
     for ( int v = 0; v < SG.N(); v++ ){
       for ( int ii = 0; ii < SG.FromEdgeObj(v).isize(); ii++ ){
	 int ei = SG.FromEdgeObj(v)[ii];
	 sepdev sp = SG.Edges()[ei];
	 e_labels[v].push_back( ToString( sp.Sep() ) + "+/-" + ToString( sp.Dev() ) );
       }
     }
     cout << "printing" << endl;
     String SG_file = sub_dir + "/" + SCAFFOLDS_OUT + ".vertex.dot"; 
     ofstream sgout( SG_file.c_str() );
     SG.DOT_vl( sgout, v_labels, "", legends, colors, e_labels );
     
     cout << "\n" << Date() << ": writing scaffolds" << endl;
     String scaffolds_out_head = sub_dir + "/" + SCAFFOLDS_OUT;
     PRINT( scaffolds_out_head );
     A.WriteAll( scaffolds_out_head ); 

     cout << "\n" << Date( ) << ": time used = " << TimeSince(clock) << endl;
     cout << Date() << ": BigMap done!" << endl;    }
