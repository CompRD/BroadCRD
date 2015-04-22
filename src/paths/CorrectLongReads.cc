///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CorrectLongReads.  Correct errors in long reads using unibases derived from short
// reads.  In the first phase, we present an error correction of a given long read
// as a sequence of unipaths, with defined overlaps between successive unipaths in
// the sequence.
//
// UNDER CONSTRUCTION.
//
// KEY THINGS TO DO:
// 1. Clean up, modularize.  Needed to proceed.
// 2. Change phase 2 so it can recognize errors in the FIRST unipath.
// 3. Handle cases where there is a true gap in the unipaths.

// OTHER TO DO:
// - Try forming kmer graph and trimming hanging ends.
// - Look at what happens for small number of reads.
// - Proof of principle by simply replacing reads by corrected ones.
//
// Note that there is a hardcoded special cheat for ecoli.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/CorrectLongReads1.h"
#include "paths/CorrectLongReadsTools.h"
#include "paths/GetNexts.h"
#include "paths/KmerPath.h"
#include "paths/LongReadPatchOptimizer.h"
#include "paths/LongReadTools.h"
#include "paths/PdfEntry.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/UnipathScaffold.h"
#include "paths/Uniseq.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;

     // Core arguments.

     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
          "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(IN_HEAD, "all_reads");
     CommandArgument_String_OrDefault(OUT_HEAD, IN_HEAD + ".long");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Int_OrDefault(KOUT, 640);
     CommandArgument_String_OrDefault(PACBIO_RUNS, "");

     // Heuristics.

     CommandArgument_Int_OrDefault(MIN_SAFE_CN1, 150);
     CommandArgument_Double_OrDefault(MIN_VOTES_TO_PROTECT, 20);
     CommandArgument_Bool_OrDefault(CORRECT_PATCHES, False);
     CommandArgument_Bool_OrDefault(CORRECT_PATCHES_NEW, False);
     CommandArgument_Bool_OrDefault(USE_SHORTEST, False);
     CommandArgument_Bool_OrDefault(CLUSTER_ALIGNS_NEW_CLEAN, False);
     CommandArgument_Bool_OrDefault(BEST_ONLY, False);
     CommandArgument_Bool_OrDefault(NEW_FILTER, False);
     CommandArgument_Bool_OrDefault(SCREEN_NEXTS, False);
     CommandArgument_Int_OrDefault(MIN_TO_PATCH, 5);
     CommandArgument_Bool_OrDefault(NEW_PLUS_RULE, False);

     // Diagnostic arguments.

     CommandArgument_Bool_OrDefault_Doc(QLT1, False,
          "align some stuff; forces single thread on main loop");
     CommandArgument_Bool_OrDefault_Doc(DOT1, False, "generate dot file");
     CommandArgument_Bool_OrDefault_Doc(DOT2, False, "generate dot file");
     CommandArgument_Bool_OrDefault_Doc(DUMP_LOCAL, False, "print local alignments");
     CommandArgument_Bool_OrDefault_Doc(VERBOSE1, False, "turns on phase 1 logging");
     CommandArgument_Bool_OrDefault_Doc(VERBOSE2, False, "turns on phase 2 logging");
     CommandArgument_Bool_OrDefault_Doc(VALIDATE1, False, "validate some stuff");
     CommandArgument_Bool_OrDefault_Doc(VALIDATE2A, False, "validate results");
     CommandArgument_Bool_OrDefault_Doc(VALIDATE2B, False, "validate results");
     CommandArgument_Bool_OrDefault_Doc(VALIDATE_PATCHES, False, "validate patches");
     CommandArgument_Bool_OrDefault_Doc(SKIP_SILENT, True, 
          "don't report results for reads that are not corrected");
     CommandArgument_Bool_OrDefault(PRINT_DISCARDS, False);
     CommandArgument_Bool_OrDefault(PRINT_LM1, False);
     CommandArgument_Bool_OrDefault_Doc(PRINT_MISSING_KMERS, False,
          "print a sample of the missing kmers");
     CommandArgument_Bool_OrDefault(PRINT_MATCHES, False);
     CommandArgument_Int_OrDefault(PATCH_VERBOSITY, 0);
     CommandArgument_Int_OrDefault(PATCH_CORRECT_VERBOSITY, 0);
     CommandArgument_Bool_OrDefault(WRITE_MODIFIED_UNIBASES, False);
     CommandArgument_String_OrDefault_Doc(READS_TO_TRACE, "",
          "if specified, trace these reads");
     CommandArgument_String_OrDefault_Doc(U2, "",
          "if specified, use only these unipaths in phase 2");
     CommandArgument_String_OrDefault_Doc(READS_TO_USE, "",
          "if specified, use just these reads");
     CommandArgument_Int_OrDefault(LOW2, -1);
     CommandArgument_Int_OrDefault(HIGH2, -1);
     CommandArgument_Int_OrDefault(EARLY_EXIT, 0);
     CommandArgument_Bool_OrDefault(ABBREVIATE_ALIGNMENTS, True);
     CommandArgument_Bool_OrDefault(PRINT_RAW_ALIGNS, False);
     CommandArgument_Bool_OrDefault(ANNOUNCE, False);
     CommandArgument_Bool_OrDefault(PATCHES, False);
     CommandArgument_Int_OrDefault(MIN_PATCH1, 1000000000);
     CommandArgument_Int_OrDefault(MIN_PATCH2, 1000000000);
     CommandArgument_Bool_OrDefault(STANDARD_ALIGNS, False);
     CommandArgument_Bool_OrDefault(FILTER0, False);
     CommandArgument_Bool_OrDefault(FILTER, True);
     CommandArgument_Bool_OrDefault(FILTER2, False);
     CommandArgument_Bool_OrDefault(CLEAN_GRAPH, False);
     CommandArgument_Bool_OrDefault(CORRECT_PATCHES_VERBOSE, False);
     CommandArgument_Bool_OrDefault(ANNOUNCE_PATCH_CORRECTION, False);
     CommandArgument_Int_OrDefault(PATCH_U1, -1);
     CommandArgument_Int_OrDefault(PATCH_U2, -1);
     CommandArgument_Bool_OrDefault(DIRECT2, False);
     CommandArgument_Int_OrDefault(SUPER_VERBOSITY2, 0);
     CommandArgument_Bool_OrDefault(PRINT_BASES2, False);
     CommandArgument_Bool_OrDefault(CYCLIC_BAIL, True);
     CommandArgument_Bool_OrDefault(KILL_INFERIOR_CLUSTERS_NEW, False);
     CommandArgument_Int_OrDefault(MAX_GRAPH_SIZE, 0);
     CommandArgument_Bool_OrDefault(SKIP_HIGHLY_COMPLEX_GRAPHS, False );
     CommandArgument_Int_OrDefault(MIN_SPREAD1, 10);
     CommandArgument_Int_OrDefault(MIN_SPREAD2, 10);

     EndCommandArguments;

     // Thread control.

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads(NUM_THREADS);

     // Define directories, etc.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String outhead = run_dir + "/" + OUT_HEAD;
     vec<int> U2list;
     if ( U2 != "" ) ParseIntSet( U2, U2list );

     // Load data.

     cout << Date( ) << ": loading data" << endl;
     vecbasevector unibases( run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K) );
     // vecbasevector longreads( data_dir + "/long_reads_orig.fastb" );
     int nuni = unibases.size( ); 

     // Load fragment reads and create searchable index.

     cout << Date( ) << ": loading fragment reads" << endl;
     vecbasevector fbases( run_dir + "/frag_reads_edit.fastb" );
     vecqualvector fquals( run_dir + "/frag_reads_edit.qualb" );
     PairsManager fpairs( run_dir + "/frag_reads_edit.pairs" );
     fpairs.makeCache( );
     const int F = 20;
     cout << Date( ) << ": building index" << endl;
     vec< kmer<F> > fheads( fbases.size( ) );
     for ( size_t id = 0; id < fbases.size( ); id++ )
          if ( fbases[id].isize( ) >= F ) fheads[id].SetToSubOf( fbases[id], 0 );
     vec<int64_t> fids( fbases.size( ), vec<int64_t>::IDENTITY );
     ParallelSortSync( fheads, fids );

     // Load long reads.

     cout << Date( ) << ": loading long reads" << endl;
     vecbasevector longreads; 
     if ( PACBIO_RUNS != "" ) 
     {    vec<int> runs;
          ParseIntSet( PACBIO_RUNS, runs );
          for ( int i = 0; i < runs.isize( ); i++ )
          {    vecbasevector x;
               String pb_pre = ToString( runs[i] );
               pb_pre.resize(2);
               pb_pre = "0" + pb_pre;
               String dir = "/seq/pacbio_results/userdata/jobs/" + pb_pre + "/0"
                   + ToString( runs[i] ) + "/data";
               String fn;
               if ( IsRegularFile( dir + "/filtered_subreads.fa" ) )
                    fn = dir + "/filtered_subreads.fa";
               else fn = dir + "/filtered_subreads.fasta";
	       FetchReads( x, 0, fn );
	       longreads.Append(x);    }    }
     else longreads.ReadAll( run_dir + "/long_reads_orig.fastb" );
     int nreads = longreads.size( );
     cout << Date( ) << ": have " << ToStringAddCommas(nreads) << " reads" << endl;

     // Create ancillary data for unibases.

     cout << Date( ) << ": making nexts" << endl;
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc );

     // Identify unipaths that may have copy number greater than one.

     vec<Bool> cn_plus;
     if (PATCHES)
     {    cn_plus.resize( nuni, False );
          String KS = ToString(K);
          vec<double> CN_raw;
          digraphE<linklet> G;
          BinaryReader::readFile( run_dir + "/" + IN_HEAD
               + ".unipaths.cn_raw.k" + KS, &CN_raw );
          vec< vec<int> > from(nuni), to(nuni);
          vec< vec<int> > from_edge_obj(nuni), to_edge_obj(nuni);
          vec<linklet> edges;
          String line;
          fast_ifstream in( run_dir + "/" + IN_HEAD + ".unibases.k" 
               + KS + ".predicted_gaps.txt" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               int u1, u2, sep, dev, nlinks;
               istringstream iline( line.c_str( ) );
               iline >> u1 >> u2 >> sep >> dev >> nlinks;
               if ( u2 == to_rc[u1] ) continue; 
               from[u1].push_back(u2), to[u2].push_back(u1);
               from_edge_obj[u1].push_back( edges.size( ) );
               to_edge_obj[u2].push_back( edges.size( ) );
               edges.push( sep, dev, nlinks, 0 );    }
          for ( int u = 0; u < nuni; u++ )
          {    SortSync( from[u], from_edge_obj[u] );
               SortSync( to[u], to_edge_obj[u] );    }
          G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );
          for ( int u = 0; u < nuni; u++ )
          {    Bool first = True;
               const int min_links = 3;
               const double min_mult = 1.8;
               for ( int j = 0; j < G.From(u).isize( ); j++ )
               {    int v = G.From(u)[j];
                    const linklet& l = G.EdgeObjectByIndexFrom( u, j );
                    if ( l.nlinks < min_links ) continue;
                    if ( unibases[v].size( ) < unibases[u].size( ) ) continue;
                    if ( CN_raw[u] < min_mult * CN_raw[v] ) continue;
                    cn_plus[u] = True;
                    Bool cn_plus_verbose = False;
                    if (cn_plus_verbose)
                    {    if (first)
                         {    cout << "\nu = " << u << "[kmers=" 
                                   << unibases[u].isize( ) - (K-1) 
                                   << ",cn=" << CN_raw[u] << "]" << endl;    }
                         first = False;
                         cout << "--(" << l.nlinks << ")--> " << v << "[kmers=" 
                              << unibases[v].isize( ) - (K-1) << ",cn=" 
                              << CN_raw[v] << "]" << endl;    }    }    }    }

     // Heuristics.

     heuristics heur;
     heur.L = 11;
     heur.min_overlap_to_see_other = 250;
     heur.max_offset_diff_for_other = 25;
     heur.min_rdist_for_other = 10;
     heur.delta_r_max = 350;
     heur.max_delta_ratio = 1.5;
     heur.delta_sub = 10;
     heur.min_spread = 10;
     heur.min_ratio_to_kill = 5.0;
     heur.max_paths = 10000;
     heur.max_iterations = 10000;
     heur.min_safe_cn1 = MIN_SAFE_CN1;
     heur.min_win_ratio = 2.0;
     heur.min_votes = 4.0;
     heur.min_votes_to_protect = MIN_VOTES_TO_PROTECT;

     // For evaluation purposes, circularize the genome and hash it.

     vecbasevector genome, genome2;
     if ( VALIDATE1 || VALIDATE2A || VALIDATE2B || VALIDATE_PATCHES ) 
          genome.ReadAll( data_dir + "/genome.fastb" );
     const int LG = 12;
     vec< vec< pair<int,int> > > Glocs;
     if ( VALIDATE1 || VALIDATE_PATCHES )
     {    genome2.resize( genome.size( ) );
          for ( size_t j = 0; j < genome.size( ); j++ )
               genome2[j] = Cat( genome[j], genome[j] );
          Glocs.resize( IPow( 4, LG ) );
          for ( size_t i = 0; i < genome2.size( ); i++ )
          {    for ( int j = 0; j <= genome2[i].isize( ) - LG; j++ )
               {    int n = KmerId( genome2[i], LG, j );
                    Glocs[n].push( i, j );    }    }    }

     // Hash the unibases.

     cout << Date( ) << ": hashing unibases" << endl;
     int L = heur.L;
     vec< vec< pair<int,int> > > Ulocs( IPow( 4, L ) );
     for ( size_t i = 0; i < unibases.size( ); i++ ) 
     {    for ( int j = 0; j <= unibases[i].isize( ) - L; j++ ) 
          {    int n = KmerId( unibases[i], L, j );
               Ulocs[n].push( i, j );    }    } 

     // Setup for phase 1.

     cout << Date( ) << ": traversing the reads" << endl;
     vec<int> reads_to_use, reads_to_trace;
     if ( READS_TO_USE == "" ) 
     {    for ( int j = 0; j < nreads; j++ )
               reads_to_use.push_back(j);    }
     else ParseIntSet( READS_TO_USE, reads_to_use );
     for ( int i = 0; i < reads_to_use.isize( ); i++ )
     {    ForceAssertGe( reads_to_use[i], 0 );
          ForceAssertLt( reads_to_use[i], (int) longreads.size( ) );    }
     ParseIntSet( READS_TO_TRACE, reads_to_trace );

     // Load filled fragments and map them back to the unibases.
     
     vec< vec<int> > fillseqs;
     vec< vec<int> > nexts_count(nuni);
     if (SCREEN_NEXTS)
     {
     const int K = 96;
     vecbasevector fills( run_dir + "/filled_reads.fastb" );
     vec< triple<kmer<K>,int,int> > kmers_plus;
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          kmer<K> x;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j );
               kmers_plus[r].first = x;
               kmers_plus[r].second = i;
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     vec< vec<int> > fillx( fills.size( ) );
     PRINT( fills.size( ) );
     #pragma omp parallel for
     for ( size_t id = 0; id < fills.size( ); id++ )
     {    const basevector& f = fills[id];
          if ( f.isize( ) <= K ) continue;
          kmer<K> x;
          x.SetToSubOf( f, 0 );
          int64_t p = BinPosition( kmers, x );
          if ( p < 0 ) continue;
          int u = kmers_plus[p].second, pos = kmers_plus[p].third;
          if ( pos + f.isize( ) <= unibases[u].isize( ) ) continue;
          vec< pair<int,int> > locs;
          locs.push( u, pos );
          for ( int j = 1; j <= f.isize( ) - K; j++ )
          {    x.SetToSubOf( f, j );
               p = BinPosition( kmers, x );
               if ( p < 0 ) break;
               u = kmers_plus[p].second, pos = kmers_plus[p].third;
               locs.push( u, pos );    }
          if ( locs.isize( ) < f.isize( ) - K + 1 ) continue;
          vec<int> us;
          Bool fail = False;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < locs.isize( ); j++ )
                    if ( locs[j].first != locs[i].first ) break;
               if ( i > 0 && locs[i].second != 0 ) fail = True;
               if ( i > 0 && locs[i-1].second
                    != unibases[ locs[i-1].first ].isize( ) - K )
               {    fail = True;    }
               for ( int k = i+1; k < j; k++ )
                    if ( locs[k].second != locs[k-1].second + 1 ) fail = True;
               if (fail) break;
               us.push_back( locs[i].first );
               i = j - 1;    }
          if (fail) continue;
          #pragma omp critical
          {    fillseqs.push_back(us);    }    }
     int nf = fillseqs.size( );
     for ( int i = 0; i < nf; i++ )
     {    vec<int> x = fillseqs[i];
          x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = to_rc[ x[j] ];
          fillseqs.push_back(x);    }
     PRINT( fillseqs.size( ) );
     Bool fverbose = False;
     if (fverbose)
     {    for ( int j = 0; j < fillseqs.isize( ); j++ )
          {    cout << "[" << j << " fill]";
               const vec<int>& x = fillseqs[j];
               for ( int i = 0; i < x.isize( ); i++ )
                    cout << " " << x[i];
               cout << "\n";    }    }

     // Assess branches using filled fragments.

     cout << Date( ) << ": creating ffpairs" << endl;
     vec< pair<int,int> > ffpairs;
     for ( int i = 0; i < fillseqs.isize( ); i++ )
     {    const vec<int>& x = fillseqs[i];
          for ( int j = 0; j < x.isize( ) - 1; j++ )
               ffpairs.push( x[j], x[j+1] );    }
     cout << Date( ) << ": sorting" << endl;
     ParallelSort(ffpairs);
     cout << Date( ) << ": cataloging" << endl;
     for ( int u = 0; u < nuni; u++ )
          nexts_count[u].resize( nexts[u].size( ), 0 );
     for ( int i = 0; i < ffpairs.isize( ); i++ )
     {    int j = ffpairs.NextDiff(i);
          int u = ffpairs[i].first, v = ffpairs[i].second;
          for ( int l = 0; l < nexts[u].isize( ); l++ )
               if ( nexts[u][l] == v ) nexts_count[u][l] = j - i;
          i = j - 1;    }
     }

     // Phase 1.  First go through the reads, running Phase1.  Then form
     // consensus patches.

     vec<uniseq> UNISEQ(nreads);
     vec< vec<int> > UNISEQ_ID(nreads);
     vec< vec< vec<int> > > BESTS(nreads);
     vec<Bool> COMPUTED( nreads, False );
     vecbasevector all;
     vec<GapPatcher> patchers;
     cout << "\n" << Date( ) << ": START PHASE 1, PASS 1" << endl;
     vec< digraphVE<int,int> > Hall(nreads);
     vec< vec< pair< int, vec< pair<int,int> > > > > ALIGNS_ALL(nreads);
     Phase1( NUM_THREADS, unibases, to_rc, nexts, nexts_count, K, L, Ulocs, 
          longreads, reads_to_use, heur, USE_SHORTEST, CLUSTER_ALIGNS_NEW_CLEAN,
          BEST_ONLY, PATCHES, MIN_PATCH1, MIN_PATCH2, STANDARD_ALIGNS, FILTER0, 
          FILTER, FILTER2, CLEAN_GRAPH, CYCLIC_BAIL, KILL_INFERIOR_CLUSTERS_NEW, 
	  SCREEN_NEXTS, MIN_SPREAD1, MAX_GRAPH_SIZE, SKIP_HIGHLY_COMPLEX_GRAPHS, reads_to_trace, 
          VERBOSE1, SKIP_SILENT, DOT1, QLT1, data_dir, run_dir, 
          ABBREVIATE_ALIGNMENTS, PRINT_RAW_ALIGNS, UNISEQ, UNISEQ_ID, BESTS, 
          COMPUTED, all, patchers, Hall, ALIGNS_ALL );
     if (WRITE)
     {    cout << Date( ) << ": writing files" << endl;
          vec< vec<int> > raw_longs(nreads);
          for ( int id = 0; id < nreads; id++ )
               raw_longs[id] = UNISEQ[id].U( );
          BinaryWriter::writeFile( outhead + ".raw_longs", raw_longs );
          BinaryWriter::writeFile( outhead + ".read_paths", BESTS );
          BinaryWriter::writeFile( outhead + ".read_graph", Hall );
          BinaryWriter::writeFile( outhead + ".read_aligns", ALIGNS_ALL );
          cout << Date( ) << ": done writing files" << endl;    }
     if ( EARLY_EXIT == 1 ) return 0;
     if (PATCHES)
     {     vec<Bool> to_delete( patchers.size( ) );
          for ( int i = 0; i < patchers.isize( ); i++ )
          {    if ( PATCH_U1 >= 0 && patchers[i].t1 != PATCH_U1 ) 
                    to_delete[i] = True;
               if ( PATCH_U2 >= 0 && patchers[i].t2 != PATCH_U2 ) 
                    to_delete[i] = True;   }
          EraseIf( patchers, to_delete );
          vec<basevector> bpatches;
          if ( K == 96 )
	  { cout << Date() << ": Building Patches <96>" << endl;   
	    BuildPatches<96>( unibases, nexts, cn_plus, L, Ulocs, patchers, 
                    CORRECT_PATCHES, CORRECT_PATCHES_NEW,
                    CORRECT_PATCHES_VERBOSE, ANNOUNCE_PATCH_CORRECTION, genome2, 
                    PATCH_VERBOSITY, VALIDATE_PATCHES, PATCH_CORRECT_VERBOSITY, LG, 
                    Glocs, data_dir, run_dir, bpatches, MIN_TO_PATCH,
                    fbases, fquals, fpairs, fheads, fids, NEW_PLUS_RULE );    
	    cout << Date() << ": Patches built <96>" << endl;
	  }
          else if ( K == 640 )
          { cout << Date() << ": Building Patches <640>" << endl;    
	    BuildPatches<640>( unibases, nexts, cn_plus, L, Ulocs, patchers, 
                    CORRECT_PATCHES, CORRECT_PATCHES_NEW,
                    CORRECT_PATCHES_VERBOSE, ANNOUNCE_PATCH_CORRECTION, genome2, 
                    PATCH_VERBOSITY, VALIDATE_PATCHES, PATCH_CORRECT_VERBOSITY, LG, 
                    Glocs, data_dir, run_dir, bpatches, MIN_TO_PATCH,
                    fbases, fquals, fpairs, fheads, fids, NEW_PLUS_RULE );    
	    cout << Date() << ": Patches built <640>" << endl;
	  }
          else ForceAssert( 0 == 1 );
          if ( EARLY_EXIT == 2 ) return 0;
          Sort(bpatches);

          // Merge patches into unibases.

	  cout << Date() << ": Merging patches into unibases" << endl;
          vecbasevector growl(unibases);
          for ( int j = 0; j < bpatches.isize( ); j++ )
               growl.push_back_reserve( bpatches[j] );
          for ( int id1 = 0; id1 < nuni; id1++ ) 
          {    for (int j = 0; j < nexts[id1].isize(); j++) 
               {    int id2 = nexts[id1][j];
                    basevector b = unibases[id1];
                    b.resize( b.size( ) + 1 );
                    b.Set( b.size( ) - 1, unibases[id2][K-1] );
                    growl.push_back_reserve(b);    }    }
          vecKmerPath newpaths, newpathsrc, newunipaths;
          vec<tagged_rpint> newpathsdb, newunipathsdb;
          ReadsToPathsCoreY( growl, K, newpaths, newpathsrc, newpathsdb,
               run_dir + "/CLR", NUM_THREADS );
          Unipath( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb );
          KmerBaseBroker newkbb( K, newpaths, newpathsrc, newpathsdb, growl );
          vecbasevector newunibases;
          for ( size_t i = 0; i < newunipaths.size( ); i++ )
               newunibases.push_back_reserve( newkbb.Seq( newunipaths[i] ) );
          cout << Date( ) << ": old unibases = " << unibases.size( )
               << ", new unibases = " << newunibases.size( ) << endl;
          unibases = newunibases;

          // Analyze patched unibases.

	  cout << Date() << ": Analyzing patched unibases" << endl;
          if ( genome.size( ) > 0 )
	    {    vec< vec<ho_interval> > cov( genome.size( ) );
	      for ( int u = 0; u < (int) unibases.size( ); u++ )
               {    vec<placementx> places 
                         = FindGenomicPlacements( unibases[u], LG, genome2, Glocs );
                    for ( int j = 0; j < places.isize( ); j++ )
                    {    const placementx& p = places[j];
                         int pos = p.pos, Pos = p.pos + unibases[u].isize( ) - (K-1);
                         if ( p.pos >= genome[p.g].isize( ) ) continue;
                         if ( Pos <= genome[p.g].isize( ) )
                              cov[p.g].push( pos, Pos );
                         else
                         {    cov[p.g].push( pos, genome[p.g].size( ) );
                              cov[p.g].push( 
                                   0, Pos - genome[p.g].isize( ));    }    }    }
               cout << "\ngaps in coverage of genome by patched unibases:\n";
               for ( int g = 0; g < (int) genome.size( ); g++ )
               {    vec<ho_interval> un;
                    Uncovered( genome[g].size( ), cov[g], un );
                    for ( int j = 0; j < un.isize( ); j++ )
                         cout << g << "." << un[j] << "\n";    }
               cout << "\n";    }

          // Write patched unibases.

          if (WRITE_MODIFIED_UNIBASES)
          {    unibases.WriteAll( run_dir + "/" + 
                    IN_HEAD + ".CLR_modified.unibases.k" + ToString(K) );    }
          return 0;    
     }

     // Phase 2.  Reanalyze data.  Pile up the corrected reads.

     cout << Date( ) << ": phase 2, reanalyzing" << endl;
     vec< vec<int> > hits( unibases.size( ) );
     for ( int id = 0; id < nreads; id++ )
     {    if ( !COMPUTED[id] ) continue;
          for ( int j = 0; j < UNISEQ[id].N( ); j++ )
               hits[ UNISEQ[id].U(j) ].push_back(id);    }
     for ( size_t i = 0; i < unibases.size( ); i++ )
          UniqueSort( hits[i] );
     uniseq dummy;
     dummy.SetUnibases(unibases);

     // Phase 2a.  Pile up along a given copy number one unipath.

     cout << Date( ) << ": start phase2a" << endl;
     vec<int> uni1;
     int kmers_low = ( LOW2 >= 0 ? LOW2 : heur.min_safe_cn1 - (K-1) );
     int kmers_high = ( HIGH2 >= 0 ? HIGH2 : 1000000000 );
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    int nkmers = unibases[i].isize( ) - (K-1);
          if ( nkmers >= kmers_low && nkmers <= kmers_high )
          {    uni1.push_back(i);    }    }
     int count2 = uni1.size( ), done2_ptr = 0;
     vec<Bool> done2( count2, False );
     vec<String> reports2(count2);
     vecbasevector new_stuff;
     vec< vec<int> > right_exts;
     for ( int ui = 0; ui < count2; ui++ )
     {    int u = uni1[ui];
          if ( U2 != "" && !BinMember( U2list, u ) ) { done2[ui] = True; continue; }
          if ( hits[u].empty( ) && hits[ to_rc[u] ].empty( ) )
          {    done2[ui] = True; continue;    }
          double clock2 = 0.0;
          if (ANNOUNCE)
          {
               #pragma omp critical
               {    clock2 = WallClockTime( );
                    cout << Date( ) << ": begin " << u << " (bases=" 
                         << unibases[u].size( ) << ")" << endl;    }    }
          ostringstream hout;
          Phase2( u, unibases, Ulocs, to_rc, longreads, K, hits, UNISEQ, UNISEQ_ID,
               ( DIRECT2 ? cout : hout ), DOT2, PRINT_MATCHES, PRINT_LM1, 
               DUMP_LOCAL, PRINT_DISCARDS, VALIDATE1, data_dir, genome2, Glocs, LG, 
               heur, new_stuff, right_exts, NEW_FILTER, ANNOUNCE, SUPER_VERBOSITY2,
               PRINT_BASES2, MIN_SPREAD2 );
          if (ANNOUNCE)
          {
               #pragma omp critical
               {    cout << Date( ) << ": end " << u 
                         << ", time used = " << TimeSince(clock2) << endl;    }    }
          reports2[ui] = hout.str( ), done2[ui] = True;
          if ( VERBOSE2 && !DIRECT2 )
          {    
               #pragma omp critical
               {    if ( done2[done2_ptr] )
                    {    while(1)
                         {    cout << reports2[done2_ptr++];
                              if ( done2_ptr == done2.isize( ) || !done2[done2_ptr] )
                                   break;    }    }    }    }    }
     Sort(right_exts);

     sort( new_stuff.begin( ), new_stuff.end( ) );
     if ( VERBOSE2 && !DIRECT2 )
     {    for ( int j = done2_ptr; j < done2.isize( ); j++ )
               cout << reports2[j];
          cout << "\n";    }

     // Validate results.

     if (VALIDATE2A) 
     {    cout << "\nvalidation of new stuff:\n";
          const int min_mult = 1;
          Validate( new_stuff, unibases, genome, min_mult, 
               PRINT_MISSING_KMERS );    }
     if (VALIDATE2B) 
     {    cout << "\nvalidation of all:\n";
          const int min_mult = 2;
          Validate( all, unibases, genome, min_mult, PRINT_MISSING_KMERS );    }

     // Integrate in new unipaths and write results.

     if (WRITE) 
     {    cout << Date( ) << ": writing right_exts" << endl;
          BinaryWriter::writeFile( ( outhead + ".right_exts" ).c_str( ), 
               right_exts );    }

     // Done.

     cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;    }
