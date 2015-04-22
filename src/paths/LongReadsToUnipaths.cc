///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LongReadsToUnipaths.  Convert long reads into sequences of unipaths.

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
     CommandArgument_Bool_OrDefault(USE_SHORTEST, True);
     CommandArgument_Bool_OrDefault(CLUSTER_ALIGNS_NEW_CLEAN, False);
     CommandArgument_Bool_OrDefault(BEST_ONLY, True);
     CommandArgument_Bool_OrDefault(NEW_FILTER, True);
     CommandArgument_Bool_OrDefault(SCREEN_NEXTS, True);

     // Diagnostic arguments.

     CommandArgument_Bool_OrDefault_Doc(QLT1, False,
          "align some stuff; forces single thread on main loop");
     CommandArgument_Bool_OrDefault_Doc(DOT1, False, "generate dot file");
     CommandArgument_Bool_OrDefault_Doc(DUMP_LOCAL, False, "print local alignments");
     CommandArgument_Bool_OrDefault_Doc(VERBOSE1, False, "turns on phase 1 logging");
     CommandArgument_Bool_OrDefault_Doc(VALIDATE1, False, "validate some stuff");
     CommandArgument_Bool_OrDefault_Doc(SKIP_SILENT, True, 
          "don't report results for reads that are not corrected");
     CommandArgument_String_OrDefault_Doc(READS_TO_TRACE, "",
          "if specified, trace these reads");
     CommandArgument_String_OrDefault_Doc(READS_TO_USE, "",
          "if specified, use just these reads");
     CommandArgument_Bool_OrDefault(ABBREVIATE_ALIGNMENTS, True);
     CommandArgument_Bool_OrDefault(PRINT_RAW_ALIGNS, False);
     CommandArgument_Bool_OrDefault(PATCHES, False);
     CommandArgument_Int_OrDefault(MIN_PATCH1, 1000000000);
     CommandArgument_Int_OrDefault(MIN_PATCH2, 1000000000);
     CommandArgument_Bool_OrDefault(STANDARD_ALIGNS, True);
     CommandArgument_Bool_OrDefault(FILTER0, True);
     CommandArgument_Bool_OrDefault(FILTER, False);
     CommandArgument_Bool_OrDefault(FILTER2, True);
     CommandArgument_Bool_OrDefault(CLEAN_GRAPH, True);
     CommandArgument_Bool_OrDefault(CYCLIC_BAIL, False);
     CommandArgument_Bool_OrDefault(KILL_INFERIOR_CLUSTERS_NEW, True);
     CommandArgument_Int_OrDefault(MIN_SPREAD1, 0);
     CommandArgument_Int_OrDefault(MAX_GRAPH_SIZE, 0);
     CommandArgument_Bool_OrDefault(SKIP_HIGHLY_COMPLEX_GRAPHS, False );

     EndCommandArguments;

     // Thread control.

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads(NUM_THREADS);

     // Define directories, etc.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String outhead = run_dir + "/" + OUT_HEAD;

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
     if ( VALIDATE1 ) genome.ReadAll( data_dir + "/genome.fastb" );
     const int LG = 12;
     vec< vec< pair<int,int> > > Glocs;
     if ( VALIDATE1 )
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

     // Phase 1.  First go through the reads, running Phase1.  

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
	  SCREEN_NEXTS, MIN_SPREAD1, MAX_GRAPH_SIZE, SKIP_HIGHLY_COMPLEX_GRAPHS, 
	  reads_to_trace, VERBOSE1, SKIP_SILENT, DOT1, 
          QLT1, data_dir, run_dir, ABBREVIATE_ALIGNMENTS, PRINT_RAW_ALIGNS, UNISEQ, 
	  UNISEQ_ID, BESTS, COMPUTED, all, patchers, Hall, ALIGNS_ALL );
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
     cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;    }
