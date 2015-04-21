///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/AddLongReads.h"
#include "paths/long/DataSpec.h"
#include "paths/long/ImproveLongHyper.h"
#include "paths/long/KmerCount.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector6.h"
#include "paths/long/fosmid/FosmidVector.h"
#include "random/Bernoulli.h"

void ImproveLongHyper( const String& SAMPLE, const String& X, 
     SupportedHyperBasevector& shb, const String& TMP, const String& READS,
     const long_data_spec& spec, const long_heuristics& heur, unsigned nThreads, 
     const long_logging_control& log_control, const long_logging& logc,
     const int start_step, bool useOldLRPMethod )
{    
     if ( start_step >= 6 && start_step <= 10 )
     {    cout << "Because of the role of shb0, it won't work to start at any\n"
               << "step between 6 and 10.  Sorry.\n";
          cout << "Abort." << endl;
          Scram(1);    }

     shb.TestValid(logc);

     vecbasevector bases;
     vecqualvector quals;
     if ( READS != "" && !READS.Contains( ".fastb", -1 ) )
     {    bases.ReadAll( TMP + "/frag_reads_orig.fastb" );
          quals.ReadAll( TMP + "/frag_reads_orig.qualb" );    }

     // Initial steps.

     const double junk_ratio = 10.0;
     const int max_del = 1000;
     vec< vec<int> > comps;
     if ( start_step <= 2 )
     {    
          // Remove Fosmid vector.

          if ( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && X == ""
               && heur.REMOVE_FOSMID_VECTOR )
          {    RemoveFosmidVector( shb, log_control, logc );    }

          // Implement REQUIRE_EDGE_MATCH.

          if ( heur.REQUIRE_EDGE_MATCH != "" )
          {    cout << Date( ) << ": implementing REQUIRE_EDGE_MATCH" << endl;
               vecbasevector rbases(heur.REQUIRE_EDGE_MATCH);
               vecbasevector edges;
               for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
                    edges.push_back( shb.EdgeObject(e) );
               vec<Bool> marked( edges.size( ), False );
               const int K = 100;
               ForceAssertEq( K, shb.K( ) );
               vecbasevector all(rbases);
               all.Append(edges);
               vec< triple<kmer<K>,int,int> > kmers_plus;
               cout << Date( ) << ": making kmer lookup" << endl;
               MakeKmerLookup2( all, kmers_plus );
               cout << Date( ) << ": marking edges" << endl;
               for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
               {    int64_t j;
                    for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    Bool valid = False;
                    for ( int64_t k = i; k < j; k++ )
                    {    if ( kmers_plus[k].second < (int) rbases.size( ) ) 
                              valid = True;    }
                    if (valid)
                    {    for ( int64_t k = i; k < j; k++ )
                         {    if ( kmers_plus[k].second >= (int) rbases.size( ) ) 
                              {    int id = kmers_plus[k].second 
                                        - (int) rbases.size( );
                                   marked[id] = True;    }    }    }
                    i = j - 1;    }
               PRINT2( shb.EdgeObjectCount( ), Sum(marked) );
               vec<int> dels;
               const int min_kmers_to_del = 100;
               for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               {    if ( !marked[e] && shb.EdgeLengthKmers(e) >= min_kmers_to_del ) 
                         dels.push_back(e);    }
               // cout << "deleting edges " << printSeq(dels) << endl;
               shb.DeleteEdges(dels);
               shb.RemoveUnneededVertices( );
               shb.RemoveDeadEdgeObjects( );
               shb.RemoveEdgelessVertices( );    }

          // Kill homopolymeric edges.

          if (heur.KILL_HOMOPOLYMER_EDGES)
          {    if (logc.STATUS_LOGGING)
                    cout << Date( ) << ": deleting homopolymeric edges" << endl;
               const int min_other = 5;
               const int max_homop_solo = 50;
               vec<int> dels;
               for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               {    vec<int> count(4,0);
                    const basevector& E = shb.EdgeObject(e);
                    for ( int i = 0; i < E.isize( ); i++ )
                         count[ E[i] ]++;
                    if ( Max(count) > E.isize( ) - min_other )
                         dels.push_back(e);
                    else
                    {    int maxh = 0;
                         for ( int i = 0; i < E.isize( ); i++ )
                         {    int j;
                              for ( j = i + 1; j < E.isize( ); j++ )
                                   if ( E[j] != E[i] ) break;
                              maxh = Max( maxh, j - i );    }
                         if ( E.isize( ) == shb.K( ) && maxh > max_homop_solo )
                              dels.push_back(e);    }    }
               shb.DeleteEdges(dels);
               shb.RemoveDeadEdgeObjects( );
               shb.RemoveEdgelessVertices( );
               shb.TestValid(logc);
               if (logc.STATUS_LOGGING)
               {    cout << Date( ) << ": now there are " << shb.EdgeObjectCount( )
                         << " edges" << endl;    }    }

          // Delete reverse complement components.

          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": initially there are " << shb.EdgeObjectCount( )
                    << " edges" << endl;   } 
          if (logc.STATUS_LOGGING)
               cout << Date( ) << ": deleting reverse complement components" << endl;
          shb.DeleteReverseComplementComponents(logc,-1);

          // Remove hanging ends.  

          shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );

          // Remove small components.

          double yclock = WallClockTime( );
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": see " << shb.EdgeObjectCount( ) << " edges, " 
                    << shb.UsedCount( ) << " used, removing small components" 
                    << endl;    }
          shb.RemoveSmallMostlyAcyclicComponents( logc );
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": now there are " << shb.EdgeObjectCount( )
                    << " edges" << endl;    }
          shb.TestValid(logc);
          REPORT_TIME( yclock, "used removing small components" );
     
          // Clean up.
     
          shb.RemoveDeadEdgeObjects( );
          double s1clock = WallClockTime( );
          shb.RemoveEdgelessVertices( );
          shb.RemoveUnneededVertices( );
          REPORT_TIME( s1clock, "used removing small components tail 1" );
          shb.RemoveDeadEdgeObjects( );
          shb.FixWeights(logc);
          shb.TestValid(logc);

          if (heur.DELETE_WEAK)
          {    shb.DeleteWeakEdges(logc);
               shb.RemoveSmallMostlyAcyclicComponents( logc );    }

          shb.DumpFilesStandard( log_control, logc, 2 );    }

     // Pop bubbles, then delete low-coverage edges.

     if ( start_step <= 3 )
     {    
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": see " << shb.EdgeObjectCount( ) 
                    << " edges, " << shb.UsedCount( ) << " used, popping bubbles" 
                    << endl;    }
          const double low_cov = 2.0;
          const double max_pop_del2 = 5.0;
          const double min_cov_ratio = 5.0;
          const int delta_kmers = 2;
          SupportedHyperBasevector::BubbleParams parms(low_cov,max_pop_del2,
                                                          min_cov_ratio,delta_kmers);
          if ( heur.POP_BUBBLES && 
               !( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && X == "" ) )
          {    shb.PopBubbles( parms, nThreads, log_control, logc, heur );    }
          shb.DumpFilesStandard( log_control, logc, 3 );    }
     if ( start_step <= 4 )
     {    shb.DeleteLowCoverage( heur, log_control, logc );
          shb.TestValid(logc);
          if (heur.EXTRA_STEPS) 
          {    shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
               shb.DeleteReverseComplementComponents(logc);
               shb.DeleteLowCoverage( heur, log_control, logc );
               shb.TestValid(logc);
               shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
               shb.DeleteLowCoverage( heur, log_control, logc );
               shb.TestValid(logc);    }
          shb.DumpFilesStandard( log_control, logc, 4 );    }

     // Unwind assembly.

     SupportedHyperBasevector shb0;
     if ( start_step <= 5 )
     {    double uclock = WallClockTime( );
          shb0 = shb;
          REPORT_TIME( uclock, "used in prelude to unwinding" );
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": before unwinding there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;    }
          shb.UnwindAssembly( log_control, logc );
          shb.Bootstrap( shb0, logc );
          shb.TestValid(logc);
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": after unwinding there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;    }
          shb.DumpFilesStandard( log_control, logc, 5 );    }

     // Unwind in reverse.

     if ( start_step <= 6 ) // CAN'T START AT THIS STEP
     {    double b2clock = WallClockTime( );
          shb.Reverse( );
          REPORT_TIME( b2clock, "used after bootstrapping" );
          shb.UnwindAssembly( log_control, logc );
          double rclock = WallClockTime( );
          shb.Reverse( );
          REPORT_TIME( rclock, "unwinding in reverse activity" );
          shb.Bootstrap( shb0, logc );
          shb.TestValid(logc);
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": after unwinding in reverse there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;    }
          shb.DumpFilesStandard( log_control, logc, 6 );    }

     // Pull apart simple branches.

     const double min_weight_split = 10.0;
     if ( start_step <= 7 ) // CAN'T START AT THIS STEP
     {    shb.PullApart( min_weight_split, logc );
          shb.DumpFilesStandard( log_control, logc, 7 );    }

     // Try unwinding again.  Note that this does not do a lot.

     if ( start_step <= 8 ) // CAN'T START AT THIS STEP
     {    if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": before unwinding again there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;    }
          if ( !( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && X == "" ) )
          {    shb.UnwindAssembly( log_control, logc );
               shb.Bootstrap( shb0, logc );    }
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": after unwinding there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;    }
          shb.TestValid(logc);
          shb.DeleteReverseComplementComponents(logc);
          shb.DumpFilesStandard( log_control, logc, 8 );    }

     // Pull apart simple branches.

     if ( start_step <= 9 ) // CAN'T START AT THIS STEP
     {    shb.PullApart( min_weight_split, logc );
          shb.DumpFilesStandard( log_control, logc, 9 );    }

     // Try unwinding again.  Note that this does not do a lot.

     if ( start_step <= 10 ) // CAN'T START AT THIS STEP
     {    
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": before unwinding again there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;    }
          if ( !( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && X == "" ) )
          {    shb.UnwindAssembly( log_control, logc );
               shb.Bootstrap( shb0, logc );    }
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": after unwinding there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;    }
          shb.TestValid(logc);
          shb.DeleteReverseComplementComponents(logc);
          shb.DumpFilesStandard( log_control, logc, 10 );    }

     // Chunk edges.

     if ( heur.CHUNK_EDGES && start_step <= 11 )
     {    shb.ChunkEdges(logc);
          shb.DumpFilesStandard( log_control, logc, 11 );    }

     // Wordify.

     if ( heur.WORDIFY && start_step <= 12 )
     {    if (logc.STATUS_LOGGING) cout << Date( ) << ": wordifying" << endl;
          int64_t kmers1 = 0, kmers2 = 0;
          for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               kmers1 += shb.EdgeLengthKmers(e);
          shb.WordifyAlt2(log_control,logc);
          for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               kmers2 += shb.EdgeLengthKmers(e);
          shb.TestValid(logc);
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": after wordifying there are "
                    << shb.EdgeObjectCount( ) << " edges" << endl;
               cout << Date( ) << ": total kmers changed from " 
                    << ToStringAddCommas(kmers1)
                    << " to " << ToStringAddCommas(kmers2) << endl;    }
          double cxclock = WallClockTime( );
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": deleting reverse complement components" 
                    << endl;    }
          shb.DeleteReverseComplementComponents(logc);
          REPORT_TIME( cxclock, "used deleting rc components" );
          shb.DumpFilesStandard( log_control, logc, 12 );    }

     // Remove Fosmid vector.

     if ( SAMPLE == "hpool2" || SAMPLE == "hpool3" && heur.REMOVE_FOSMID_VECTOR )
          RemoveFosmidVector( shb, log_control, logc );

     // Make loops and unroll them.

     if ( heur.MAKE_LOOPS && start_step <= 13 )
     {    shb.MakeLoops(logc);
          shb.UnrollLoops(logc);
          shb.DumpFilesStandard( log_control, logc, 13 );    }

     // Delete weakly competing edges.

     if ( start_step <= 14 )
     {    shb.KillWeakExits2( heur, logc );
          if (logc.STATUS_LOGGING)
               cout << Date( ) << ": deleting reverse complement components" << endl;
          shb.DeleteReverseComplementComponents(logc);
          shb.DumpFilesStandard( log_control, logc, 14 );    }

     // Pull apart simple branches.

     if ( start_step <= 15 )
     {    shb.PullApart( min_weight_split, logc );
          shb.DumpFilesStandard( log_control, logc, 15 );    }

     // Clean up steps.

     if ( start_step <= 16 )
     {    shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
          shb.RemoveSmallMostlyAcyclicComponents( logc );
          shb.TestValid(logc);
          shb.DumpFilesStandard( log_control, logc, 16 );    }

     // Remove weakly supported loops.

     if ( heur.WEAK_LOOPS && start_step <= 17 )
     {    double wclock = WallClockTime( );
          if (logc.STATUS_LOGGING)
               cout << Date( ) << ": removing weakly supported loops" << endl;
          vec<vec<pair<int, int> > > paths_index( shb.EdgeObjectCount( ) );
          for ( int i = 0; i < shb.Paths( ).isize( ); i++ )
          for ( int j = 0; j < shb.Path(i).isize( ); j++ ) 
               paths_index[shb.Path(i,j)].push(i,j);
          const fix64_6 weak = 2;
          const int strong_mult = 10;
          vec<int> dels;
          for ( int v = 0; v < shb.N( ); v++ )
          {    if ( shb.To(v).size( ) != 2 || shb.From(v).size( ) != 2 ) continue;
               int i1 = 0, i2 = 1, j1 = 0, j2 = 1;
               if ( shb.To(v)[i1] != v ) swap(i1, i2);
               if ( shb.To(v)[i1] != v ) continue;
               if ( shb.From(v)[j1] != v ) swap(j1, j2);
               if ( shb.From(v)[j1] != v ) continue;
               int u = shb.To(v)[i2], w = shb.From(v)[j2];
               if ( u == v || u == w || v == w ) continue;
               int e = shb.EdgeObjectIndexByIndexTo( v, i2 );
               int f = shb.EdgeObjectIndexByIndexTo( v, i1 );
               int g = shb.EdgeObjectIndexByIndexFrom( v, j2 );
               fix64_6 eg_mult = 0, ef_mult = 0, fg_mult = 0;
               for ( int i = 0; i < paths_index[f].isize( ); i++ )
               {    int id = paths_index[f][i].first, p = paths_index[f][i].second;
                    if ( p > 0 && shb.Path(id)[p-1] == e ) ef_mult += shb.Weight(id);
                    if ( p < shb.Path(id).isize( ) - 1 && shb.Path(id)[p+1] == g )
                         fg_mult += shb.Weight(id);     }
               for ( int i = 0; i < paths_index[g].isize( ); i++ )
               {    int id = paths_index[g][i].first, p = paths_index[g][i].second;
                    if ( p > 0 && shb.Path(id)[p-1] == e ) 
                         eg_mult += shb.Weight(id);    }
               if ( eg_mult == 0 ) continue;
               if ( logc.verb[ "WEAK_LOOPS" ] >= 1 )
                    PRINT6( e, f, g, eg_mult, ef_mult, fg_mult );
               if ( ( ef_mult <= weak && eg_mult >= strong_mult * ef_mult )
                    || ( fg_mult <= weak && eg_mult >= strong_mult * fg_mult ) )
               {    dels.push_back(f);    }    }
          UniqueSort(dels);
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": deleting " << dels.size( ) << " weak loops" 
                    << endl;    }
          shb.DeleteEdges(dels);
          shb.RemoveDeadEdgeObjects( );
          shb.RemoveEdgelessVertices( );
          shb.TestValid(logc);
          REPORT_TIME( wclock, "used removing weakly supported loops" );
          shb.DumpFilesStandard( log_control, logc, 17 );    }

     // More pulling apart.

     if ( heur.PULL_APART2 && start_step <= 18 )
     {    shb.PullApart2( min_weight_split, logc );
          shb.DumpFilesStandard( log_control, logc, 18 );    }

     // More weak exit killing.

     if ( start_step <= 19 )
     {    shb.KillWeakExits2( heur, logc );
          shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
          shb.RemoveSmallMostlyAcyclicComponents(logc);
          shb.DumpFilesStandard( log_control, logc, 19 );    }

     // Divine bubbles.

     if ( start_step <= 20 && bases.size( ) > 0 )
     {    shb.DivineBubbles( heur.DIVINE_K, bases, quals, heur, logc );
          shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
          shb.RemoveSmallMostlyAcyclicComponents(logc);
          shb.RemoveEdgelessVertices( );
          if ( heur.PACBIO_PATCH && ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) )
          {    vecbasevector pb;
               pb.ReadAll( TMP + "/pb.fastb" );
               UnrollWithPacbioReads( pb, &shb, useOldLRPMethod, heur.CORRECT_PACBIO_PATCH );
               if ( heur.REMOVE_FOSMID_VECTOR )
                    RemoveFosmidVector( shb, log_control, logc, False );    }
          shb.DumpFilesStandard( log_control, logc, 20 );    }

     // Gulp edges.

     if ( heur.GULP && start_step <= 21 )
     {    shb.Gulp(log_control,logc);
          shb.Ungulp(logc);
          shb.DeleteReverseComplementComponents(logc);
          shb.DumpFilesStandard( log_control, logc, 21 );    }

     // Orient to reference.

     if ( start_step <= 22 && heur.ORIENT_TO_REFERENCE && log_control.G != 0 )
     {    OrientToReference( shb, *log_control.G, logc );
          shb.DumpFilesStandard( log_control, logc, 22 );    }

     // Hookup.

     if ( heur.HOOKUP && start_step <= 23 )
     {    shb.Hookup(logc);
          shb.DumpFilesStandard( log_control, logc, 23 );    }
     
     if ( heur.DIVINE_SINGLE_MUTATION && start_step <= 24 && bases.size( ) > 0 )
     {    shb.DivineSingleMutationBubbles( heur.DIVINE_K, bases, quals, heur, logc );
     /*
          shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
          shb.RemoveSmallMostlyAcyclicComponents(logc);
          */
          shb.RemoveEdgelessVertices( );
          shb.DumpFilesStandard( log_control, logc, 24 );    }
}
