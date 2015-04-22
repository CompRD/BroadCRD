///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/Heuristics.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"

int main( )
{    RunTime( );

     // Load shb.

     SupportedHyperBasevector shb;
     BinaryReader::readFile( "gamma.shbv", &shb );

     // Set up logging etc.

     long_logging logc( "" );
     long_heuristics heur( "" );
     ref_data ref;
     vec<ref_loc> readlocs;
     String VERB = "", OUT_INT_HEAD = "";
     long_logging_control log_control( ref, &readlocs, OUT_INT_HEAD, VERB );

     // Define heuristics.

     const double min_weight_split = 10.0;
     const double junk_ratio = 10.0;
     const int max_del = 100;

     // Simplify.

     PRINT( shb.EdgeObjectCount( ) );

     // Delete short alpha satellite edges. (Ends up eliminating 11 edges)

     cout << Date( ) << ": deleting alpha satellite" << endl;
     vec<int> dels2;
     double max_score = 20.0;
     basevector alpha( 
"AGCATTCTCAGAAACTTCTTTGTGATGTGTGTATTCAACTCACAGAGTTGAACATTTCTTTTGATAGAGCAGTTTGG"
"AAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCGCTTTGAGGATTATGGTGGAAAAGGGAATATCTTCA"
"TATAAAAACTAGACAGA" );
     basevector alphan;
     for ( int j = 1; j <= 4; j++ )
          alphan = Cat( alphan, alpha );
     vec<Bool> looks_alpha( shb.EdgeObjectCount( ), False );
     #pragma omp parallel for
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    if ( shb.EdgeObject(e).isize( ) > 400 ) continue;
          basevector E = shb.EdgeObject(e);
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) E.ReverseComplement( );
               int best_loc;
               alignment a;
               SmithWatFree( E, alphan, best_loc, a );
               vec<int> mgg = a.MutationsGap1Gap2( E, alphan );
               double score 
                    = 100.0 * ( mgg[0] + 3 * ( mgg[1] + mgg[2] ) ) / E.isize( );
               if ( score <= max_score )
               {    looks_alpha[e] = True;
                    break;    }    }    }
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          if ( looks_alpha[e] ) dels2.push_back(e);
     shb.DeleteEdges(dels2);
     shb.RemoveUnneededVertices( );
     shb.RemoveDeadEdgeObjects( );
     shb.RemoveEdgelessVertices( );
     cout << Date( ) << ": alpha done" << endl;

     // Delete components having no edge of at least 1000 kmers.

     vec<int> e_to_delete;
     vec< vec<int> > comps;
     shb.Components(comps);
     for ( size_t i = 0; i < comps.size( ); i++ )
     {    const vec<int>& o = comps[i];
          vec<int> e;
          for ( size_t j = 0; j < o.size( ); j++ )
          {    int v = o[j];
               {    for ( size_t t = 0; t < shb.From(v).size( ); t++ )
                    {    e.push_back( shb.EdgeObjectIndexByIndexFrom( 
                              v, t ) );    }    }    }
          int maxe = 0;
          for ( int j = 0; j < e.isize( ); j++ )
               maxe = Max( maxe, shb.EdgeLengthKmers( e[j] ) );
          if ( maxe < 1000 ) e_to_delete.append(e);    }
     shb.DeleteEdges(e_to_delete);
     shb.TruncatePaths(logc);
     shb.RemoveDeadEdgeObjects( );
     shb.RemoveEdgelessVertices( );
     shb.FixWeights(logc);
     shb.TestValid(logc);

     for ( int pass = 1; pass <= 3; pass++ )
     {    shb.PullApart( min_weight_split, logc );
          shb.PullApart2( min_weight_split, logc );
          // shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
          shb.RemoveSmallMostlyAcyclicComponents(logc);
          shb.DeleteLowCoverage( heur, log_control, logc );
          shb.RemoveSmallMostlyAcyclicComponents(logc);
          // shb.KillWeakExits2( heur, logc );
          shb.DeleteReverseComplementComponents(logc);
          // shb.Hookup(logc);
          PRINT( shb.EdgeObjectCount( ) );    }

     shb.WordifyAlt2(log_control,logc);
     shb.Ungulp(logc);

     for ( int pass = 1; pass <= 1; pass++ )
     {    shb.PullApart( min_weight_split, logc );
          shb.PullApart2( min_weight_split, logc );
          // shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
          shb.RemoveSmallMostlyAcyclicComponents(logc);
          shb.DeleteLowCoverage( heur, log_control, logc );
          shb.RemoveSmallMostlyAcyclicComponents(logc);
          // shb.KillWeakExits2( heur, logc );
          shb.DeleteReverseComplementComponents(logc);
          // shb.Hookup(logc);
          PRINT( shb.EdgeObjectCount( ) );    }

     shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );

     // does bad things
     /*
     vecbasevector bases;
     vecqualvector quals;
     String TMP = "/wga/scr4/jaffe/Glumper/tmp.xxx.21";
     bases.ReadAll( TMP + "/frag_reads_orig.fastb" );
     quals.ReadAll( TMP + "/frag_reads_orig.qualb" );
     shb.DivineBubbles( heur.DIVINE_K, bases, quals, heur, logc );
     */

     PRINT( shb.EdgeObjectCount( ) );

     logc.USE_GENOME_FOR_DUMP = False;
     shb.DumpFiles( "delta", log_control, logc );    }
