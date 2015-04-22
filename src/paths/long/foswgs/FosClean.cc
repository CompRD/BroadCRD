///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/Heuristics.h"
#include "paths/long/Logging.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/fosmid/Fosmids.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault(INSTANCE, "");
     CommandArgument_Bool_OrDefault(FOS_FILTER, False);
     CommandArgument_String_OrDefault(N, "");
     EndCommandArguments;

     if ( INSTANCE != "" ) INSTANCE += ".";

     cout << Date( ) << ": loading Fosmid reads" << endl;
     vecbasevector bases( 
          "/wga/scr4/jaffe/fos_filter/tmp.fos/frag_reads_orig.fastb" );

     cout << Date( ) << ": loading selected Fosmid reads" << endl;
     vecbasevector bases1;
     vecqualvector quals1;
     if ( N == "all" )
     {    bases1 = bases;
          quals1.ReadAll( 
               "/wga/scr4/jaffe/fos_filter/tmp.fos/frag_reads_orig.fastb" );    }
     else
     {    vec<int> fos;
          // if ( N == "all" ) fos = AllFosmids( );
          // else 
          ParseIntSet( N, fos );
          for ( int i = 0; i < fos.isize( ); i++ )
          {    int n = fos[i];
               bases1.ReadAll( "/wga/scr4/jaffe/CompareVars/" + ToString(n)
                    + "/tmp.fos/frag_reads_orig.fastb", True );
               quals1.ReadAll( "/wga/scr4/jaffe/CompareVars/" + ToString(n)
                    + "/tmp.fos/frag_reads_orig.qualb", True );    }    }

     cout << Date( ) << ": getting assembly" << endl;
     SupportedHyperBasevector shb;
     BinaryReader::readFile( "sss." + INSTANCE + "final.shbv", &shb );
     vecbasevector edges;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          edges.push_back( shb.EdgeObject(e) );

     if (FOS_FILTER)
     {    cout << Date( ) << ": building all" << endl;
          vec<int> marked( edges.size( ), False );
          const int K = 100;
          ForceAssertEq( K, shb.K( ) );
          vecbasevector all(bases);
          all.Append(edges);
          vec< triple<kmer<K>,int,int> > kmers_plus;
          cout << Date( ) << ": making kmer lookup" << endl;
          MakeKmerLookup2( all, kmers_plus );
          cout << Date( ) << ": marking pairs" << endl;
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               Bool valid = False;
               int count = 0;
               for ( int64_t k = i; k < j; k++ )
                    if ( kmers_plus[k].second < (int) bases.size( ) ) count++;
               if ( count >= 2 ) valid = True;
               if (valid)
               {    for ( int64_t k = i; k < j; k++ )
                    {    if ( kmers_plus[k].second >= (int) bases.size( ) ) 
                         {    int id = kmers_plus[k].second - (int) bases.size( );
                              marked[id] = True;    }    }    }
               i = j - 1;    }
          PRINT2( shb.EdgeObjectCount( ), Sum(marked) );
          vec<int> dels;
          for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               if ( !marked[e] ) dels.push_back(e);
          shb.DeleteEdges(dels);
          shb.RemoveUnneededVertices( );
          shb.RemoveDeadEdgeObjects( );
          shb.RemoveEdgelessVertices( );    }

     long_logging logc( "" );

     // Delete short alpha satellite edges.

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

     // Clean up.

     const double junk_ratio = 10.0;
     const int max_del = 1000;
     long_heuristics heur( "" );
     shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
     shb.RemoveSmallMostlyAcyclicComponents(logc);
     const int max_small_comp = 1500;
     shb.RemoveSmallComponents2( logc, max_small_comp );
     shb.DeleteReverseComplementComponents(logc);

     // Map selected Fosmid reads to assembly and use them to kill edges.

     while(1)
     {    cout << Date( ) << ": setting up locs" << endl;
          HyperBasevector hb_fw(shb), hb_rc(shb);
          hb_rc.Reverse( );
          vec<int> to_right_fw, to_right_rc;
          hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
          vecbasevector x_fw, x_rc;
          for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
               x_fw.push_back( hb_fw.EdgeObject(i) );
          for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
               x_rc.push_back( hb_rc.EdgeObject(i) );
          VecIntPairVec locs_fw, locs_rc;
          const int L = 15;
          CreateGlocs( x_fw, L, locs_fw );
          CreateGlocs( x_rc, L, locs_rc );
          const int max_locs = 100;
          for ( int64_t i = 0; i < (int64_t) locs_fw.size( ); i++ )
          {    if ( (int64_t) locs_fw[i].size( ) > max_locs )
                    locs_fw[i].resize(0);    }
          for ( int64_t i = 0; i < (int64_t) locs_rc.size( ); i++ )
          {    if ( (int64_t) locs_rc[i].size( ) > max_locs )
                    locs_rc[i].resize(0);    }
     
          vec<fix64_6> support( shb.EdgeObjectCount( ), 0 );
          cout << Date( ) << ": aligning" << endl;
          const int batch = 10000;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t id0 = 0; id0 < (int64_t) bases1.size( ); id0 += batch )
          {    vec<vec<read_place>> PLACES;
               for ( int64_t id = id0; 
                    id < Min( id0 + batch, (int64_t) bases1.size( ) ); id++ )
               {    vec<read_place> places;
                    int n = KmerId( bases1[id], L, 0 );
                    const int infinity = 1000000000;
                    int qual_sum = infinity;
                    FindPlaces( bases1[id], quals1[id], n, hb_fw, hb_rc, to_right_fw,
                         to_right_rc, locs_fw, locs_rc, places, qual_sum );
                    if ( places.size( ) <= 2 ) PLACES.push_back(places);    }
               #pragma omp critical
               for ( int l = 0; l < PLACES.isize( ); l++ )
               {    int np = PLACES[l].size( );
                    for ( int i = 0; i < np; i++ )
                    {    if ( PLACES[l][i].Qsum( ) > 20000 ) continue;
                         if ( PLACES[l][i].N( ) == 1 ) continue;
                         for ( int j = 0; j < PLACES[l][i].N( ); j++ )
                         {    support[ PLACES[l][i].E(j) ] += fix64_6(1,np);
                              /*
                              cout << "read " << id << " supports edge " 
                                   << PLACES[l][i].E(j) 
                                   << " with weight " << fix64_6(1,np) << endl;    
                              */
                                   }    }    }    }

          vec<int> dels;
          const int min_mult = 5;
          cout << Date( ) << ": identifying edges to delete" << endl;
          for ( int v = 0; v < shb.N( ); v++ )
          {    {    int n = shb.From(v).size( );
                    vec<fix64_6> s( n, 0 );
                    for ( int j = 0; j < n; j++ )
                         s[j] = support[ shb.EdgeObjectIndexByIndexFrom( v, j ) ];
                    vec<int> id( n, vec<int>::IDENTITY );
                    ReverseSortSync( s, id );
                    if ( s.size( ) >= 2 && s[0] >= min_mult 
                         && s[0] >= min_mult * s[1] )
                    {    for ( int j = 1; j < n; j++ )
                         {    int e = shb.EdgeObjectIndexByIndexFrom( v, id[j] );
                              dels.push_back(e);
                              if ( shb.Inv(e) >= 0 ) 
                                   dels.push_back( shb.Inv(e) );    }    }    }
               {    int n = shb.To(v).size( );
                    vec<fix64_6> s( n, 0 );
                    for ( int j = 0; j < n; j++ )
                         s[j] = support[ shb.EdgeObjectIndexByIndexTo( v, j ) ];
                    vec<int> id( n, vec<int>::IDENTITY );
                    ReverseSortSync( s, id );
                    if ( s.size( ) >= 2 && s[0] >= min_mult 
                         && s[0] >= min_mult * s[1] )
                    {    for ( int j = 1; j < n; j++ )
                         {    int e = shb.EdgeObjectIndexByIndexTo( v, id[j] );
                              dels.push_back(e);
                              if ( shb.Inv(e) >= 0 ) 
                                   dels.push_back( shb.Inv(e) );    }    }    }    }
          PRINT( dels.size( ) );
          shb.DeleteEdges(dels);
          shb.DeleteUnusedPaths( );
          shb.RemoveUnneededVertices( );
          shb.RemoveDeadEdgeObjects( );
          shb.RemoveEdgelessVertices( );
          long_logging logc;
          shb.DeleteReverseComplementComponents(logc);

          // Dump support.  Not really right since we've since edited the graph.
     
          if ( dels.empty( ) )
          {    Ofstream( out, "/wga/dev/jaffe/BroadCRD/belch." + INSTANCE
                    + "support" );
               for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
                    out << e << " --> " << support[e] << endl;
               break;    }    }

     // Clean up again.

     shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
     shb.RemoveSmallMostlyAcyclicComponents(logc);
     shb.RemoveSmallComponents2( logc, max_small_comp );
     shb.DeleteReverseComplementComponents(logc);

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

     // Summarize assembly.

     int nkmers = 0;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          nkmers += shb.EdgeLengthKmers(e);
     cout << "edges = " << ToStringAddCommas( shb.EdgeObjectCount( ) ) << endl;
     vec<int> len;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          len.push_back( shb.EdgeLengthKmers(e) );
     Sort(len);
     cout << "N50 = " << N50(len) << endl;
     cout << "components = " << shb.NComponents( ) << endl;
     cout << "kmers = " << ToStringAddCommas(nkmers) << endl;
     BinaryWriter::writeFile( "/wga/dev/jaffe/BroadCRD/belch." + INSTANCE 
          + "shbv", shb );
     Ofstream( out, "/wga/dev/jaffe/BroadCRD/belch." + INSTANCE + "dot" );
     shb.PrintSummaryDOT0w( out, True, False, True );
     Ofstream( bout, "/wga/dev/jaffe/BroadCRD/belch." + INSTANCE + "fasta" );
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          shb.EdgeObject(e).Print( bout, e );
     cout << Date( ) << ": done" << endl;    }
