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

int main( )
{
     RunTime( );

     cout << Date( ) << ": loading Fosmid reads" << endl;
     vecbasevector bases1( 
          "/wga/scr4/jaffe/fos_filter/tmp.fos/frag_reads_orig.fastb" );
     vecqualvector quals1( 
          "/wga/scr4/jaffe/fos_filter/tmp.fos/frag_reads_orig.qualb" );

     cout << Date( ) << ": getting assembly" << endl;
     SupportedHyperBasevector shb;
     BinaryReader::readFile( "beta.shbv", &shb );

     long_logging logc( "" );

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

          if ( dels.empty( ) ) break;    }

     // Clean up again.

     shb.RemoveSmallMostlyAcyclicComponents(logc);
     const int max_small_comp = 1500;
     shb.RemoveSmallComponents2( logc, max_small_comp );
     shb.DeleteReverseComplementComponents(logc);

     PRINT( shb.EdgeObjectCount( ) );

     BinaryWriter::writeFile( "/wga/dev/jaffe/BroadCRD/gamma.shbv", shb );
     cout << Date( ) << ": done" << endl;    }
