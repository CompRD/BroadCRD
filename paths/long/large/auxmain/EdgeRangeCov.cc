///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// EdgeRangeCov.  Given two edges e1 and e2, compute the coverage for the part of
// the graph between e1 and e2.  Note that bad things will happen if you incorrectly
// specify e1 and e2.  Also we do not allow cycles between e1 and e2.

#include "Intvector.h"
#include "MainTools.h"
#include "Set.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", 
          "looks for DIR/a.{hbv,inv,paths.inv}");
     CommandArgument_String_Doc(E, 
          "expression of form e1..e2, or comma-separated list of such" );
     CommandArgument_Bool_OrDefault_Doc(VERBOSE, False, "say what we're doing");
     EndCommandArguments;

     // Load data.

     if (VERBOSE) cout << Date( ) << ": loading data" << endl;
     HyperBasevector hb;
     BinaryReader::readFile( DIR + "/a.hbv", &hb );
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Go through the pairs.

     vec<String> es;
     ParseStringSet( "{" + E + "}", es );
     vec<vec<String>> rows;
     vec<String> row = { "edges", "length", "+/-", "pairs", "cov", "dev" };
     rows.push_back(row);
     for ( int ei = 0; ei < es.isize( ); ei++ )
     {    int e1 = es[ei].Before( ".." ).Int( ), e2 = es[ei].After( ".." ).Int( );

          // Find edges between e1 and e2, and get their inverses too.

          if (VERBOSE) cout << Date( ) << ": finding edges between" << endl;
          vec<int> edges;
          if ( e1 == e2 ) edges.push_back(e1);
          else edges = hb.EdgesBoundedBy( e1, e2, to_left, to_right );
          vec<int> edges2(edges);
          int fails = 0;
          for ( int i = 0; i < edges.isize( ); i++ )
          {    if ( inv[ edges[i] ] >= 0 ) edges2.push_back( inv[ edges[i] ] );
               else fails++;    }
          if ( fails > 0 ) 
          {    cout << "Warning: failed to find inverse of " << fails << " edges." 
                    << endl;    }
          UniqueSort(edges2);

          // Find the associated vertices.

          vec<int> verts;
          for ( int i = 0; i < edges.isize( ); i++ )
               verts.push_back( to_left[ edges[i] ], to_right[ edges[i] ] );
          UniqueSort(verts);
          int v1 = BinPosition( verts, to_right[e1] ); 
          int v2 = BinPosition( verts, to_left[e2] );

          // Compute length of region.
     
          int low, high;
          if ( e1 == e2 ) low = high = hb.EdgeLengthBases(e1);
          else
          {    digraphE<basevector> GX( 
                    digraphE<basevector>::COMPLETE_SUBGRAPH, hb, verts );
               if ( !GX.Acyclic( ) )
               {    cout << "Specified part of graph has cycle." << endl;
                    Scram(1);    }
               vec<int> L;
               for ( int i = 0; i < GX.EdgeObjectCount( ); i++ )
                    L.push_back( GX.EdgeObject(i).isize( ) - hb.K( ) + 1 );
               digraphE<int> G( GX, L );
               if (VERBOSE) 
                    cout << Date( ) << ": computing length of region" << endl;
               vec<int> p;
               G.ShortestPath( v1, v2, p );
               low = hb.EdgeLengthBases(e1) + hb.EdgeLengthKmers(e2);
               for ( int j = 0; j < p.isize( ) - 1; j++ )
               {    int v1 = p[j], v2 = p[j+1];
                    int m = 1000000000;
                    for ( int l = 0; l < G.From(v1).isize( ); l++ )
                    {    if ( G.From(v1)[l] == v2 )
                              m = Min( m, G.EdgeObjectByIndexFrom( v1, l ) );    }
                    low += m;    }
               for ( int i = 0; i < G.EdgeObjectCount( ); i++ )
                    G.EdgeObjectMutable(i) = -G.EdgeObject(i);
               high = hb.EdgeLengthBases(e1) + hb.EdgeLengthKmers(e2);
               for ( int j = 0; j < p.isize( ) - 1; j++ )
               {    int v1 = p[j], v2 = p[j+1];
                    int m = 1000000000;
                    for ( int l = 0; l < G.From(v1).isize( ); l++ )
                    {    if ( G.From(v1)[l] == v2 )
                              m = Min( m, G.EdgeObjectByIndexFrom( v1, l ) );    }
                    high -= m;    }    }

          // Compute coverage.

          if (VERBOSE) cout << Date( ) << ": loading pid info" << endl;
          VecULongVec P;
          P.Read( DIR + "/a.paths.inv", edges2 );
          vec<int> pids;
          for ( int i = 0; i < (int) P.size( ); i++ )
          for ( int j = 0; j < (int) P[i].size( ); j++ )
               pids.push_back( P[i][j] );
          UniqueSort(pids);
          if (VERBOSE) cout << Date( ) << ": computing coverage" << endl;
          int npids = pids.size( );
          double cov = double( npids ) / double( (low+high)/2 );
          double dev = cov / sqrt(npids);
          vec<String> row = { es[ei], ToString( (low+high)/2 ),
               ToString( (high-low)/2 ),
               ToString(npids), ToString(cov,4), ToString(dev,4) };
          rows.push_back(row);    }

     // Print report.

     PrintTabular( cout, rows, 2 );
     Scram(0);    }
