///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Equiv.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/KmerBasket.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"

void KmerBasket( const vecbasevector& R, const vecqualvector& Q,      
     HyperBasevector& hb, const int K, const int max_mismatches, 
     const int verbosity, ostream& out )
{
     // Define constants.

     const int M = 12;
     const int min_floor = 9;
     const int min_depth = 2;

     // Sanity checks.

     ForceAssertLe( M, 16 );

     // Define the kmers.

     double clock = WallClockTime( );
     int KP = K + 1;
     vec<basevector> kmers;
     vec<qualvector> qmers;
     for ( size_t i = 0; i < R.size( ); i++ )
     {    for ( int j = 0; j <= R[i].isize( ) - KP; j++ )
          {    basevector b;
               b.SetToSubOf( R[i], j, KP );
               kmers.push_back(b);
               qualvector q;
               q.SetToSubOf( Q[i], j, KP );
               qmers.push_back(q);    }    }
     if ( verbosity >= 1 )
          out << TimeSince(clock) << " used to define kmers" << endl;

     // Build equivalence relation.

     equiv_rel e( kmers.size( ) );
     if ( verbosity >= 1 )
     {    out << "\n";
          PRINT_TO( out, e.OrbitCount( ) );    }
     clock = WallClockTime( );

     // Slower, more sensitive:

     if ( M == 0 )
     {    for ( size_t i1 = 0; i1 < kmers.size( ); i1++ )
          {    for ( size_t i2 = i1 + 1; i2 < kmers.size( ); i2++ )
               {    const basevector &b1 = kmers[i1], &b2 = kmers[i2];
                    int mis = 0;
                    for ( int j = 0; j < KP; j++ )
                    {    if ( b1[j] != b2[j] )
                         {    ++mis;
                              if ( mis > max_mismatches ) break;    }    }
                    if ( mis <= max_mismatches ) e.Join( i1, i2 );    }    }    }

     // Faster, less sensitive, seeds on M-mers:

     if ( M > 0 )
     {    for ( int r = 0; r <= K - M; r++ )
          {    vec<int32_t> mmers( kmers.size( ) );
               for ( size_t j = 0; j < kmers.size( ); j++ )
                    mmers[j] = kmers[j].extractKmer( r, M );
               vec<int> ids( kmers.size( ), vec<int>::IDENTITY );
               SortSync( mmers, ids );
               for ( int l1 = 0; l1 < mmers.isize( ); l1++ )
               {    int l2;
                    for ( l2 = l1 + 1; l2 < mmers.isize( ); l2++ )
                         if ( mmers[l2] != mmers[l1] ) break;
                    for ( int z1 = l1; z1 < l2; z1++ )
                    {    for ( int z2 = l1 + 1; z2 < l2; z2++ )
                         {    int i1 = ids[z1], i2 = ids[z2];
                              if ( e.Equiv( i1, i2 ) ) continue;
                              const basevector &b1 = kmers[i1], &b2 = kmers[i2];
                              int mis = 0;
                              for ( int j = 0; j < KP; j++ )
                              {    if ( b1[j] != b2[j] )
                                   {    ++mis;
                                        if ( mis > max_mismatches ) break;    }    }
                              if ( mis <= max_mismatches ) 
                                   e.Join( i1, i2 );    }    }
                    l1 = l2 - 1;    }    }    }

     if ( verbosity >= 1 )
     {    out << TimeSince(clock) << " used to find equivalence relation" << endl;
          PRINT_TO( out, e.OrbitCount( ) );    }
     clock = WallClockTime( );

     // Find orbits and consensus.

     if ( verbosity >= 2 ) out << "\n" << String( KP, '=' ) << "\n\n";
     vec<int> reps;
     e.OrbitReps(reps);
     vecbasevector unibases;
     for ( int i = 0; i < reps.isize( ); i++ )
     {    vec<int> o;
          e.Orbit( reps[i], o );
          if ( o.isize( ) < min_depth ) continue;
          if ( verbosity >= 2 )
          {    out << "orbit " << i << "\n";
               for ( int j = 0; j < o.isize( ); j++ )
                    out << kmers[ o[j] ].ToString( ) << "\n";
               out << String( KP, '-' ) << "\n";    }
          basevector con(KP);
          qualvector edge(KP);
          int floor = 1000000000;
          int depth_floor = 1000000000;
          for ( int u = 0; u < KP; u++ )
          {    vec<int> score( 4, 0 ), ids( 4, vec<int>::IDENTITY );
               for ( int j = 0; j < o.isize( ); j++ )
               {    int x = o[j];
                    score[ kmers[x][u] ] += qmers[x][u];    }
               ReverseSortSync( score, ids );
               con.Set( u, ids[0] );
               edge[u] = Min( 99, score[0] - score[1] );
               floor = Min( floor, (int) edge[u] );    }
          if ( verbosity >= 2 )
          {    PrintReadWithScores( con, edge, out, KP );
               out << "floor = " << floor;    }
          if ( floor >= min_floor ) 
          {    unibases.push_back(con);
               if ( verbosity >= 2 ) out << ", ACCEPT!";    }
          if ( verbosity >= 2 ) out << "\n\n" << String( KP, '=' ) << "\n\n";    }

     // Generate HyperKmerPath.

     vecbasevector unibases_plus(unibases);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    basevector b = unibases[i];
          b.ReverseComplement( );
          unibases_plus.push_back(b);    }
     vecKmerPath paths, pathsrc, unipaths;
     vec<tagged_rpint> pathsdb, unipathsdb;
     ReadsToPathsCoreY( unibases_plus, K, paths, pathsrc, pathsdb );
     Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
     digraph A;
     BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, unipathsdb, A );
     KmerBaseBroker( K, paths, pathsrc, pathsdb, unibases_plus );
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
     KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, unibases_plus );

     // Keep only the forward edges of the HyperKmerPath h.

     vec<basevector> unibasesv;
     for ( size_t i = 0; i < unibases.size( ); i++ )
          unibasesv.push_back( unibases[i] );
     Sort(unibasesv);
     vec<int> to_delete;
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
     {    Bool good = False;
          basevector e = kbb.Seq( h.EdgeObject(i) );
          for ( int j = 0; j <= e.isize( ) - KP; j++ )
          {    basevector b; 
               b.SetToSubOf( e, j, KP );
               if ( BinMember( unibasesv, b ) ) good = True;    }
          if ( !good ) to_delete.push_back(i);    }
     h.DeleteEdges(to_delete);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.Reverse( );

     // Finish up.

     hb.Initialize( h, kbb );
     if ( verbosity >= 1 )
          out << TimeSince(clock) << " used for the rest" << endl;    }
