/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: ScaffoldComponents

   Take connected components of an assembly, and scaffold them into
   scaffolds.
*/

#include "Basevector.h"
#include "MainTools.h"
#include "ReadPairing.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/PdfEntry.h"
#include "paths/SeqOnHyper.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "feudal/BinaryStream.h"

void ReadsHittingUniqueKmers2( const HyperKmerPath& h, const HyperBasevector& hb,
			      const vec<tagged_rpint>& uniqdb20,
			      const vec<CompressedSeqOnHyper>& csaligns,
			      vec<Bool>& hits_unique ) {

  // Mark (putatively) unique kmers on edges.
  
  int nedges = h.EdgeObjectCount( );
  vec<nkmers_t> L(nedges); 
  for ( int i = 0; i < nedges; i++ )
    L[i] = h.EdgeLength(i);

  vec< vec<ho_interval> > huniq(nedges);
  for ( int i = 0; i < nedges; i++ ) {
    vec<ho_interval> huniq0;
    const KmerPath& e = h.EdgeObject(i);
    int pos = 0;
    for ( int j = 0; j < e.NSegments( ); j++ ) {
      vec<path_interval_id_t> places;
      const KmerPathInterval& J = e.Segment(j);
      Contains( uniqdb20, J, places );
      for ( int u = 0; u < places.isize( ); u++ ) {
	const tagged_rpint& t = uniqdb20[ places[u] ];
	int start = Max( 0, int( t.Start( ) - J.Start( ) ) );
	int stop = Min( J.Length( ), int( t.Stop( ) - J.Start( ) ) );
	huniq0.push_back( ho_interval( pos+start, pos+stop ) );
      }
      pos += J.Length( );
    }
    ExtractGivenCoverage( L[i], 1, huniq0, huniq[i] );
  }

  // Determine which read alignments hit unique kmers.
  
  hits_unique.resize_and_set( csaligns.size( ), False );
  for ( int i = 0; i < csaligns.isize( ); i++ ) {
    SeqOnHyper ap;
    csaligns[i].DecompressInto(ap, hb);
    for ( int j = 0; j < ap.N( ); j++ ) {
      int id2 = ap.Id2(j);
      ho_interval H( ap.pos2(j), ap.Pos2(j) - h.K( ) + 1 );
      if ( Overlap( H, huniq[id2] ) > 0 )
	hits_unique[i] = True;
    }
  }
}

/* Function: FindLinksBetweenComponents2
   Find links between components from read pair information
   Does NOT do anything to the HyperKmerPath, just displays a list of links
   (This now attempts to build scaffolds.)
*/

void FindLinksBetweenComponents2(const HyperKmerPath& h, const HyperBasevector& hb,
     const vec<tagged_rpint>& uniqdb20, const vec<CompressedSeqOnHyper>& csaligns,
     const vec<read_pairing>& pairs, size_t nreads)
{
     double clock = WallClockTime( );
     cout << "\nLinks between components:\n";
     equiv_rel e;
     h.ComponentRelation(e);
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     vec< vec<int> > saligns_index(nreads);
     for ( int i = 0; i < csaligns.isize( ); i++ )
          saligns_index[ csaligns[i].Id1( ) ].push_back(i);
     vec<int> to_right_vertex;
     h.ToRight(to_right_vertex);
     vec<Bool> hits_unique;
     ReadsHittingUniqueKmers2( h, hb, uniqdb20, csaligns, hits_unique );
     int n = reps.size( );

     // Define a graph G having 2n vertices 0rc, 0fw, 1rc, 1fw, ... where 0, 1, ...
     // are the components of the HyperKmerPath h.  We add edges as we find links.

     // Local type: component_idx_t
     // The vertex id of a component in the 
     typedef int component_idx_t;

     vec< vec<int> > from(2*n), to(2*n), from_edge_obj(2*n), to_edge_obj(2*n);
     vec<int> edges;
     digraphE<int> G( from, to, edges, to_edge_obj, from_edge_obj );

     // Go through the read pairs in three passes, using longer and longer links.

     for ( int pass = 1; pass <= 3; pass++ )
     {
          // Build links.

          vec< pair<component_idx_t,component_idx_t> > links;
          for ( int i = 0; i < pairs.isize( ); i++ ) 
          {    int id1 = pairs[i].id1, id2 = pairs[i].id2, sep = pairs[i].sep;
               if ( pass == 1 && sep > 1000 ) continue;
               if ( pass == 2 && ( sep <= 1000 || sep > 10000 ) ) continue;
               if ( pass == 3 && sep <= 10000 ) continue;
               if ( !saligns_index[id1].solo( ) ) continue;
               if ( !saligns_index[id2].solo( ) ) continue;
               int i1 = saligns_index[id1][0], i2 = saligns_index[id2][0];
               if ( !hits_unique[i1] || !hits_unique[i2] ) continue;
               const CompressedSeqOnHyper &cs1 = csaligns[i1], &cs2 = csaligns[i2];
               int j1 = cs1.Part0( ).Id2( ), j2 = cs2.Part0( ).Id2( );
               int v1 = to_right_vertex[j1], v2 = to_right_vertex[j2];
               int c1 = BinPosition( reps, e.ClassId(v1) );
               int c2 = BinPosition( reps, e.ClassId(v2) );

               // Record links that go between sources and sinks.

               if ( c1 == c2 ) continue;
               int v = 2*c1 + ( cs1.Rc1( ) ? 1 : 0 );
               int w = 2*c2 + ( cs2.Rc1( ) ? 0 : 1 );
               if ( !G.Sink(v) || !G.Source(w) ) continue;
               if ( !Member( G.From(v), w ) ) links.push( v, w );    }
          Sort(links);

          // Add edge to graph for each link which occurs at least MIN_LINKS times.

          const int MIN_LINKS = 60;
          for ( int i = 0; i < links.isize( ); i++ ) 
          {    int j = links.NextDiff(i);
               if ( j - i >= MIN_LINKS ) 
                    G.AddEdge( links[i].first, links[i].second, j-i );
               i = j - 1;    }    }

     // Print the links.

     vec< pair<int,int> > used_links;
     for ( int x = 0; x < 2*n; x++ )
     {    int v = x;
          if ( !G.Source(v) || G.Sink(v) ) continue;
          vec<int> line, counts;
          line.push_back(v);
          Bool fail = False;
          while(1)
          {    if ( G.Sink(v) ) break;
               if ( !G.From(v).solo( ) )
               {    fail = True;
                    break;    }
               int v0 = v;
               v = G.From(v)[0];
               if ( Member( line, v ) )
               {    fail = True;
                    break;    }
               line.push_back(v);
               counts.push_back( G.EdgeObjectByIndexFrom( v0, 0 ) );    }
          if ( !fail )
          {    cout << "\n";
               for ( int j = 0; j < line.isize( ); j++ )
               {    if ( j > 0 ) 
                    {    cout << " --(" << counts[j-1] << ")--> ";
                         used_links.push( line[j-1], line[j] );    }
                    cout << line[j]/2 << ( line[j] % 2 == 0 ? "rc" : "fw" );    }
               cout << "\n";    }    }
     Sort(used_links);
     for ( int v = 0; v < 2*n; v++ )
     {    for ( int j = 0; j < G.From(v).isize( ); j++ )
          {    int w = G.From(v)[j];
               if ( BinMember( used_links, make_pair(v,w) ) ) continue;
               cout << "\n" << v/2 << ( v/2 % 2 == 0 ? "rc" : "fw" ) 
                    << " --(" << G.EdgeObjectByIndexFrom( v, j ) << ")--> "
                    << w/2 << ( w/2 % 2 == 0 ? "rc" : "fw" ) << "\n";    }    }
     for ( int j = 0; j < n; j++ )
     {    if ( G.Source(2*j) && G.Sink(2*j) && G.Source(2*j+1) && G.Sink(2*j+1) )
               cout << "\n" << j << "\n";    }
     cout << "\n" << TimeSince(clock) << " used finding links\n" << endl;    }



int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String(SUBDIR);
  EndCommandArguments;

  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  HyperKmerPath h;
  HyperBasevector hb;
  vec<tagged_rpint> uniqdb20;
  vec<CompressedSeqOnHyper> csaligns;
  cout << Date() << " Loading reads..." << endl;
  vecbasevector reads( sub_dir + "/for_scaffolding.reads.fastb" );
  cout << Date() << " Loading components as HyperKmerPath..." << endl;
  BinaryReader::readFile( sub_dir + "/for_scaffolding.hyper", &h );
  cout << Date() << " Loading components as HyperBasevector..." << endl;
  BinaryReader::readFile( sub_dir + "/for_scaffolding.hyperbasevector", &hb );
  cout << Date() << " loading uniqdb20..." << endl;
  BinaryReader::readFile( sub_dir + "/for_scaffolding.uniqdb20.pathsdb", &uniqdb20 );
  cout << Date() << " Loading csaligns..." << endl;
  BinaryRead( sub_dir + "/for_scaffolding.csaligns", csaligns, h, hb, reads );
  cout << Date() << " Loading pairs..." << endl;
  vec<read_pairing> pairs;
  ReadPairsFile( run_dir + "/reads.pairto", pairs );
  cout << " Finding links between components..." << endl;
  FindLinksBetweenComponents2( h, hb, uniqdb20, csaligns, pairs, reads.size() );
  
}
