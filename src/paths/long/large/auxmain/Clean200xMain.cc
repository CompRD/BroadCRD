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
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Clean200.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", 
          "looks for DIR/a.{fastb,inv,paths,paths.inv}");
     CommandArgument_String_Doc(E, "edge to examine");
     EndCommandArguments;

     const int version = 2;

     // Heuristics.

     const int max_del = 15;
     const int max_exts = 10;
     const int npasses = 2;
     const int max_rl = 250;

     // Load assembly.

     HyperBasevectorX hb;
     BinaryReader::readFile( DIR + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );

     // Test edge.

     int e = E.Int( );
     int v = hb.ToRight(e);
     if ( hb.To(v).size( ) == 0 || hb.From(v).size( ) <= 1 )
     {    cout << "edge doesn't meet requirements" << endl;
          Scram(0);    }

     // Find extensions of e.

     int n = hb.From(v).size( );
     int depth = max_rl;
     vec< vec<int> > exts;
     GetExtensions( hb, v, max_exts, exts, depth );
     if ( exts.isize( ) > max_exts ) 
     {    cout << "too many extensions" << endl;
          Scram(0);    }
     int N = exts.size( );
     vec<int> ei(N);
     for ( int i = 0; i < N; i++ )
     for ( int j = 0; j < n; j++ )
          if ( exts[i][0] == hb.IFrom( v, j ) ) ei[i] = j;

     // Convert to basevectors.

     vec<basevector> bexts, rbexts;
     for ( int i = 0; i < N; i++ )
     {    const vec<int>& x = exts[i];
          basevector b = hb.Cat(x);
          bexts.push_back(b);
          b.ReverseComplement( );
          rbexts.push_back(b);    }

     // Look up edges.

     vec<int> ed;
     for ( int u = 0; u < (int) hb.To(v).size( ); u++ )
     {    int e = hb.ITo( v, u );
          int re = inv[e];
          ed.push_back( e, re );    }
     for ( int m = 0; m < n; m++ )
     {    int ep = hb.IFrom( v, m );
          int rep = inv[ hb.IFrom( v, m ) ];
          ed.push_back( ep, rep );    }
     UniqueSort(ed);
     VecULongVec ind;
     ind.Read( DIR + "/a.paths.inv", ed );

     // Look up paths.

     vec<int> ids;
     for ( int j = 0; j < (int) ind.size( ); j++ )
     for ( int i = 0; i < (int) ind[j].size( ); i++ )
     {    int id = ind[j][i];
          ids.push_back(id);    }
     UniqueSort(ids);
     ReadPathVec paths;
     paths.Read( DIR + "/a.paths", ids );
     vecbasevector bases;
     bases.Read( DIR + "/../data/frag_reads_orig.fastb", ids );
     vecqualvector quals;
     quals.Read( DIR + "/../data/frag_reads_orig.qualb", ids );

     // Score each read path containing e or a successor of it.

     vec<vec<int>> scores(n);
     vec< pair<int,int> > pi; // {(id,start)}
     vec<int> es;
     for ( int u = 0; u < (int) hb.To(v).size( ); u++ )
     {    int e = hb.ITo( v, u );
          es.push_back(e);
          int eid = BinPosition( ed, e );
          for ( int i = 0; i < (int) ind[eid].size( ); i++ )
          {    int id = ind[eid][i];
               const ReadPath& p = paths[ BinPosition( ids, id ) ];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == e ) 
                    {    int start = p.getOffset( );
                         for ( int l = 0; l <= j; l++ )
                              start -= hb.Kmers( p[l] );
                         pi.push( id, start );    }    }    }    }
     for ( int m = 0; m < n; m++ )
     {    int ep = hb.IFrom( v, m );
          int eid = BinPosition( ed, ep );
          for ( int i = 0; i < (int) ind[eid].size( ); i++ )
          {    int id = ind[eid][i];
               const ReadPath& p = paths[ BinPosition( ids, id ) ];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == ep ) 
                    {    if ( j > 0 && Member( es, p[j-1] ) ) continue;
                         int start = p.getOffset( );
                         for ( int l = 0; l < j; l++ )
                              start -= hb.Kmers( p[l] );
                         pi.push( id, start );    }    }    }    }
     for ( int i = 0; i < pi.isize( ); i++ )
     {    int id = pi[i].first, start = pi[i].second;
          int idx = BinPosition( ids, id );
          const ReadPath& p = paths[idx];
          vec<int> q( N, 0 );
          for ( int pos = 0; pos < depth + hb.K( ) - 1; pos++ )
          {    int rpos = pos - start;
               if ( rpos < 0 || rpos >= bases[idx].isize( ) ) continue;
               for ( int l = 0; l < N; l++ )
               {    if ( bexts[l][pos] != bases[idx][rpos] )
                         q[l] += quals[idx][rpos];    }    }
          vec<int> qq( n, 1000000000 );
          for ( int l = 0; l < N; l++ )
               qq[ ei[l] ] = Min( qq[ ei[l] ], q[l] );
          vec<int> idy( n, vec<int>::IDENTITY );
          SortSync( qq, idy );
          if ( qq[0] < qq[1] ) scores[ idy[0] ].push_back( qq[1] - qq[0] );    }

     // Score each read path containing rc of e or its successor.

     vec< pair<int,int> > rpi; // {(id, start of read rel edge re)}
     vec<int> res;
     for ( int u = 0; u < (int) hb.To(v).size( ); u++ )
     {    int e = hb.ITo( v, u );
          int re = inv[e];
          res.push_back(re);
          int eid = BinPosition( ed, re );
          for ( int i = 0; i < (int) ind[eid].size( ); i++ )
          {    int id = ind[eid][i];
               const ReadPath& p = paths[ BinPosition( ids, id ) ];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == re ) 
                    {    int start = p.getOffset( );
                         for ( int l = 0; l < j; l++ )
                              start -= hb.Kmers( p[l] );
                         rpi.push( id, start );    }    }    }    }
     for ( int m = 0; m < n; m++ )
     {    int rep = inv[ hb.IFrom( v, m ) ];
          int eid = BinPosition( ed, rep );
          for ( int i = 0; i < (int) ind[eid].size( ); i++ )
          {    int id = ind[eid][i];
               const ReadPath& p = paths[ BinPosition( ids, id ) ];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == rep ) 
                    {    if ( j < (int) p.size( ) - 1 && Member( res, p[j+1] ) )
                              continue;
                         int start = p.getOffset( );
                         for ( int l = 0; l <= j; l++ )
                              start -= hb.Kmers( p[l] );
                         rpi.push( id, start );    }    }    }    }
     for ( int i = 0; i < rpi.isize( ); i++ )
     {    int id = rpi[i].first, start = rpi[i].second;
          int idx = BinPosition( ids, id );
          const ReadPath& p = paths[idx];
          vec<int> q( N, 0 );
          for ( int pos = 0; pos < depth + hb.K( ) - 1; pos++ )
          {    int rpos = hb.K( ) - 2 - pos - start;
               if ( rpos < 0 || rpos >= bases[idx].isize( ) ) continue;
               for ( int l = 0; l < N; l++ )
               {    int s = bexts[l].size( );
                    if ( rbexts[l][s-pos-1] != bases[idx][rpos] )
                         q[l] += quals[idx][rpos];    }    }
          vec<int> qq( n, 1000000000 );
          for ( int l = 0; l < N; l++ )
               qq[ ei[l] ] = Min( qq[ ei[l] ], q[l] );
          vec<int> idy( n, vec<int>::IDENTITY );
          SortSync( qq, idy );
          if ( qq[0] < qq[1] ) scores[ idy[0] ].push_back( qq[1] - qq[0] );    }

     // Analyze scores.

     for ( int j = 0; j < n; j++ )
          ReverseSort( scores[j] );
     {    vec<int> qsum(n);
          for ( int j = 0; j < n; j++ )
          {    for ( int i = 0; i < scores[j].isize( ); i++ )
                    qsum[j] += scores[j][i];    }
          vec<int> ids( n, vec<int>::IDENTITY );
          ReverseSortSync( qsum, ids );
          for ( int j = 0; j < n; j++ )
          {    cout << "e" << j+1 << ": ";
               cout << printSeq( scores[ ids[j] ] ) << endl;    }
          cout << endl;    }
     vec<int> to_delete;
     const int verbosity = 1;
     AnalyzeScores( hb, inv, v, scores, to_delete, 1, verbosity, version );
     Scram(0);    }
