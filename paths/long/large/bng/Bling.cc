///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Questions:
//
// 1. Can we make alignment of two reads symmetric, so that it does not depend on 
//    the order and orientation of the reads?
//
// 2. What is the right way to simplify the graph?
//
// 3. Can memory requirement be reduced?
//
// 4. Is there a better way to find nearby points in five-space?
//
// 5. How can we simplify a standard bubble?
//
// 6. In the case where a partial restriction site is observed, how often is the
//    site correct?
//
// 7. Can the algorithm be speeded up?
//
// 8. How can contiguity be increased?

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Equiv.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "TokenizeString.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/bng/BngAlign.h"
#include "paths/long/large/bng/BngAlignCore.h"
#include "paths/long/large/bng/BngEdge.h"
#include "paths/long/large/bng/Texps.h"
#include "paths/long/large/bng/Zmatches.h"

template<class T> Bool IsBubbleStart( const digraphE<T>& G, const int v, int& w  )
{    if ( G.From(v).size( ) < 2 || !G.To(v).solo( ) ) return False;
     vec<int> ws = G.From(v);
     for ( int i = 0; i < ws.isize( ); i++ )
     {    while( G.To( ws[i] ).solo( ) && G.From( ws[i] ).solo( ) )
               ws[i] = G.From( ws[i] )[0];    }
     UniqueSort(ws);
     if ( !ws.solo( ) ) return False;
     w = ws[0];
     if ( !G.From(w).solo( ) ) return False;
     if ( G.From(v).size( ) != G.To(w).size( ) ) return False;
     int x = G.To(v)[0], y = G.From(w)[0];
     vec<int> s = { x, y, v, w };
     UniqueSort(s);
     return s.size( ) == 4;    }

// Splay all vertices at branches except those in bubbles that look something like 
// this:
//
//  -------->    ------> -----> ------->  ----------->
//               ----------> ---------->
//               ----> ---------------->

template<class T> void SplayAll( digraphE<T>& G )
{    vec<Bool> bubbly( G.N( ), False );
     for ( int v = 0; v < G.N( ); v++ )
     {    int w;
          if ( IsBubbleStart( G, v, w ) ) bubbly[v] = bubbly[w] = True;    }
     for ( int v = 0; v < G.N( ); v++ )
     {    if ( bubbly[v] ) continue;
          if ( G.From(v).size( ) > 1 || G.To(v).size( ) > 1 )
               G.SplayVertex(v);    }    }

void Compress( const digraphE<bng_edge>& G, digraphE< vec<int> >& I )
{    vec< vec<int> > vedges;
     for ( int e = 0; e < G.E( ); e++ )
     {    vec<int> v = { G.EdgeObject(e).Len( ) };
          vedges.push_back(v);    }
     I.Initialize( G, vedges );
     for ( int v = 0; v < I.N( ); v++ )
     {    if ( !I.To(v).solo( ) || !I.From(v).solo( ) ) continue;
          if ( I.To(v)[0] == v || I.From(v)[0] == v
               || I.To(v)[0] == I.From(v)[0] )
          {    continue;    }
          int e1 = I.ITo( v, 0 ), e2 = I.IFrom( v, 0 );
          vec<int> x = I.EdgeObject(e1);
          x.append( I.EdgeObject(e2) );
          I.JoinEdges( v, x );    }
     I.RemoveDeadEdgeObjects( );
     I.RemoveEdgelessVertices( );    }

int main(int argc, char *argv[])
{
     RunTime( );
     double all_clock = WallClockTime( );

     BeginCommandArguments;
     CommandArgument_Int_OrDefault_Doc(REGION, 1,
          "1 = 10 Mb, 2 = 100 Mb, 3 = all");
     CommandArgument_Int_OrDefault_Doc(NREADS, 1000000,
          "use this number of reads; default is all");
     CommandArgument_Int_OrDefault_Doc(PRINT_FAIL, 0,
          "print this number of failed alignments");
     CommandArgument_Bool_OrDefault_Doc(ONE_SEED, False,
          "use just one seed per read pair");
     CommandArgument_Bool_OrDefault_Doc(SHOW_ASYMMETRY, False,
          "indicate if we found an alignment of read i1 to i2, but "
          "not from i2 to i1");
     CommandArgument_Bool_OrDefault_Doc(ALIGN_LOGGING, False,
          "show info for all alignments");
     CommandArgument_String_Doc(OUT_DIR, "directory for output files");
     CommandArgument_Bool_OrDefault_Doc(OLD, True, "use old data");
     CommandArgument_Double_OrDefault_Doc(MAX_SCORE, 9.5,
          "max score in % for an alignment");
     CommandArgument_Bool_OrDefault_Doc(LINK, False, "find connections");
     CommandArgument_Bool_OrDefault_Doc(START_GRAPH, False, "start with graph");
     CommandArgument_Bool_OrDefault_Doc(SPLAY, False, "splay vertices");
     CommandArgument_Bool_OrDefault_Doc(REMOVE_RC, True, "remove rc components");
     CommandArgument_Bool_OrDefault_Doc(DISCO, False, 
          "align disco and bng assemblies");
     EndCommandArguments;

     // Bool mode = 1; // look for false positives
     Bool mode = 2; // look for true positives
     Mkdir777(OUT_DIR);

     // Heuristics.

     const double max_score = MAX_SCORE/100;   // maximum score for an alignment
     const int min_direct = 10;
     const int min_seed = 9000;          // minimum seed size
     const int min_edges0 = 5;           // min edges in component to keep (initial)
     const int min_edges = 15;           // min edges in component to keep
     const double max_total_err = 0.06;  // max end-to-end distance fraction
     const int seed_radius = 500;
     const int min_seeds = 2;
     const int max_dist_base = 1000;
     const int max_dist_frac = 0.05;
     const int max_diff_zipper = 1000;   // max separation to zipper
     const double max_frac_zipper = 0.5; // max fraction of adjacent to zipper
     const int max_partial_diff = 500;
     const int hang_floor = 3;
     const int hang_mult = 4;
     const int max_vertex = 10;          // max edges incident upon a vertex

     // Algorithm.

     const Bool RFILTER = True;

     // Experiment control.

     const int other = 0;

     // Output control.

     // Bool EDGE_LABEL_MODE = 1; // dot: edges labeled len
     Bool EDGE_LABEL_MODE = 2; // dot: edges labeled len:{r1,...,rn}
     // Bool EDGE_LABEL_MODE = 3; // dot: edges labeled len(l1,...,ln)
     Bool ALIGN_DETAILS = True;
     Bool SHOW_ALL_SEEDS = False;

     // Load BNG map for region.

     vec<vec<double>> X;
     {    cout << "\n" << Date( ) << ": loading map" << endl;

          String fn;
          if ( REGION == 1 ) // map from 10 Mb region
          {    fn = 
               "hcc1143_negativesel_vs_grch38_trimmed_chr4_65000000_75000000.bnx";
                    }
          if ( REGION == 2 ) // map from 100 Mb region
          {    fn =
               "hcc1143_negativesel_vs_grch38_trimmed_chr6_65000000_165000000.bnx";
                    }
          if ( REGION == 3 ) // entire genome
          {    if (OLD) fn = "HCC1143_82x.bnx";
               else fn = "hcc1143AllCondensed.bnx";    }
          fast_ifstream in( "/wga/scr4/vendor/bng/HCC1143/BNX/" + fn );

          String line;
          vec<double> x;
          vec<String> parts;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "1", 0 ) ) continue;
               Tokenize( line, '\t', parts );
               x.clear( );
               for ( int j = 2; j < parts.isize( ); j++ )
                    x.push_back( parts[j].Double( ) - parts[j-1].Double( ) );
               X.push_back(x);
               x.ReverseMe( );
               X.push_back(x);    
               if ( X.isize( )/2 == NREADS ) 
               {    cout << "TRUNCATING" << endl;
                    break;    }    }
          PRINT( X.size( ) );    }
     int nx = X.size( );

     // Calculate coverage.

     int64_t total = 0;
     for ( int i = 0; i < X.isize( ); i++ )
          total += Sum( X[i] );
     int64_t genome_size = 10000000;
     cout << "coverage = " << setprecision(3) << double(total) / genome_size << endl;

     // Load part of full map.

     vec<vec<double>> X2;
     // if ( mode == 1 )
     {    cout << Date( ) << ": loading map" << endl;

          // Full map.
          fast_ifstream in( "/wga/scr4/vendor/bng/HCC1143/BNX/HCC1143_82x.bnx" );

          String line;
          vec<double> x;
          vec<String> parts;
          while(1)
          {    if ( X2.isize( ) == other ) break;
               getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "1", 0 ) ) continue;
               Tokenize( line, '\t', parts );
               x.clear( );
               for ( int j = 2; j < parts.isize( ); j++ )
                    x.push_back( parts[j].Double( ) - parts[j-1].Double( ) );
               X2.push_back(x);    }
          PRINT( X2.size( ) );    }

     {    cout << Date( ) << ": identifying duplicates" << endl;
          vec<Bool> to_delete( X2.isize( ), False );
          #pragma omp parallel for
          for ( int i2 = 0; i2 < X2.isize( ); i2++ )
          for ( int i = 0; i < X.isize( ); i++ )
          {    vec< pair<int,int> > d;
               for ( int j = 0; j < X[i].isize( ); j++ )
                    d.push( X[i][j], 0 );
               for ( int j = 0; j < X2[i2].isize( ); j++ )
                    d.push( X2[i2][j], 1 );
               Sort(d);
               int eq = 0;
               for ( int j = 0; j < d.isize( ) - 1; j++ )
               {    if ( d[j].first == d[j+1].first && d[j].second != d[j+1].second )
                         eq++;    }
               double sim = double(eq) / Max( X[i].size( ), X2[i2].size( ) );
               if ( eq >= 2 && sim >= 0.5 ) to_delete[i2] = True;    }
          int pcount = Sum(to_delete);
          PRINT2( X.size( ), pcount );
          EraseIf( X2, to_delete );
          PRINT( X2.size( ) );
          cout << Date( ) << ": done" << endl;    }
     X.append(X2);

     // Start with all-vs-all alignment, assuming 50 kb seed matching +/- 10%.

     const vec<vec<double>>& Y = ( mode == 1 ? X2 : X );

     vec<int> S;

     /*
     vec<int> SS = { 342,385,392,405 };
     Bool Spartial = True;
     for ( int j = 0; j < SS.isize( ); j++ )
          S.push_back( 2 * (SS[j]/2), ( 2 * (SS[j]/2) ) + 1 );
     */

     for ( int j = 0; j < X.isize( ); j++ )
          S.push_back(j);
     Bool Spartial = False;

     UniqueSort(S);



#include "paths/long/large/bng/SuperScore.h"

/*
auto fsuper = []( int x )
{    if ( x <= 400 ) return 2.5;
     else if ( x > 400 && x <= 1000 ) return 1000.0/x;
     else if ( x > 1000 && x <= 2000 ) return 4 - 3*x/1000.0;
     else return -2.0;    };

// Note asymmetric:

auto gsuper = [&]( int i1, int i2, int j1, int j2 )
{    vec< pair<int,int> > r, l;
     int s1 = 0, s2 = 0;
     for ( int k = j1 + 1; k < X[i1].isize( ); k++ )
     {    s1 += X[i1][k];
          r.push( s1, 1 );    }
     for ( int k = j2 + 1; k < X[i2].isize( ); k++ )
     {    s2 += X[i2][k];
          if ( s2 > s1 ) break;
          r.push( s2, 2 );    }
     s1 = 0, s2 = 0;
     for ( int k = j1 - 1; k >= 0; k-- )
     {    s1 += X[i1][k];
          l.push( s1, 1 );    }
     for ( int k = j2 - 1; k >= 0; k-- )
     {    s2 += X[i2][k];
          if ( s2 > s1 ) break;
          l.push( s2, 2 );    }
     Sort(r), Sort(l);
     double score = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    const vec< pair<int,int> >& x = ( pass == 1 ? r : l );
          for ( int i = 0; i < x.isize( ); i++ )
          {    if ( x[i].second != 2 ) continue;
               int best = 2000;
               for ( int j = i - 1; j >= 0; j-- )
               {    if ( x[i].first - x[j].first >= best ) break;
                    if ( x[j].second != 1 ) continue;
                    best = x[i].first - x[j].first;    }
               for ( int j = i + 1; j < x.isize( ); j++ )
               {    if ( x[j].first - x[i].first >= best ) break;
                    if ( x[j].second != 1 ) continue;
                    best = x[j].first - x[i].first;    }
               score += fsuper(best);    }    }
     return score;    };
*/

     // Build alignment seeds.

     vec< quad<int,int,int,int> > zmatches;
     if ( !START_GRAPH ) 
     {    const int w = 4;
          const int cutoff = 1000;
          Zmatches( X, S, zmatches, True, -1, w, cutoff );    }

     // Generate alignments.

     digraphE<bng_edge> G;
     if (START_GRAPH) BinaryReader::readFile( OUT_DIR + "/1.graph", &G );
     else
     {    cout << Date( ) << ": generating alignments" << endl;
          int align_calls = 0;
          vec< triple< int, int, vec<int> > > aligns;
          vec<double> scores;
          double pclock = WallClockTime( );
          double ac1 = 0, ac2 = 0, ac3 = 0, ac4 = 0, ac5 = 0;

          vec<int64_t> zstarts( X.size( ) + 1, -1 );
          for ( int64_t i = zmatches.jsize( ) - 1; i >= 0; i-- )
               zstarts[ zmatches[i].first ] = i;
          zstarts[ X.size( ) ] = zmatches.size( );
          for ( int j = X.isize( ) - 1; j >= 0; j-- )
               if ( zstarts[j] < 0 ) zstarts[j] = zstarts[j+1];

          #pragma omp parallel for schedule(dynamic, 1)
          for ( int s1 = 0; s1 < S.isize( ); s1++ )
          {    
               int i1 = S[s1];
     
               int64_t low = zstarts[i1], high = zstarts[i1+1];
               if ( low == high ) continue;

               BngAlignCore( X, zmatches, low, high, S, ALIGN_LOGGING, 
                    SHOW_ALL_SEEDS, ALIGN_DETAILS, ONE_SEED, RFILTER, PRINT_FAIL,
                    align_calls, aligns, scores, ac1, ac2, ac3, ac4, ac5,
                    min_direct, max_score, max_total_err, s1, i1 );     }

          UniqueSortSync( aligns, scores );
          cout << "\n";
          int naligns = aligns.size( );
          PRINT2( naligns, align_calls );
          PRINT5( ac1, ac2, ac3, ac4, ac5 );
          cout << "\n" << TimeSince(pclock) << " used in main loop" << endl;
     
          // Check for alignment asymmetry.

          if (SHOW_ASYMMETRY)
          {    cout << "\n";
               vec<vec<int>> a( X.size( ) );
               for ( int i = 0; i < aligns.isize( ); i++ )
                    a[ aligns[i].first ].push_back( aligns[i].second );
               for ( int i = 0; i < aligns.isize( ); i++ )
               {    int id1 = aligns[i].first, id2 = aligns[i].second;
                    if ( !Member( a[id2], id1 ) )
                    {    cout << "found alignment from " << id1 << " to " << id2
                              << ", but not from " << id2 << " to " 
                              << id1 << endl;    }    }    }

          // Symmetrize alignments.  We might try to determine why this is necessary.
          // (See the commented-out block of code below to see the issues.)  The
          // approach here is to delete half the alignments, then when we form the 
          // graph, add two sets of connections for each alignment.

          /*
          for ( int i = 0; i < 10; i++ )
          {    cout << "\n";
               PRINT(i);
               int r = aligns[i].first, s = aligns[i].second;
               const vec<int>& p = aligns[i].third;
               int rx = ( r % 2 == 0 ? r + 1 : r - 1 );
               int sx = ( s % 2 == 0 ? s + 1 : s - 1 );
               cout << "p =";
               for ( int j = 0; j < p.isize( ); j++ )
               {    int rp = p[j] / X[s].size( ), sp = p[j] % X[s].size( );
                    cout << " (" << rp << "," << sp << ")";    }
               cout << "\n";
               for ( int i2 = 0; i2 < aligns.isize( ); i2++ )
               {    int r2 = aligns[i2].first, s2 = aligns[i2].second;
                    if ( r2 != rx || s2 != sx ) continue;
                    const vec<int>& p2 = aligns[i2].third;
                    cout << "p2 =";
                    for ( int j = p2.isize( ) - 1; j >= 0; j-- )
                    {    int rp = p2[j] / X[sx].size( ), sp = p2[j] % X[sx].size( );
                         rp = X[rx].isize( ) - rp - 1;
                         sp = X[sx].isize( ) - sp - 1;
                         cout << " (" << rp << "," << sp << ")";    }
                    cout << "\n";    }    }
          */
          int nal = aligns.size( );
          vec<Bool> adel( nal, False );
          vec< pair<int,int> > aligns12(nal);
          for ( int i = 0; i < nal; i++ )
               aligns12[i] = make_pair( aligns[i].first, aligns[i].second );
          for ( int i = 0; i < nal; i++ )
          {    int r = aligns[i].first, s = aligns[i].second;
               int i2;
               for ( i2 = i + 1; i2 < aligns.isize( ); i2++ )
                    if ( aligns[i2].first != r || aligns[i2].second != s ) break;
               const vec<int>& p = aligns[i].third;
               int rx = ( r % 2 == 0 ? r + 1 : r - 1 );
               int sx = ( s % 2 == 0 ? s + 1 : s - 1 );
               int64_t low = LowerBound( aligns12, make_pair( rx, sx ) );
               int64_t high = UpperBound( aligns12, make_pair( rx, sx ) );
               if ( low < high && make_pair( rx, sx ) < make_pair( r, s ) ) 
               {    for ( int64_t m = i; m < i2; m++ )
                         adel[m] = True;    }
               i = i2 - 1;    }
          EraseIf( aligns, adel );

          // Use the alignments to define an equivalence relation on the reads.

          cout << Date( ) << ": defining equivalence relation" << endl;
          int M = 0;
          for ( int i = 0; i < S.isize( ); i++ )
               M += X[ S[i] ].size( ) + 1;
          equiv_rel e(M);
          PRINT(M);

          // Not right unless S = X....
          vec<int64_t> starts( X.size( ), 0 );
          for ( int i = 0; i < S.isize( ) - 1; i++ )
               starts[i+1] = starts[i] + X[ S[i] ].size( ) + 1;

          vec< pair<int,int> > joins;
          #pragma omp parallel for
          for ( int j = 0; j < aligns.isize( ); j++ )
          {    
               int r = aligns[j].first, s = aligns[j].second;
               const vec<int>& p = aligns[j].third;
     
               // p is an alignment of read r(y1) to read s(y2)

               vec<int> y1, y2;
               for ( int j = 0; j < X[r].isize( ); j++ )
                    y1.push_back( X[r][j] );
               for ( int j = 0; j < X[s].isize( ); j++ )
                    y2.push_back( X[s][j] );

               // Two passes: one for the alignments between the reads, and one for
               // the inferred alignment between their reverses.

               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 )
                    {    r = ( r % 2 == 0 ? r + 1 : r - 1 );
                         s = ( s % 2 == 0 ? s + 1 : s - 1 );    }
     
                    // Get the start positions in e.
     
                    int n1 = starts[r], n2 = starts[s];

                    // Join.  

                    for ( int i = 0; i < p.isize( ) - 1; i++ )
                    {    int r1 = p[i] / y2.size( ), r2 = p[i+1] / y2.size( );
                         int s1 = p[i] % y2.size( ), s2 = p[i+1] % y2.size( );
                         if ( r2 == r1 + 1 && s2 == s1 + 1 )
                         {    int dist1 = Abs( y1[r1] - y2[s1] );
                              int dist2 = Abs( y1[r2] - y2[s2] );
                              double denom1 = Min( y1[r1], y2[s1] );
                              double denom2 = Min( y1[r2], y2[s2] );
                              if ( (dist1 - max_dist_base) / denom1 <= max_dist_frac
                                   && (dist2 - max_dist_base) 
                                        / denom2 <= max_dist_frac )
                              {    
                                   #pragma omp critical
                                   {    if ( pass == 1 )
                                             joins.push( n1 + r2, n2 + s2 );    
                                        else 
                                        {    joins.push( n1 + ( X[r].isize( ) - r2 ),
                                                  n2 + ( X[s].isize( ) - s2 ) );    }
                                             }    }    }    }    }    }
          ParallelUniqueSort(joins);
          cout << Date( ) << ": making joins" << endl;
          for ( int64_t i = 0; i < joins.jsize( ); i++ )
          {    int m1 = joins[i].first, m2 = joins[i].second;
               e.Join( m1, m2 );    }
     
          // Build graph from aligns.
     
          cout << Date( ) << ": building graph" << endl;
          vec<int> o;
          e.OrbitRepsAlt(o);
          int N = o.size( );
          vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
          vec<bng_edge> edges;
          int nn = 0;
          for ( int j = 0; j < S.isize( ); j++ )
          {    int s = S[j];
               vec<int> y;
               for ( int j = 0; j < X[s].isize( ); j++ )
                    y.push_back( X[s][j] );
               for ( int i = 0; i < y.isize( ); i++ )
               {    int v = BinPosition( o, e.ClassId( nn + i ) );
                    int w = BinPosition( o, e.ClassId( nn + i + 1 ) );
                    int m = y[i];
                    from[v].push_back(w), to[w].push_back(v);
                    from_edge_obj[v].push_back( edges.size( ) );
                    to_edge_obj[w].push_back( edges.size( ) );
                    vec< pair<int,int> > p;
                    p.push( s, i );
                    edges.push( m, p );    }
               nn += y.size( ) + 1;    }
          for ( int v = 0; v < N; v++ )
          {    SortSync( from[v], from_edge_obj[v] );
               SortSync( to[v], to_edge_obj[v] );    }
          G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );
          PRINT2( G.N( ), G.E( ) );
          BinaryWriter::writeFile( OUT_DIR + "/1.graph", G );    }

     // Combine parallel edges.

     cout << Date( ) << ": combining parallel edges" << endl;
     vec<int> dels;
     for ( int v = 0; v < G.N( ); v++ )
     {    vec<int> ws = G.From(v);
          vec<int> ids = vec<int>( G.From(v).size( ), vec<int>::IDENTITY );
          SortSync( ws, ids );
          for ( int j = 0; j < ws.isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < ws.isize( ); k++ )
                    if ( ws[k] != ws[j] ) break;
               if ( k - j >= 2 )
               {    vec<int> e;
                    vec< pair<int,int> > f;
                    for ( int l = j; l < k; l++ )
                    {    e.push_back( G.EdgeObject( G.IFrom( v, ids[l] ) ).Len( ) );
                         f.append( 
                              G.EdgeObject( G.IFrom( v, ids[l] ) ).IdPos( ) );    }
                    Sort(f);
                    G.EdgeObjectMutable( G.IFrom( v, ids[j] ) ).SetLen( Mean(e) );
                    G.EdgeObjectMutable( G.IFrom( v, ids[j] ) ).SetIdPos(f);
                    for ( int l = j + 1; l < k; l++ )
                         dels.push_back( G.IFrom( v, ids[l] ) );    }
               else dels.push_back( G.IFrom( v, ids[j] ) );
               j = k - 1;    }    }
     G.DeleteEdges(dels);
     G.RemoveDeadEdgeObjects( );
     G.RemoveEdgelessVertices( );
     BinaryWriter::writeFile( OUT_DIR + "/2.graph", G );

     // Remove rc components, a bit probabilistically.  Slightly dangerous.
     // Two components are deemed rc to each other if they have the same edge
     // lengths and "pair" ids.  Also it would be better to move this up in the
     // code, to right after the graph creation step, but there are
     // computational performance issues.

     if (REMOVE_RC)
     {    vec<vec<int>> compe;
          G.ComponentsE(compe);
          int nc = compe.size( );
          vec<vec<int>> lens(nc), pids(nc);
          for ( int i = 0; i < nc; i++ )
          {    for ( int j = 0; j < compe[i].isize( ); j++ )
               {    int e = compe[i][j];
                    lens[i].push_back( G.EdgeObject(e).Len( ) );
                    for ( int j = 0; j < G.EdgeObject(e).Rids( ).isize( ); j++ )
                         pids[i].push_back( G.EdgeObject(e).Rids( )[j] / 2 );    }
               Sort(lens[i]), Sort(pids[i]);    }
          vec<int> ids( nc, vec<int>::IDENTITY ), dels;
          SortSync( lens, pids, ids );
          for ( int i = 1; i < nc; i++ )
          {    if ( lens[i] == lens[i-1] && pids[i] == pids[i-1] ) 
                    dels.append( compe[ ids[i] ] );    }
          G.DeleteEdges(dels);
          G.RemoveDeadEdgeObjects( );
          G.RemoveEdgelessVertices( );    }
     BinaryWriter::writeFile( OUT_DIR + "/3.graph", G );

     // Explode overloaded vertices.  Don't know how it happens that a vertex
     // can have a huge number of edges incident upon it, but it is bad for the
     // downstream code.  So we kill these vertices.

     cout << Date( ) << ": exploding overloaded vertices" << endl;
     for ( int v = 0; v < G.N( ); v++ )
     {    if ( G.From(v).isize( ) + G.To(v).isize( ) > max_vertex )
               G.SplayVertex(v);    }
     int nover = 0;
     for ( int v = 0; v < G.N( ); v++ )
          nover = Max( nover, G.From(v).isize( ) + G.To(v).isize( ) );
     cout << "now most overloaded vertex has " << nover << " edges" << endl;

     // Remove small components.

     cout << Date( ) << ": removing small components" << endl;
     dels.clear( );
     vec<vec<int>> comp;
     G.ComponentsE(comp);
     for ( int i = 0; i < comp.isize( ); i++ )
          if ( comp[i].isize( ) < min_edges0 ) dels.append( comp[i] );
     G.DeleteEdges(dels);
     G.RemoveDeadEdgeObjects( );
     G.RemoveEdgelessVertices( );
     BinaryWriter::writeFile( OUT_DIR + "/4.graph", G );

     // Zipper up the graph.

     PRINT2( G.N( ), G.E( ) );
     cout << Date( ) << ": zippering" << endl;
     vec<vec<int>> compx;
     G.Components(compx);
     int MC = 0;
     for ( int c = 0; c < compx.isize( ); c++ )
          MC = Max( MC, compx[c].isize( ) );
     cout << "largest component has " << MC << " vertices" << endl;
     const int max_zip_depth = 1000;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int c = 0; c < compx.isize( ); c++ )
     {    while(1)
          {    int zips = 0;
               for ( int pass = 1; pass <= 2; pass++ )
               {    for ( int vi = 0; vi < compx[c].isize( ); vi++ )
                    {    int v = compx[c][vi];
                         int64_t nv = G.From(v).size( ), ntv = G.To(v).size( );
                         if ( nv * nv * ntv > max_zip_depth ) continue;
                         for ( int j1 = 0; j1 < G.From(v).isize( ); j1++ )
                         {    for ( int j2 = j1 + 1; j2 < G.From(v).isize( ); j2++ )
                              {    int w1 = G.From(v)[j1], w2 = G.From(v)[j2];
                                   if ( w1 == v || w2 == v ) continue;
                                   if ( G.HasEdge( w1, w2 ) || G.HasEdge( w2, w1 ) ) 
                                        continue;
                                   int e1 = G.IFrom( v, j1 ), e2 = G.IFrom( v, j2 );
                                   int diff = Abs( G.O(e1).Len( ) - G.O(e2).Len( ) );
                                   if ( diff > max_diff_zipper ) continue;
                                   int min_adj = 1000000000;
                                   for ( int i = 0; i < G.To(v).isize( ); i++ )
                                   {    min_adj = Min( min_adj, 
                                             G.O( G.ITo(v,i) ).Len( ) );    }
                                   for ( int i = 0; i < G.From(w1).isize( ); i++ )
                                   {    min_adj = Min( min_adj, 
                                             G.O( G.IFrom(w1,i) ).Len( ) );    }
                                   for ( int i = 0; i < G.From(w2).isize( ); i++ )
                                   {    min_adj = Min( min_adj, 
                                             G.O( G.IFrom(w2,i) ).Len( ) );    }
                                   if ( diff > max_frac_zipper * min_adj ) continue;
                                   zips++;
                                   vec< pair<int,int> > x = G.O(e1).IdPos( );
                                   x.append( G.O(e2).IdPos( ) );
                                   Sort(x);
                                   if ( x.size( ) > 1000 ) // XXXXXXXXXXXXXXXXXXXXXX
                                        PRINT( x.size( ) ); // XXXXXXXXXXXXXXXXXXXXX
                                   G.EdgeObjectMutable(e1).SetIdPos(x);
                                   vec<int> len( x.size( ) );
                                   for ( int i = 0; i < x.isize( ); i++ )
                                   {    int rid = x[i].first, rpos = x[i].second;
                                        len[i] = X[rid][rpos];    }
                                   G.EdgeObjectMutable(e1).SetLen( Mean(len) );
                                   G.DeleteEdgeFrom( v, j2 );
                                   if ( w1 != w2 ) G.TransferEdges( w2, w1 );    
                                   goto next_vertex;    }    }    
                         next_vertex: continue;    }
     
                    // More zippering.  Note embedded heuristics.
     
                    for ( int vi = 0; vi < compx[c].isize( ); vi++ )
                    {    int v = compx[c][vi];
                         for ( int j1 = 0; j1 < G.From(v).isize( ); j1++ )
                         {    int w1 = G.From(v)[j1], f1 = G.IFrom( v, j1 );
                              if ( w1 == v ) continue;
                              int64_t nv = G.From(v).size( );
                              int64_t nw1 = G.From(w1).size( );
                              if ( nv * nv * nw1 > max_zip_depth ) continue;
                              for ( int m = 0; m < G.From(w1).isize( ); m++ )
                              {    int w2 = G.From(w1)[m], f2 = G.IFrom( w1, m );
                                   if ( w2 == w1 || w2 == v ) continue;
                                   for ( int j2 = 0; j2 < G.From(v).isize( ); j2++ )
                                   {    int x = G.From(v)[j2], e = G.IFrom( v, j2 );
                                        if ( x == v || x == w1 || x == w2 ) continue;
                                        int delta = G.O(e).Len( ) 
                                             - G.O(f1).Len( ) - G.O(f2).Len( );
                                        if ( Abs(delta) > 500 ) continue;
                                        zips++;
                                        G.TransferEdges( x, w2 );    
                                        goto next_vertex2;    }    }    }
                         next_vertex2: continue;    }

                    // Reverse component.

                    const vec<int>& o = compx[c];
                    for ( int j = 0; j < o.isize( ); j++ )
                         G.ReverseVertex( o[j] );    }
               
               if ( zips == 0 ) break;    }    }
     G.RemoveDeadEdgeObjects( );
     G.RemoveEdgelessVertices( );

     // Delete some hanging ends.

     {    vec<int> dels;
          for ( int pass = 1; pass <= 2; pass++ )
          {    for ( int v = 0; v < G.N( ); v++ )
               for ( int j1 = 0; j1 < G.From(v).isize( ); j1++ )
               for ( int j2 = 0; j2 < G.From(v).isize( ); j2++ )
               {    if ( j1 == j2 ) continue;
                    int e1 = G.IFrom( v, j1 ), e2 = G.IFrom( v, j2 );
                    int w1 = G.From(v)[j1], w2 = G.From(v)[j2];
                    if ( G.From(w2).nonempty( ) || !G.To(w2).solo( ) ) continue;
                    if ( G.O(e2).Weight( ) > hang_floor ) continue;
                    if ( G.O(e1).Weight( ) < hang_mult * G.O(e2).Weight( ) ) 
                         continue;
                    dels.push_back(e2);    }
               G.DeleteEdges(dels);
               G.Reverse( );    }
          G.RemoveDeadEdgeObjects( );
          G.RemoveEdgelessVertices( );    }

     // Delete some more hanging ends.  Not clear that this correctly handles cycles.

     const int max_dist = 20;
     vec<int> edgelens( G.E( ), 1 ), tod;
     for ( int pass = 1; pass <= 2; pass++ )
     {    vec<int> D;
          DistancesToEndArr( G, edgelens, max_dist, True, D );
          for ( int v = 0; v < G.N( ); v++ )
          for ( int j1 = 0; j1 < G.From(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < G.From(v).isize( ); j2++ )
          {    int w1 = G.From(v)[j1], w2 = G.From(v)[j2];
               int d1 = 1 + D[w1], d2 = 1 + D[w2];
               if ( d1 >= 5 * d2 && d2 <= 2 ) tod.push_back( G.IFrom( v, j2 ) );    }
          G.Reverse( );    }
     G.DeleteEdges(tod);
     G.RemoveDeadEdgeObjects( );
     G.RemoveEdgelessVertices( );

     // Resolve some partial restriction sites.  Note embedded heuristics.
     // Note that this may be too aggressive, and thus not properly reflecting
     // true polymorphic sites.

     cout << Date( ) << ": removing some partial restriction sites" << endl;
     while(1)
     {    vec<int> dels;
          for ( int v = 0; v < G.N( ); v++ )
          for ( int j1 = 0; j1 < G.From(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < G.From(v).isize( ); j2++ )
          {    if ( j1 == j2 ) continue;
               int w = G.From(v)[j1], x = G.From(v)[j2];
               if ( !G.From(x).solo( ) || G.From(x)[0] != w ) continue;
               int e = G.IFrom( v, j1 );
               int f1 = G.IFrom( v, j2 ), f2 = G.IFrom( x, 0 );
               int le = G.O(e).Len( ), lf1 = G.O(f1).Len( ), lf2 = G.O(f2).Len( );
               int diff = Abs( le - lf1 - lf2 );
               int min_flank = 1000000000;
               for ( int j = 0; j < G.To(v).isize( ); j++ )
                    min_flank = Min( min_flank, G.O( G.ITo(v,j) ).Len( ) );
               for ( int j = 0; j < G.From(w).isize( ); j++ )
                    min_flank = Min( min_flank, G.O( G.IFrom(w,j) ).Len( ) );
               if ( diff > max_diff_zipper || diff > min_flank * max_frac_zipper )
                    continue;
               int ne = G.O(e).Weight( );
               int nf = Max( G.O(f1).Weight( ), G.O(f2).Weight( ) );

               // Delete non-restriction edges.

               // if ( ne == 2 && nf >= 4 ) dels.push_back(e);
               // More aggressive:
               if ( nf >= 3 || ne == 2 ) dels.push_back(e);

               // Delete restriction edges.

               if ( ne >= 4 && nf == 2 && G.To(x).solo( ) ) 
                    dels.push_back( f1, f2 );    }
          if ( dels.empty( ) ) break;
          G.DeleteEdges(dels);
          G.RemoveDeadEdgeObjects( );
          G.RemoveEdgelessVertices( );    }

     // Splay vertices.

     if (SPLAY) SplayAll(G);

     // Compute component links.

     if (LINK)
     {    
          // Heuristics.

          const double max_max_over = 2.0;
          const int max_dist_link = 8;
          const int min_link_len = 5;
          const int link_passes = 3;

          // Do the linking.

          for ( int lpass = 1; lpass <= link_passes; lpass++ )
          {

          cout << "\nLINKING PASS " << lpass << endl;

          cout << "\n";
          vec<vec<int>> comp;
          G.Components(comp);
          vec<int> to_comp( G.N( ) );
          for ( int i = 0; i < comp.isize( ); i++ )
          for ( int j = 0; j < comp[i].isize( ); j++ )
               to_comp[ comp[i][j] ] = i;
          int ncomps = comp.size( );
          vec<int> sources, sinks;
          G.Sources(sources), G.Sinks(sinks);
          int nsources = sources.size( ), nsinks = sinks.size( );
          vec<int> source_to_comp(nsources), sink_to_comp(nsinks);
          vec<int> comp_to_source(ncomps), comp_to_sink(ncomps);
          for ( int i = 0; i < nsources; i++ )
          {    int c = to_comp[ sources[i] ];
               source_to_comp[i] = c;
               comp_to_source[c] = i;    }
          for ( int i = 0; i < nsinks; i++ )
          {    int c = to_comp[ sinks[i] ];
               sink_to_comp[i] = c;
               comp_to_sink[c] = i;    }

          // Compute component lengths.  We assume splaying.

          cout << Date( ) << ": computing component lengths" << endl;
          vec<int> complen( comp.size( ), -1 );
          for ( int c = 0; c < comp.isize( ); c++ )
          {    
               digraphE<bng_edge> M = G.Subgraph( comp[c] );
               vec<int> sources, sinks;
               M.Sources(sources), M.Sinks(sinks);
               Bool done = False;
               if ( sources.solo( ) && sinks.solo( ) )
               {    vec<vec<int>> paths;
                    // Don't know if 10000 cap needed at all.
                    Bool ok = M.EdgePaths( 
                         sources[0], sinks[0], paths, -1, -1, 10000 );
                    if (ok)
                    {    int M = 0;
                         for ( int j = 0; j < paths.isize( ); j++ )
                              M = Max( M, paths[j].isize( ) );
                         complen[c] = M;
                         done = True;    }    }
               if ( !done ) // assume only simple bubbles, not really right
               {    complen[c] = -1;
                    for ( int j = 0; j < comp[c].isize( ); j++ )
                    {    int v = comp[c][j];
                         complen[c]++;
                         complen[c] 
                              -= Max( 0, G.From(v).isize( ) - 1 );    }    }    }

          cout << Date( ) << ": done" << endl;
          vec< vec< pair<int,int> > > sourcex(nsources), sinkx(nsinks);
          for ( int i = 0; i < sources.isize( ); i++ )
          {    int v = sources[i];
               if ( G.From(v).solo( ) )
               {    int e = G.IFrom(v,0);
                    sourcex[i] = G.O(e).IdPos( );    }    }
          for ( int i = 0; i < sinks.isize( ); i++ )
          {    int v = sinks[i];
               if ( G.To(v).solo( ) )
               {    int e = G.ITo(v,0);
                    sinkx[i] = G.O(e).IdPos( );    
                    for ( int j = 0; j < sinkx[i].isize( ); j++ )
                         sinkx[i][j].second++;    }    }
          vec< quad<int,int,int,int> > s4; //(rid,rpos,source-sink-id,source or sink)
          for ( int i = 0; i < sourcex.isize( ); i++ )
          for ( int j = 0; j < sourcex[i].isize( ); j++ )
               s4.push( sourcex[i][j].first, sourcex[i][j].second, i, 0 );
          for ( int i = 0; i < sinkx.isize( ); i++ )
          for ( int j = 0; j < sinkx[i].isize( ); j++ )
               s4.push( sinkx[i][j].first, sinkx[i][j].second, i, 1 );
          ParallelSort(s4);
          vec< vec< quad<int,int,int,int> > > lfrom( sinks.size( ) ); 
          vec< vec< quad<int,int,int,int> > > lto( sources.size( ) ); 
          for ( int i1 = 0; i1 < s4.isize( ); i1++ )
          {    int rid = s4[i1].first, i2;
               for ( i2 = i1 + 1; i2 < s4.isize( ); i2++ )
                    if ( s4[i2].first != rid ) break;
               for ( int j1 = i1; j1 < i2; j1++ )
               {    int rpos1 = s4[j1].second, ssid1 = s4[j1].third;
                    Bool source1 = ( s4[j1].fourth == 0 );
                    if (source1) continue;
                    for ( int j2 = i1; j2 < i2; j2++ )
                    {    int rpos2 = s4[j2].second, ssid2 = s4[j2].third;
                         Bool source2 = ( s4[j2].fourth == 0 );
                         if ( !source2 ) continue;
                         if ( !( rpos1 <= rpos2 ) ) continue;
                         lfrom[ssid1].push( ssid2, rid, rpos1, rpos2 );
                         lto[ssid2].push( ssid1, rid, rpos1, rpos2 );    }    }
               i1 = i2 - 1;    }
          for ( int i = 0; i < sinks.isize( ); i++ )
               Sort( lfrom[i] );
          for ( int i = 0; i < sources.isize( ); i++ )
               Sort( lto[i] );

          cout << "\nCONNECTIONS:\n";
          for ( int i1 = 0; i1 < sinks.isize( ); i1++ )
          for ( int j1 = 0; j1 < lfrom[i1].isize( ); j1++ )
          {    int i2 = lfrom[i1][j1].first, j2;
               for ( j2 = j1 + 1; j2 < lfrom[i1].isize( ); j2++ )
                    if ( lfrom[i1][j2].first != i2 ) break;
               if ( j2 - j1 > 1 ) // ignore singleton links
               {    int c1 = sink_to_comp[i1], c2 = source_to_comp[i2];
                    int len1 = complen[c1], len2 = complen[c2];
                    // if ( len1 > 1 && len2 > 1 )
                    {    cout << "from " << c1 << "[l=" << len1 << "]"
                              << " to " << c2 << "[l=" << len2 << "]" << " =";
                         for ( int k = j1; k < j2; k++ )
                         {    int rid = lfrom[i1][k].second;
                              int rpos1 = lfrom[i1][k].third; 
                              int rpos2 = lfrom[i1][k].fourth;
                              cout << " read_" << rid << "." << rpos1 
                                   << "-" << rpos2;    }
                         cout << endl;    }    }
               j1 = j2 - 1;    }

          // Data at this point:
          // - for each sink, lfrom[sink] = {(source,rid,rpos1,rpos2)}
          // - for each source, lto[source] = {(sink,rid,rpos1,rpos2)}.
          //
          // Create
          // Lfrom[sink] = {source,{rid,rpos1,rpos2}}
          // Lto[source] = {sink,{rid,rpos1,rpos2}}

          vec< vec< pair< int, vec< triple<int,int,int> > > > > Lfrom(nsinks);
          vec< vec< pair< int, vec< triple<int,int,int> > > > > Lto(nsources);
          for ( int i1 = 0; i1 < sinks.isize( ); i1++ )
          for ( int j1 = 0; j1 < lfrom[i1].isize( ); j1++ )
          {    int i2 = lfrom[i1][j1].first, j2;
               for ( j2 = j1 + 1; j2 < lfrom[i1].isize( ); j2++ )
                    if ( lfrom[i1][j2].first != i2 ) break;
               if ( j2 - j1 > 1 )
               {    vec< triple<int,int,int> > rx;
                    for ( int k = j1; k < j2; k++ )
                    {    int rid = lfrom[i1][k].second;
                         int rpos1 = lfrom[i1][k].third, rpos2 = lfrom[i1][k].fourth;
                         rx.push( rid, rpos1, rpos2 );    }
                    Lfrom[i1].push( i2, rx );
                    Lto[i2].push( i1, rx );    }
               j1 = j2 - 1;    }

          // Find things to link.  Criteria:
          // 1. Two things of reasonable size.
          // 2. Reasonably close.
          // 3. All other players can be positioned between or beyond.
          // 4. Nothing in between of good size.
          // 5. Resulting link is uncontested.

          cout << "\n";
          vec< pair<int,int> > accept;
          for ( int s1 = 0; s1 < nsinks; s1++ )
          for ( int j = 0; j < Lfrom[s1].isize( ); j++ )
          {    int s2 = Lfrom[s1][j].first;
               int c1 = sink_to_comp[s1], c2 = source_to_comp[s2];
               int len1 = complen[c1], len2 = complen[c2];
               if ( len1 < min_link_len || len2 < min_link_len ) continue;
               vec<int> d2;
               for ( int k = 0; k < Lfrom[s1][j].second.isize( ); k++ )
               {    d2.push_back( Lfrom[s1][j].second[k].third 
                         - Lfrom[s1][j].second[k].second );    }
               Sort(d2);
               int max_dist = Max(d2);
               double dist2 = Mean(d2);
               if ( max_dist > max_dist_link ) continue;
               cout << "\nexamining link from " << c1 << "[l=" << len1
                    << "] to " << c2 << "[l=" << len2 << "]" << endl; // XXX
               double max_over = 0;
               Bool intermediate = False;
               for ( int k = 0; k < Lfrom[s1].isize( ); k++ )
               {    if ( k == j ) continue;
                    int s3 = Lfrom[s1][k].first;
                    int c3 = source_to_comp[s3];
                    int len3 = complen[c3];
                    if ( len3 < min_link_len ) continue;
                    vec<int> d3;
                    for ( int j = 0; j < Lfrom[s1][k].second.isize( ); j++ )
                    {    d3.push_back( Lfrom[s1][k].second[j].third
                              - Lfrom[s1][k].second[j].second );    }
                    Sort(d3);
                    double dist3 = Mean(d3);
                    cout << "comparing " << c2 << " at " << dist2 << "-"
                         << dist2 + len2 << " to " << c3 << "[l=" << len3
                         << "] at " << dist3
                         << "-" << dist3 + len3 << endl;    
                    double left2 = dist2, right2 = dist2 + len2;
                    double left3 = dist3, right3 = dist3 + len3;
                    max_over = Max( max_over, 
                         IntervalOverlap( left2, right2, left3, right3 ) );
                    if ( len3 >= min_link_len && left3 >= -1 && right3 <= left2 + 1 )
                    {    cout << "looks like an intermediate" << endl;
                         intermediate = True;    }    }
               for ( int k = 0; k < Lto[s2].isize( ); k++ )
               {    int s3 = Lto[s2][k].first;
                    if ( s3 == s1 ) continue;
                    int c3 = sink_to_comp[s3];
                    int len3 = complen[c3];
                    if ( len3 < min_link_len ) continue;
                    vec<int> d3;
                    for ( int j = 0; j < Lto[s2][k].second.isize( ); j++ )
                    {    d3.push_back( Lto[s2][k].second[j].third
                              - Lto[s2][k].second[j].second );    }
                    Sort(d3);
                    double dist3 = Mean(d3);
                    cout << "comparing " << c2 << " at " << dist2 << "-"
                         << dist2 + len2 << " to " << c3 << "[l=" << len3
                         << "] at " 
                         << dist2 - dist3 - len3 << "-" << dist2 - dist3 << endl;    
                    max_over = Max( max_over, IntervalOverlap( dist2, dist2 + len2,
                         dist2 - dist3 - len3, dist2 - dist3 ) );    }
               cout << "max interval overlap = " << max_over << endl;
               if (intermediate) cout << "have intermediate" << endl;
               if ( max_over <= max_max_over && !intermediate ) 
               {    accept.push( s1, s2 );
                    cout << "accepting" << endl;    }    }

          // Remove contested links.

          vec<Bool> adel( accept.size( ), False );
          Sort(accept);
          for ( int pass = 1; pass <= 2; pass++ )
          {    for ( int i = 0; i < accept.isize( ); i++ )
                    swap( accept[i].first, accept[i].second );
               SortSync( accept, adel );
               for ( int i = 0; i < accept.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < accept.isize( ); j++ )
                         if ( accept[j].first != accept[i].first ) break;
                    if ( j - i > 1 )
                    {    for ( int k = i; k < j; k++ )
                              adel[k] = True;    }
                    i = j - 1;    }    }
          EraseIf( accept, adel );
          cout << "\n";
          DPRINT( accept.size( ) );

          // On last pass, make any strong join.

          if ( lpass == link_passes )
          {    for ( int s1 = 0; s1 < nsinks; s1++ )
               for ( int j = 0; j < Lfrom[s1].isize( ); j++ )
               {    int s2 = Lfrom[s1][j].first;
                    int c1 = sink_to_comp[s1], c2 = source_to_comp[s2];
                    int len1 = complen[c1], len2 = complen[c2];
                    if ( len1 < 10 || len2 < 10 ) continue; // length at last 10
                    vec<int> d2;
                    for ( int k = 0; k < Lfrom[s1][j].second.isize( ); k++ )
                    {    d2.push_back( Lfrom[s1][j].second[k].third 
                              - Lfrom[s1][j].second[k].second );    }
                    Sort(d2);
                    if ( d2.size( ) < 5 ) continue; // require 5 links
                    int max_dist = Max(d2);
                    if ( max_dist > 5 ) continue; // distance at most 5
                    accept.push( s1, s2 );    }    }
          UniqueSort(accept);

          // Join accepted links.

          for ( int i = 0; i < accept.isize( ); i++ )
          {    int s1 = accept[i].first, s2 = accept[i].second;
               int v1 = sinks[s1], v2 = sources[s2];
               for ( int j = 0; j < Lfrom[s1][0].second.isize( ); j++ )
               {    int rid = Lfrom[s1][0].second[j].first;
                    int rpos1 = Lfrom[s1][0].second[j].second;
                    int rpos2 = Lfrom[s1][0].second[j].third;
                    if ( rpos1 == rpos2 )
                    {    int len = 0;
                         vec< pair<int,int> > id_pos;
                         G.AddEdge( v1, v2, bng_edge( len, id_pos ) );    }
                    else
                    {    int v = v1;
                         for ( int m = rpos1; m < rpos2; m++ )
                         {    int len = X[rid][m];
                              vec< pair<int,int> > id_pos;
                              id_pos.push( rid, m );
                              int n;
                              if ( m < rpos2 - 1 )
                              {    n = G.N( );
                                   G.AddVertices(1);    }
                              else n = v2;
                              G.AddEdge( v, n, bng_edge( len, id_pos ) );
                              v = n;    }    }    }    }    }    }

     // Eliminate some zero-length edges.  Should understand how these arise.

     dels.clear( );
     for ( int v = 0; v < G.N( ); v++ )
     {    vec<int> x = G.From(v);
          UniqueSort(x);
          for ( int i = 0; i < x.isize( ); i++ )
          {    int w = x[i];
               Bool zero = True;
               for ( int j = 0; j < G.From(v).isize( ); j++ )
               {    if ( G.From(v)[j] == w && G.O( G.IFrom(v,j) ).Len( ) > 0 ) 
                         zero = False;    }
               if ( !zero || w == v || G.HasEdge( w, v ) ) continue;
               for ( int j = 0; j < G.From(v).isize( ); j++ )
               {    if ( BinMember( x, G.From(v)[j] ) )
                         dels.push_back( G.IFrom(v,j) );    }
               G.TransferEdges( w, v );    }    }

     // Flatten certain bubbles.

     cout << Date( ) << ": flattening" << endl;
     for ( int vv = 0; vv < G.N( ); vv++ )
     {    int v = vv, w;
          if ( !IsBubbleStart( G, v, w ) ) continue;

          // Find length in edges of each bubble branch.

          vec<int> ne( G.From(v).size( ), 0 );
          for ( int i = 0; i < G.From(v).isize( ); i++ )
          {    int t = G.From(v)[i], e = G.IFrom(v,i);
               if ( G.O(e).Len( ) == 0 && t == w );
               else
               {    while(1)
                    {    ne[i]++;
                         if ( t == w ) break;
                         t = G.From(t)[0];    }    }    }
          int low = Min(ne), high = Max(ne);

          // Deal with zero edges.  Need to investigate how these are created.

          /*
          if ( low == 0 && high == 0 )
          {    for ( int j = 0; j < G.From(v).isize( ); j++ )
                    dels.push_back( G.IFrom( v, j ) );
               G.TransferEdges( v, w );
               continue;    }
          */

          // Try to combine on the left.

          while( low >= 1 && ( low == high || low >= 2 ) )
          {    vec<int> lens( G.From(v).size( ) );
               for ( int j = 0; j < G.From(v).isize( ); j++ )
                    lens[j] = G.O( G.IFrom(v,j) ).Len( );
               Sort(lens);
               if ( lens.front( ) == 0 || lens.back( ) - lens.front( ) > 1000 ) 
                    break;
               int len = Mean(lens);
               vec< pair<int,int> > id_pos;
               for ( int j = 0; j < G.From(v).isize( ); j++ )
               {    int e = G.IFrom(v,j);
                    id_pos.append( G.O(e).IdPos( ) );
                    dels.push_back(e);    }
               if ( high == 1 )
               {    G.AddEdge( v, w, bng_edge( len, id_pos ) );
                    break;    }
               int vnew = G.N( );
               G.AddVertices(1);
               vec<int> f = G.From(v);
               for ( int j = 0; j < f.isize( ); j++ )
                    G.TransferEdges( f[j], vnew );
               G.AddEdge( v, vnew, bng_edge( len, id_pos ) );
               v = vnew;
               low--;
               high--;    }    }
     G.DeleteEdges(dels);
     G.RemoveDeadEdgeObjects( );
     G.RemoveEdgelessVertices( );

     // Remove small components.

     cout << Date( ) << ": removing small components" << endl;
     {    dels.clear( );
          vec<vec<int>> comp;
          G.ComponentsE(comp);
          for ( int i = 0; i < comp.isize( ); i++ )
               if ( comp[i].isize( ) < min_edges ) dels.append( comp[i] );
          G.DeleteEdges(dels);
          G.RemoveDeadEdgeObjects( );
          G.RemoveEdgelessVertices( );    }
     BinaryWriter::writeFile( OUT_DIR + "/final.graph", G );

     // Print graph.

     Ofstream( gout, OUT_DIR + "/graph.dot" );
     vec<String> edge_labels;
     for ( int e = 0; e < G.E( ); e++ )
     {    ostringstream out;
          if ( EDGE_LABEL_MODE == 1 ) out << G.EdgeObject(e).Len( );
          if ( EDGE_LABEL_MODE == 2 )
          {    out << G.EdgeObject(e).Len( ) << ":{" 
                    << printSeq( G.EdgeObject(e).Rids( ) ) << "}";    }
          if ( EDGE_LABEL_MODE == 3 )
          {    vec<int> len;
               vec< pair<int,int> > x = G.EdgeObject(e).IdPos( );
               for ( int i = 0; i < x.isize( ); i++ )
               {    int rid = x[i].first, rpos = x[i].second;
                    len.push_back( X[rid][rpos] );    }
               Sort(len);
               out << G.EdgeObject(e).Len( ) << "(" << printSeq(len) << ")";    }
          edge_labels.push_back( out.str( ) );    }
     digraphE<String> H( G, edge_labels );
     H.DOT0(gout);    

     // Generate another more compressed version of the graph.  

     digraphE< vec<int> > I;
     Compress( G, I );

     // Build disco and bng assemblies.

     vec<vec<int>> disco_locs( I.E( ) );
     if (DISCO)
     {    cout << Date( ) << ": loading disco assembly" << endl;
          String dir = "/wga/scr4/jaffe/GapToy/50921.HCC1143+BL/a.final";
          HyperBasevectorX hb;
          BinaryReader::readFile( dir + "/a.hbx", &hb );
          vec<vec<vec<vec<int>>>> lines;
          BinaryReader::readFile( dir + "/a.lines", &lines );
          const int min_line = 10000;
          const int max_ignored_indel = 120;
          String cut = "GCTCTTC";
          vec<vec<vec<int>>> texps;
          Bool verbose = True;
          Ofstream( out, OUT_DIR + "/disco.lines" );
          Texps( hb, lines, cut, min_line, max_ignored_indel, texps, verbose, out );
          vec<vec<double>> X;
          vec<String> ids;
          cout << Date( ) << ": building X" << endl;
          for ( int i = 0; i < texps.isize( ); i++ )
          for ( int j = 0; j < texps[i].isize( ); j++ )
          {    if ( texps[i][j].empty( ) ) continue;
               vec<double> x( texps[i][j].size( ) );
               for ( int l = 0; l < texps[i][j].isize( ); l++ )
                    x[l] = texps[i][j][l];
               X.push_back(x);    
               ids.push_back( "L" + ToString(i) );    }
          cout << "X size from disco = " << X.size( ) << endl;
          int div = X.size( );
          for ( int e = 0; e < I.E( ); e++ )
          {    if ( I.O(e).size( ) < 4 ) continue;
               vec<double> x( I.O(e).size( ) );
               for ( int i = 0; i < x.isize( ); i++ )
                    x[i] = I.O(e)[i];
               X.push_back(x);
               ids.push_back( "B" + ToString(e) );    }
          cout << "X size from disco and bng = " << X.size( ) << endl;
          vec<int> S( X.size( ), vec<int>::IDENTITY );
          cout << Date( ) << ": finding matches from disco to bng assembly" << endl;
          vec< quad<int,int,int,int> > zmatches, zmatches2;
          {    const int w = 4;
               const int cutoff = 1000;
               Zmatches( X, S, zmatches, False, div, w, cutoff );    }
          {    const int w = 3;
               const int cutoff = 500;
               Zmatches( X, S, zmatches2, False, div, w, cutoff );    }
          zmatches.append(zmatches2);
          ParallelUniqueSort(zmatches);
          cout << "DONE" << endl;

          // Make the actual disco/bng alignments.

          cout << Date( ) << ": generating disco/bng alignments" << endl;
          int align_calls = 0;
          vec< triple< int, int, vec<int> > > aligns;
          vec<double> scores;
          double pclock = WallClockTime( );
          double ac1 = 0, ac2 = 0, ac3 = 0, ac4 = 0, ac5 = 0;
          vec<int64_t> zstarts( X.size( ) + 1, -1 );
          for ( int64_t i = zmatches.jsize( ) - 1; i >= 0; i-- )
               zstarts[ zmatches[i].first ] = i;
          zstarts[ X.size( ) ] = zmatches.size( );
          for ( int j = X.isize( ) - 1; j >= 0; j-- )
               if ( zstarts[j] < 0 ) zstarts[j] = zstarts[j+1];
          // ALIGN_LOGGING = True;

          // Different heuristics.

          const int min_direct = 8;
     
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int s1 = 0; s1 < S.isize( ); s1++ )
          {    int i1 = S[s1];
               int64_t low = zstarts[i1], high = zstarts[i1+1];
               if ( low == high ) continue;
               BngAlignCore( X, zmatches, low, high, S, ALIGN_LOGGING, 
                    SHOW_ALL_SEEDS, ALIGN_DETAILS, ONE_SEED, RFILTER, PRINT_FAIL,
                    align_calls, aligns, scores, ac1, ac2, ac3, ac4, ac5,
                    min_direct, max_score, max_total_err, s1, i1 );     }
          UniqueSortSync( aligns, scores );
          cout << "\n";
          int naligns = aligns.size( );
          PRINT2( naligns, align_calls );
          PRINT5( ac1, ac2, ac3, ac4, ac5 );
          cout << "\n" << TimeSince(pclock) << " used in main loop" << endl;

          // Unwind the alignments.

          for ( int i = 0; i < aligns.isize( ); i++ )
          {    int i1 = aligns[i].first, i2 = aligns[i].second;
               cout << ids[i1] <<  " aligned to " << ids[i2] << endl;    
               disco_locs[ ids[i1].After( "B" ).Int( ) ].push_back(
                    ids[i2].After( "L" ).Int( ) );    }
          for ( int e = 0; e < I.E( ); e++ )
               UniqueSort( disco_locs[e] );    }

     // Compute top sizes.

     digraphE< vec<int> > B(I);
     cout << Date( ) << ": splaying" << endl;
     SplayAll(B);
     B.RemoveEdgelessVertices( );
     vec<int64_t> clens;
     vec<vec<int>> compy;
     cout << Date( ) << ": getting components" << endl;
     B.Components(compy);
     PRINT( compy.size( ) );
     cout << Date( ) << ": computing top sizes" << endl;
     vec<vec<int>> clens_all( compy.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < compy.isize( ); i++ )
     {    digraphE< vec<int> > J = B.Subgraph( compy[i] );
          vec<int> sources, sinks;
          J.Sources(sources), J.Sinks(sinks);
          Bool OK = False;
          if ( SPLAY && sources.solo( ) && sinks.solo( ) && J.Acyclic( ) )
          {    vec<vec<int>> paths;
               Bool ok = J.EdgePaths( sources[0], sinks[0], paths, -1, -1, 10000 );
               if (ok)
               {    vec<int> lens;
                    for ( int j = 0; j < paths.isize( ); j++ )
                    {    int L = 0;
                         for ( int m = 0; m < paths[j].isize( ); m++ )
                              L += Sum( J.O( paths[j][m] ) );
                         lens.push_back(L);    }
                    Sort(lens);
                    if ( lens.nonempty( ) ) // not sure if needed, or why it would be
                    {    clens_all[i].push_back( Median(lens) );
                         OK = True;    }    }    }
          if ( !OK )
          {    for ( int e = 0; e < J.E( ); e++ )
                    clens_all[i].push_back( Sum( J.O(e) ) );    }    }
     for ( int i = 0; i < compy.isize( ); i++ )
          clens.append( clens_all[i] );
     cout << Date( ) << ": sorting lens" << endl;
     ReverseSort(clens);
     cout << "top compressed sizes =";
     for ( int j = 0; j < Min( 4, clens.isize( ) ); j++ )
          cout << " " << clens[j];

     // Compute total size.

     int64_t totalx = 0;
     for ( int e = 0; e < I.E( ); e++ )
          totalx += Sum( I.O(e) );
     cout << "\ntotal compressed length = " << ToStringAddCommas(totalx) << endl;

     // Create graph I2 parallel to I that carries G edge ids.

     digraphE< vec<int> > I2;
     {    vec< vec<int> > vedges;
          for ( int e = 0; e < G.E( ); e++ )
          {    vec<int> v = {e};
               vedges.push_back(v);    }
          I2.Initialize( G, vedges );
          for ( int v = 0; v < I2.N( ); v++ )
          {    if ( !I2.To(v).solo( ) || !I2.From(v).solo( ) ) continue;
               if ( I2.To(v)[0] == v || I2.From(v)[0] == v
                    || I2.To(v)[0] == I2.From(v)[0] )
               {    continue;    }
               int e1 = I2.ITo( v, 0 ), e2 = I2.IFrom( v, 0 );
               vec<int> x = I2.EdgeObject(e1);
               x.append( I2.EdgeObject(e2) );
               I2.JoinEdges( v, x );    }
          I2.RemoveDeadEdgeObjects( );
          I2.RemoveEdgelessVertices( );    }

     // Generate compressed dot files.

     cout << "\ncompressed graph has " << I.E( ) << " edges" << endl;
     vec<String> vedge_labels;
     for ( int e = 0; e < I.E( ); e++ )
     {    ostringstream out;
          out << e << ":";
          const vec<int>& x = I.O(e);
          const int show = 4;
          const int mshow = 10;
          if ( x.isize( ) <= mshow ) out << printSeq(x);
          else
          {    for ( int j = 0; j < show; j++ )
               {    if ( j > 0 ) out << ",";
                    out << x[j];    }
               out << ",[" << x.isize( ) - 2*show << "]";
               for ( int j = x.isize( ) - show; j < x.isize( ); j++ )
                    out << "," << x[j];    }
          if ( x.size( ) > 1 ) out << " = " << Sum(x);

          // Show number of reads that support edge.

          vec<int> rids;
          const vec<int>& y = I2.O(e);
          for ( int i = 0; i < y.isize( ); i++ )
               rids.append( G.O( y[i] ).Rids( ) );
          UniqueSort(rids);
          out << " (" << rids.size( ) << ")";

          // Add disco locs.

          for ( int j = 0; j < disco_locs[e].isize( ); j++ )
               out << " L" << disco_locs[e][j];

          // Save label.

          vedge_labels.push_back( out.str( ) );    }
     {    vec<vec<int>> comp;
          I.Components(comp);
          const int comp_per_dot = 200;
          SystemSucceed( "/bin/rm -f " + OUT_DIR + "/comp.*.dot" );
          for ( int i = 0; i < comp.isize( ); i += comp_per_dot )
          {    int start = i, stop = Min( comp.isize( ), i + comp_per_dot );
               vec<int> v;
               for ( int j = start; j < stop; j++ )
                    v.append( comp[j] );
               Sort(v);
               int id = ( i / comp_per_dot ) + 1;
               digraphE< vec<int> > T = I.Subgraph(v);
               Ofstream( hout, OUT_DIR + "/comp." + ToString(id) + ".dot" );
               vec<String> vel;
               for ( int j = 0; j < v.isize( ); j++ )
               for ( int k = 0; k < I.From( v[j] ).isize( ); k++ )
                    vel.push_back( vedge_labels[ I.IFrom( v[j], k ) ] );
               digraphE<String> IS( T, vel );
               IS.DOT0(hout);    }

          // Create dot files of special components.
          
          vec<int> specials;
          for ( int i = 0; i < comp.isize( ); i++ )
          {    digraphE< vec<int> > J = I.Subgraph( comp[i] );
               vec<int> sources, sinks;
               J.Sources(sources), J.Sinks(sinks);
               if ( !sources.solo( ) || !sinks.solo( ) ) specials.push_back(i);    }
          const int specials_per_dot = 20;
          for ( int i = 0; i < specials.isize( ); i += specials_per_dot )
          {    int start = i, stop = Min( specials.isize( ), i + specials_per_dot );
               vec<int> v;
               for ( int j = start; j < stop; j++ )
                    v.append( comp[ specials[j] ] );
               Sort(v);
               int id = ( i / specials_per_dot ) + 1;
               digraphE< vec<int> > T = I.Subgraph(v);
               Ofstream( sout, OUT_DIR + "/comp.special." + ToString(id) + ".dot" );
               vec<String> vel;
               for ( int j = 0; j < v.isize( ); j++ )
               for ( int k = 0; k < I.From( v[j] ).isize( ); k++ )
                    vel.push_back( vedge_labels[ I.IFrom( v[j], k ) ] );
               digraphE<String> IS( T, vel );
               IS.DOT0(sout);    }    }

     // Create "fasta" file.

     {    Ofstream( hout2, OUT_DIR + "/comp.fasta" );
          for ( int e = 0; e < I.E( ); e++ )
          {    hout2 << ">" << e << endl;
               hout2 << printSeq( I.EdgeObject(e) ) << endl;    }    }

     // Look for outside reads.

     /*
     vec<int> bones;
     for ( int e = 0; e < G.E( ); e++ )
     {    const vec<int>& x = G.EdgeObject(e).Rids( );
          for ( int j = 0; j < x.isize( ); j++ )
               if ( x[j] >= nx ) bones.push_back( x[j] );    }
     PRINT2( X.size( ), X2.size( ) );
     cout << "outside reads = " << printSeq(bones) << endl;
     cout << "\n";
     */

     // Look for undetected duplications.

     vec<int> all;
     for ( int e = 0; e < G.E( ); e++ )
          all.push_back( G.O(e).Len( ) );
     UniqueSort(all);
     int dups = G.E( ) - all.isize( );
     cout << PERCENT_RATIO( 3, dups, G.E( ) ) << " length duplication" << endl;

     // Print summary stats.

     int nedges = G.EdgeObjectCount( ), ncomponents = G.NComponents( );
     int nloops = 0, nbranches = 0;
     for ( int v = 0; v < G.N( ); v++ )
     {    for ( int j = 0; j < G.From(v).isize( ); j++ )
               if ( G.From(v)[j] == v ) nloops++;
          if ( G.From(v).size( ) > 1 ) nbranches += G.From(v).size( ) - 1;
          if ( G.To(v).size( ) > 1 ) nbranches += G.To(v).size( ) - 1;    }
     PRINT4( nedges, ncomponents, nloops, nbranches );
     cout << TimeSince(all_clock) << " used" << endl;
     cout << "peak mem usage = " << PeakMemUsageGB( ) << " GB\n";
     cout << "\n" << Date( ) << ": done" << endl;
     Scram(0);    }
