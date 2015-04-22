///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DisplayLinear: display a linearized version of a SupportedHyperBasevector
// and compute some stats for it.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

// MakeDepend: dependency QueryLookupTable

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperEfasta.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"

void GetCells2( const HyperBasevector& he, vec<vec<int>>& cells,
     const int max_internal_edge )
{    cells.clear( );

     // Go through the vertices.

     #pragma omp parallel for
     for ( int x = 0; x < he.N( ); x++ )
     for ( int i = 0; i < he.From(x).isize( ); i++ )
     {    int v = he.From(x)[i];
          int e = he.EdgeObjectIndexByIndexFrom( x, i );

          //   x ---e---> v bounds on the left

          vec< vec<int> > ecells;
          vec<int> fs;

          // Find a successor w of v that is not v.

          vec<int> suc, pre;
          he.GetSuccessors1( v, suc );
          suc.EraseValue(v);
          for ( int wi = 0; wi < suc.isize( ); wi++ )
          {    int w = suc[wi];
               if ( w == v ) continue;

               // Let w ---f---> y bound on the right.

               for ( int j = 0; j < he.From(w).isize( ); j++ )
               {    int y = he.From(w)[j];
                    int f = he.EdgeObjectIndexByIndexFrom( w, j );

                    Bool bad = False;

                    // Now we posit
                    //   x ---e---> v    // CELL //     w ---f---> y
                    
                    // Start from v and keep walking.  The only thing you
                    // can't do is walk through f.

                    vec<int> sucv;
                    sucv.push_back(v);
                    while(1)
                    {    if (bad) break;
                         Bool progress = False;
                         for ( int l = 0; l < sucv.isize( ); l++ )
                         {    if (bad) break;
                              int s = sucv[l];
                              for ( int m = 0; m < he.From(s).isize( ); m++ )
                              {    int t = he.From(s)[m];
                                   if ( Member( sucv, t ) ) continue;
                                   int g = he.EdgeObjectIndexByIndexFrom( s, m );
                                   if ( g == f ) continue;
                                   if ( he.EdgeLengthKmers(g) > max_internal_edge )
                                   {    bad = True;
                                        break;    }
                                   sucv.push_back(t);
                                   progress = True;    }    }
                         if ( !progress ) break;    }

                    if (bad) continue;

                    Sort(sucv);

                    // Go back the other way.

                    vec<int> prew;
                    prew.push_back(w);
                    while(1)
                    {    Bool progress = False;
                         for ( int l = 0; l < prew.isize( ); l++ )
                         {    int s = prew[l];
                              for ( int m = 0; m < he.To(s).isize( ); m++ )
                              {    int t = he.To(s)[m];
                                   if ( Member( prew, t ) ) continue;
                                   if ( he.EdgeObjectIndexByIndexTo( s, m ) == e )
                                        continue;
                                   prew.push_back(t);
                                   progress = True;    }    }
                         if ( !progress ) break;    }
                    Sort(prew);

                    // Both views should be the same.

                    if ( sucv != prew ) continue;

                    // Now we know that v and w bound a cell.

                    vec<int> between = sucv;
                    between.EraseValue(v);
                    between.EraseValue(w);

                    vec<int> cell;
                    cell.push_back(v);
                    cell.append(between);
                    cell.push_back(w);

                    ecells.push_back(cell);
                    fs.push_back(f);    }    }

          int m = 1000000000;
          for ( int z = 0; z < ecells.isize( ); z++ )
               m = Min( m, ecells[z].isize( ) );
          for ( int z = 0; z < ecells.isize( ); z++ )
          {    if ( ecells[z].isize( ) > m ) continue;
               const vec<int>& cell = ecells[z];
               int f = fs[z];

               #pragma omp critical
               {    cells.push_back(cell);    

                    /*
                    if ( cells.size( ) % 100 == 0 ) PRINT( cells.size( ) );
                    vec<int> edges;
                    for ( int j = 0; j < cell.isize( ); j++ )
                    {    int v = cell[j];
                         for ( int l = 0; l < he.From(v).isize( ); l++ )
                         {    int w = he.From(v)[l];
                              if ( !Member( cell, w ) ) continue;
                              int e = he.EdgeObjectIndexByIndexFrom( v, l );
                              edges.push_back(e);    }    }
                    Sort(edges);
                    cout << "\ncell: " << printSeq(cell) << endl;
                    cout << "edges: " << printSeq(edges) << endl;
                    PRINT2( e, f );
                    */

                              }    }    }
     Sort(cells);    }

vec<int> CloseVertices( const SupportedHyperBasevector& shb,
     const int v, const int max_dist )
{    vec< triple<int,int,Bool> > X;
     if ( max_dist < 0 ) return vec<int>( );
     X.push( v, 0, False );
     while(1)
     {    Bool changed = False;
          for ( int j = 0; j < X.isize( ); j++ )
          {    if ( X[j].third ) continue;
               X[j].third = True;
               int w = X[j].first;
               for ( int l = 0; l < shb.From(w).isize( ); l++ )
               {    int e = shb.EdgeObjectIndexByIndexFrom( w, l );
                    int d = X[j].second + shb.EdgeLengthKmers(e);
                    if ( d > max_dist ) continue;
                    int z = shb.From(w)[l];
                    Bool found = False;
                    for ( int k = 0; k < X.isize( ); k++ )
                    {    if ( X[k].first == z )
                         {    found = True;
                              if ( X[k].second > d )
                              {    X[k].second = d;
                                   X[k].third = False;
                                   changed = True;    }
                              break;    }    }
                    if ( !found )
                    {    X.push( z, d, False );
                         changed = True;    }    }    }
          if ( !changed ) break;    }
     vec<int> answer;
     for ( int j = 0; j < X.isize( ); j++ )
          answer.push_back( X[j].first );
     Sort(answer);
     return answer;    }

int main(int argc, char *argv[])
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD, "looks for IN_HEAD.shbv");
     CommandArgument_String_Doc(OUT_HEAD, "generate OUT_HEAD.dot");
     EndCommandArguments;

     // Load assembly.

     SupportedHyperBasevector shb;
     BinaryReader::readFile( IN_HEAD + ".shbv", &shb );

     // Define heuristics.

     const int min_split = 2000;
     const int max_to_cut = 700;
     const int min_solo = 500;

     // Splay 'branch' vertices.

     int nv = shb.N( );
     vec<Bool> splay( nv, False );
     for ( int pass = 1; pass <= 2; pass++ )
     {    shb.Reverse( );
          for ( int v = 0; v < nv; v++ )
          {    for ( int i1 = 0; i1 < shb.From(v).isize( ); i1++ )
               {    int w1 = shb.From(v)[i1];
                    int e1 = shb.EdgeObjectIndexByIndexFrom( v, i1 );
                    vec<int> c1 = CloseVertices( shb, w1, 
                         min_split - shb.EdgeLengthKmers(e1) );
                    for ( int i2 = i1 + 1; i2 < shb.From(v).isize( ); i2++ )
                    {    int w2 = shb.From(v)[i2];
                         int e2 = shb.EdgeObjectIndexByIndexFrom( v, i2 );
                         vec<int> c2 = CloseVertices( shb, w2, 
                              min_split - shb.EdgeLengthKmers(e2) );
                         if ( !Meet(c1, c2) && !BinMember(c1, v) 
                              && !BinMember(c2, v) )
                         {    splay[v] = True;
                              break;    }    }
                    if ( splay[v] ) break;    }    }    }
     for ( int v = 0; v < nv; v++ )
          if ( splay[v] ) shb.SplayVertex(v);

     // Cut off bits of fuzz at the ends of components.

     while(1)
     {    vec<int> L( shb.EdgeObjectCount( ) );
          for ( int i = 0; i < shb.EdgeObjectCount( ); i++ )
               L[i] = shb.EdgeLengthKmers(i);
          digraphE<int> G( shb, L );
          vec<Bool> cut( shb.N( ), False );
          for ( int pass = 1; pass <= 2; pass++ )
          {    G.Reverse( );
               for ( int v = 0; v < shb.N( ); v++ )
               {    vec< pair<int,int> > s;
                    G.GetSuccessors( {v}, s );
                    int M = 0;
                    for ( int j = 0; j < s.isize( ); j++ )
                         M = Max( M, s[j].second );
                         if ( M <= max_to_cut ) cut[v] = True;    }    }
          int cuts = 0;
          for ( int v = 0; v < cut.isize( ); v++ )
          {    if ( shb.From(v).solo( ) && shb.To(v).solo( )
                    && shb.From(v)[0] == v )
               {    continue;    }
               if ( cut[v] && shb.From(v).size( ) + shb.To(v).size( ) > 1 )
               {    cuts++;
                    shb.SplayVertex(v);    }    }
          if ( cuts == 0 ) break;    }

     // Delete tiny solo edges.

     vec<int> dels;
     for ( int v = 0; v < shb.N( ); v++ )
     {    if ( !shb.From(v).solo( ) || !shb.To(v).empty( ) ) continue;
          int w = shb.From(v)[0];
          if ( w == v ) continue;
          int e = shb.EdgeObjectIndexByIndexFrom( v, 0 );
          if ( !shb.To(w).solo( ) || !shb.From(w).empty( ) ) continue;
          if ( shb.EdgeLengthKmers(e) > min_solo ) continue;
          dels.push_back(e);    }

     // Delete reverse complement components.

     vec< vec<int> > comp;
     shb.ComponentsE(comp);
     vec<Bool> used( shb.EdgeObjectCount( ), False );
     for ( int c = 0; c < comp.isize( ); c++ )
     {    if ( comp[c].empty( ) ) continue;
          int e = comp[c][0];
          if ( used[ shb.Inv(e) ] )
          {    for ( int j = 0; j < comp[c].isize( ); j++ )
                    dels.push_back( comp[c][j] );    }
          else
          {    for ( int j = 0; j < comp[c].isize( ); j++ )
                    used[ comp[c][j] ] = True;    }    }
     shb.DeleteEdges(dels);    

     // Compute stats for the linearized assembly.

     vec< vec<int> > comp2;
     shb.ComponentsE(comp2);
     vec<int> edge_sizes, comp_sizes;
     for ( int c = 0; c < comp2.isize( ); c++ )
     {    if ( comp2[c].empty( ) ) continue;
          int nkmers = 0;
          for ( int j = 0; j < comp2[c].isize( ); j++ )
          {    int n = shb.EdgeLengthKmers( comp2[c][j] );
               edge_sizes.push_back(n);
               nkmers += n;    }
          comp_sizes.push_back(nkmers);    }
     Sort(edge_sizes), Sort(comp_sizes);
     cout << "N50 edge size = " << N50(edge_sizes) << endl;
     int e1 = 0, e2 = 0, e3 = 0, e4 = 0;
     for ( int i = 0; i < edge_sizes.isize( ); i++ )
     {    int n = edge_sizes[i];
          if ( n < 1000 ) e1++;
          else if ( n < 10000 ) e2++;
          else if ( n < 100000 ) e3++;
          else e4++;    }
     cout << e1 << " edges < 1 kb\n";
     cout << e2 << " edges between 1 and 10 kb\n";
     cout << e3 << " edges between 10 and 100 kb\n";
     cout << e4 << " edges > 100 kb\n";
     cout << "\nN50 chunk size = " << N50(comp_sizes) << endl;
     int c1 = 0, c2 = 0, c3 = 0, c4 = 0;
     for ( int i = 0; i < comp_sizes.isize( ); i++ )
     {    int n = comp_sizes[i];
          if ( n < 1000 ) c1++;
          else if ( n < 10000 ) c2++;
          else if ( n < 100000 ) c3++;
          else c4++;    }
     cout << c1 << " chunks < 1 kb\n";
     cout << c2 << " chunks between 1 and 10 kb\n";
     cout << c3 << " chunks between 10 and 100 kb\n";
     cout << c4 << " chunks > 100 kb\n";

     // Display the linearized assembly.  Likely shortenable.

     vec<Bool> xused;
     shb.Used(xused);
     vec<Bool> invisible( shb.EdgeObjectCount( ), False );
     vec<String> edge_color( shb.EdgeObjectCount( ) );
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          if ( !xused[e] ) invisible[e] = True;
     long_logging logc( "" );
     vec<String> edge_names( shb.EdgeObjectCount( ) );
     for ( int i = 0; i < shb.EdgeObjectCount( ); i++ )
          edge_names[i] = ToString(i);
     shb.DumpDot( OUT_HEAD, invisible, edge_color, logc, False, edge_names );

     // Find cells.

     HyperBasevector he(shb);
     vec<vec<int>> cells, ecells;
     const int max_internal_edge = 5000;
     GetCells2( he, cells, max_internal_edge );
     ecells.resize( cells.size( ) );
     for ( int i = 0; i < cells.isize( ); i++ )
     {    const vec<int>& cell = cells[i];
          for ( int j = 0; j < cell.isize( ); j++ )
          {    int v = cell[j];
               for ( int l = 0; l < he.From(v).isize( ); l++ )
               {    int w = he.From(v)[l];
                    if ( !Member( cell, w ) ) continue;
                    int e = he.EdgeObjectIndexByIndexFrom( v, l );
                    ecells[i].push_back(e);    }    }
          Sort( ecells[i] );    }
     vec<Bool> toobig( cells.size( ), False );
     for ( int i = 0; i < cells.isize( ); i++ )
     {    for ( int j = i - 1; j >= 0; j-- )
          {    if ( cells[j][0] != cells[i][0] ) break;
               if ( cells[j].isize( ) < cells[i].isize( ) )
               {    toobig[i] = True;
                    break;    }    }
          for ( int j = i + 1; j < cells.isize( ); j++ )
          {    if ( cells[j][0] != cells[i][0] ) break;
               if ( cells[j].isize( ) < cells[i].isize( ) )
               {    toobig[i] = True;
                    break;    }    }    }
     vec<Bool> crush( cells.size( ), False );
     vec<int> in( cells.size( ) ), out( cells.size( ) );
     for ( int i = 0; i < cells.isize( ); i++ )
     {    if ( toobig[i] ) continue;
          const vec<int> &cell = cells[i], &edges = ecells[i];

          // Define in edge and out edge.

          vec<int> ins, outs;
          int v = cell.front( ), w = cell.back( );
          for ( int j = 0; j < he.To(v).isize( ); j++ )
          {    if ( Member( cell, he.To(v)[j] ) ) continue;
               ins.push_back( he.EdgeObjectIndexByIndexTo( v, j ) );    }
          for ( int j = 0; j < he.From(w).isize( ); j++ )
          {    if ( Member( cell, he.From(w)[j] ) ) continue;
               outs.push_back( he.EdgeObjectIndexByIndexFrom( w, j ) );    }
          if ( ins.size( ) != 1 || outs.size( ) != 1 ) continue;

          const int max_crush = 5;
          if ( edges.isize( ) <= max_crush ) crush[i] = True;
          in[i] = ins[0], out[i] = outs[0];

          cout << "\nexamining cell " << i << "\n";
          cout << "vertices: " << printSeq(cell) << endl;
          cout << "incoming edge: " << printSeq(ins) << endl;
          cout << "outgoing edge: " << printSeq(outs) << endl;
          cout << "edges out: " << printSeq( he.FromEdgeObj( cell[0] ) ) << endl;
          cout << "edges: " << printSeq(edges) << endl;    

          vec< vec<int> > paths;
          const int max_iterations = 1000;
          const int max_paths = 15;
          Bool OK = he.EdgePaths( cell.front( ), cell.back( ), paths,
                    -1, max_paths, max_iterations );
          if ( !OK ) continue;
          if ( paths.isize( ) > max_paths ) continue;
          if ( paths.size( ) <= 1 ) continue;
     
          vecbasevector bpaths;
          for ( int i = 0; i < paths.isize( ); i++ )
               bpaths.push_back( shb.Cat( paths[i] ) );

          if ( paths.size( ) == 2 )
          {    alignment al;
               SmithWatAffine( bpaths[0], bpaths[1], al );
               align a(al);
               vec<int> mgg = a.MutationsGap1Gap2( bpaths[0], bpaths[1] );
               if ( mgg[0] == 1 && mgg[1] == 0 && mgg[2] == 0 )
               {    cout << "tame: substitution" << endl;
                    continue;    }
               if ( mgg[0] == 0 && a.Nblocks( ) == 2 )
               {    cout << "tame: indel" << endl;
                    continue;    }    }

          if ( paths.size( ) == 3 && bpaths[0].size( ) == bpaths[1].size( )
               && bpaths[1].size( ) == bpaths[2].size( ) )
          {    int count = 0;
               for ( int j = 0; j < bpaths[0].isize( ); j++ )
               {    if ( bpaths[0][j] != bpaths[1][j]
                         || bpaths[1][j] != bpaths[2][j] )
                    {    count++;    }    }
               if ( count == 2 )
               {    cout << "tame: two substitutions" << endl;
                    continue;    }    }

          Bool homopolymer = True;
          vec<Bool> hins( 4, False ); // homopolymer insertions
          cout << "\n" << paths.size( ) << " paths" << endl;
          for ( int i1 = 0; i1 < paths.isize( ); i1++ )
          for ( int i2 = i1 + 1; i2 < paths.isize( ); i2++ )
          {    cout << "\npath " << i1+1 << " vs " << i2+1 << endl;
               alignment al;
               SmithWatAffine( bpaths[i1], bpaths[i2], al );
               align a(al);

               const basevector &b1 = bpaths[i1], &b2 = bpaths[i2];
               int answer = 0, j, p1 = a.pos1( ), p2 = a.pos2( );
               int loc1 = p1, loc2 = p2;
               for ( j = 0; j < a.Nblocks( ); j++ )
               {    if ( a.Gaps(j) > 0 ) 
                    {    cout << loc1 << " insertion " << loc2 << " ";
                         for (int i = 0; i < a.Gaps(j); ++i)
	                      cout << as_base( bpaths[i2][p2+i] );
                         cout << endl;
                         for ( int i = 1; i < a.Gaps(j); i++ )
                              if ( b2[p2+i] != b2[p2] ) homopolymer = False;
                         hins[ b2[p2] ] = True;
                         p2 += a.Gaps(j);
                         loc2 += a.Gaps(j);    }
                    if ( a.Gaps(j) < 0 ) 
                    {    cout << loc1 << " ";
                         for ( int i = 0; i < -a.Gaps(j); ++i )
	                      cout << as_base( bpaths[i1][p1+i] );
                         cout << " "  << loc2 << " deletion" << endl;
                         for ( int i = 1; i < -a.Gaps(j); i++ )
                              if ( b1[p1+i] != b1[p1] ) homopolymer = False;
                         hins[ b1[p1] ] = True;
                         p1 -= a.Gaps(j);
                         loc1 -= a.Gaps(j);     }
                    for ( int x = 0; x < a.Lengths(j); x++ )
                    {    if ( bpaths[i1][p1] != bpaths[i2][p2] ) 
                         {    cout << loc1 << " " << as_base( bpaths[i1][p1] )
	                           << " " << loc2 << " " 
                                   << as_base( bpaths[i2][p2] ) << endl;
                              homopolymer = False;    }
                         ++p1; 
                         ++p2;
                         ++loc1;
                         ++loc2;    }    }

               PrintVisualAlignment( 
                    True, cout, bpaths[i1], bpaths[i2], a );    }

          if ( homopolymer && Sum(hins) == 1 )
          {    cout << "\nlooks like a homopolymer ambiguity" << endl;
               crush[i] = True;    }    }

     // Crush cells.  Creates nonsense.

     for ( int i = 0; i < cells.isize( ); i++ )
     {    if ( !crush[i] ) continue;
          he.DeleteEdges( ecells[i] );
          basevector b( 2 * ( shb.K( ) - 1 ) );
          const basevector& e = shb.EdgeObject( in[i] );
          const basevector& f = shb.EdgeObject( out[i] );
          for ( int j = 0; j < shb.K( ) - 1; j++ )
          {    b.Set( j, e[ e.isize( ) - shb.K( ) + j + 1 ] );
               b.Set( j + shb.K( ) - 1, f[j] );    }
          he.AddEdge( cells[i].front( ), cells[i].back( ), b );    }
     he.RemoveUnneededVertices( );
     he.RemoveDeadEdgeObjects( );
     edge_sizes.clear( );
     for ( int e = 0; e < he.EdgeObjectCount( ); e++ )
          edge_sizes.push_back( he.EdgeLengthKmers(e) );
     Sort(edge_sizes);
     cout << "\nN50 edge size after crushing = " << N50(edge_sizes) << endl;    }
