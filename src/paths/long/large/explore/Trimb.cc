///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// cd GapToy/52278.mouse5/a.final
// Subproject E=3628043,2101651,2101650,2101649,3169509 OUT_HEAD= ~/lcrd/sample
// GapToy READS=sample.fastb SAVE_60=True INSTANCE=99

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Improve60.h"

int main( )
{    RunTime( );

     // Define directory.

     String work_dir = "/wga/scr4/jaffe/GapToy/99";

     // Control.

     Bool verbose = True;
     String dir = work_dir + "/a.60";

     // Load assembly.

     HyperBasevector hb;
     BinaryReader::readFile( dir + "/a.hbv", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );
     ReadPathVec paths( dir + "/a.paths" );
     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     vecbasevector bases( work_dir + "/data/frag_reads_orig.fastb" );
     VecPQVec qualsp( work_dir + "/data/frag_reads_orig.qualp" );

     // Go through two passes.

     const int modes = 1; // FOR NOW JUST ONE PASS!
     for ( int mode = 1; mode <= modes; mode++ )
     {
          // Improve assembly.

          cout << "\nmode =  " << mode << endl;
          Bool verbose;
          if ( mode == 1 ) verbose = False;
          if ( mode == 2 ) verbose = True;
          Improve60( hb, inv, paths, bases, qualsp, verbose );

          // Write files.

          String dir2 = dir + "clean";
          if ( mode == 2 ) dir2 += "clean";
          Mkdir777(dir2);
          BinaryWriter::writeFile( dir2 + "/a.hbv", hb );
          BinaryWriter::writeFile( dir2 + "/a.inv", inv );
          paths.WriteAll( dir2 + "/a.paths" );
          HyperBasevectorX hbx(hb);
          BinaryWriter::writeFile( dir2 + "/a.hbx", hbx );
          invert( paths, paths_index, hb.EdgeObjectCount( ) );
          paths_index.WriteAll( dir2 + "/a.paths.inv" );    }

     // Delete small components having at most one read pair.

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     const int max_del = 100;
     vec<vec<int>> comp;
     hb.ComponentsE(comp);
     vec<int> dels;
     for ( int i = 0; i < comp.isize( ); i++ )
     {    Bool OK = False;
          vec<int64_t> pids;
          for ( int j = 0; j < comp[i].isize( ); j++ )
          {    int e = comp[i][j];
               int re = inv[e];
               for ( int l = 0; l < (int) paths_index[e].size( ); l++ )
                    pids.push_back( paths_index[e][l]/2 );
               for ( int l = 0; l < (int) paths_index[re].size( ); l++ )
                    pids.push_back( paths_index[re][l]/2 );
               if ( hb.Kmers(e) > max_del ) OK = True;    }
          UniqueSort(pids);
          if ( pids.size( ) > 1 ) OK = True;
          if ( !OK ) dels.append( comp[i] );    }
     hb.DeleteEdges(dels);
     Cleanup( hb, inv, paths );
     invert( paths, paths_index, hb.EdgeObjectCount( ) );

     // Again write files.

     {    String dir2 = dir + "clean";
          BinaryWriter::writeFile( dir2 + "/a.hbv", hb );
          BinaryWriter::writeFile( dir2 + "/a.inv", inv );
          paths.WriteAll( dir2 + "/a.paths" );
          HyperBasevectorX hbx(hb);
          BinaryWriter::writeFile( dir2 + "/a.hbx", hbx );
          invert( paths, paths_index, hb.EdgeObjectCount( ) );
          paths_index.WriteAll( dir2 + "/a.paths.inv" );
          vecbasevector edges( hb.E( ) );
          for ( int e = 0; e < hb.E( ); e++ )
               edges[e] = hb.O(e);
          edges.WriteAll( dir2 + "/a.fastb" );    }

     // Find unsatisfied pairs.

     hb.ToLeft(to_left), hb.ToRight(to_right);
     int64_t N = paths.size( );
     cout << "\n";
     int unsat = 0;
     for ( int i = 0; i < N/2; i++ )
     {    
          int64_t pid = i;
          const ReadPath &p1 = paths[2*pid], &p2 = paths[2*pid+1];
          if ( p1.size( ) == 0 || p2.size( ) == 0 ) continue;

          // if ( p1.size( ) == 1 && p2.size( ) == 1 && p1[0] == inv[ p2[0] ] )
          //      continue;

          vec<int> x2;
          for ( int j = 0; j < (int) p2.size( ); j++ )
               x2.push_back( inv[ p2[j] ] );
          x2.ReverseMe( );
          vec<int> y1, y2(x2);
          for ( int j = 0; j < (int) p1.size( ); j++ )
               y1.push_back( p1[j]);
          UniqueSort(y1), UniqueSort(y2);
          if ( Meet( y1, y2 ) ) continue;

          const int max_depth = 2;
          vec<int> es;
          es.push_back( p1.back( ) );
          Bool con = False;
          for ( int d = 1; d <= max_depth; d++ )
          {    vec<int> es2;
               for ( int j = 0; j < es.isize( ); j++ )
               {    int w = to_right[ es[j] ];
                    for ( int l = 0; l < hb.From(w).isize( ); l++ )
                         es2.push_back( hb.IFrom(w,l) );    }
               UniqueSort(es2);
               if ( BinMember( es2, inv[ p2.back( ) ] ) )
               {    con = True;
                    break;    }
               es = es2;    }
          if (con) continue;

          cout << "\n[" << i+1 << "=" << pid << "] "
               << printSeq(p1) << ".." << printSeq(x2) << endl;
          unsat++;    }
     cout << endl;

     // Find tree-ends.  A right tree-end is defined by an initial edge e: v --> w, 
     // having the property that e is the only edge entering w, that the graph
     // downstream of w is a tree, and that none of its edges are longer than a
     // fixed length mn = 200 kmers.  Moreover we only consider maximal right 
     // tree-ends.  Left tree-ends are defined symmetrically.

     const int mn = 200;
     cout << "\nfinding right tree-ends" << endl;
     vec<vec<int>> rte, lte;
     for ( int e = 0; e < hb.E( ); e++ )
     {    int w = to_right[e];
          vec<int> ws = {w}, es = {e};
          for ( int i = 0; i < ws.isize( ); i++ )
          {    int x = ws[i];
               if ( !hb.To(x).solo( ) ) goto next1;
               for ( int j = 0; j < hb.From(x).isize( ); j++ )
               {    int y = hb.From(x)[j], f = hb.IFrom(x,j);
                    if ( Member( ws, y ) || hb.Kmers(f) > mn ) goto next1;
                    ws.push_back(y), es.push_back(f);    }    }
          rte.push_back(es);
          next1: continue;    }
     vec<Bool> drte( rte.size( ), False );
     for ( int i1 = 0; i1 < rte.isize( ); i1++ )
     {    Bool maximal = True;
          for ( int i2 = 0; i2 < rte.isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               if ( Subset( rte[i1], rte[i2] ) )
               {    maximal = False;
                    break;    }    }
          Bool has_reads = False;
          for ( int j = 0; j < rte[i1].isize( ); j++ )
          {    if ( paths_index[ rte[i1][j] ].size( ) > 0
                    || paths_index[ inv[ rte[i1][j] ] ].size( ) > 0 )
               {    has_reads = True;    }    }
          if ( maximal && has_reads ) 
               cout << "right tree-end: " << printSeq( rte[i1] ) << endl;
          else drte[i1] = True;    }
     EraseIf( rte, drte );
     cout << "\nfinding left tree-ends" << endl;
     for ( int e = 0; e < hb.E( ); e++ )
     {    int w = to_left[e];
          vec<int> ws = {w}, es = {e};
          for ( int i = 0; i < ws.isize( ); i++ )
          {    int x = ws[i];
               if ( !hb.From(x).solo( ) ) goto next2;
               for ( int j = 0; j < hb.To(x).isize( ); j++ )
               {    int y = hb.To(x)[j], f = hb.ITo(x,j);
                    if ( Member( ws, y ) || hb.Kmers(f) > mn ) goto next2;
                    ws.push_back(y), es.push_back(f);    }    }
          lte.push_back(es);
          next2: continue;    }
     vec<Bool> dlte( lte.size( ), False );
     for ( int i1 = 0; i1 < lte.isize( ); i1++ )
     {    Bool maximal = True;
          for ( int i2 = 0; i2 < lte.isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               if ( Subset( lte[i1], lte[i2] ) )
               {    maximal = False;
                    break;    }    }
          Bool has_reads = False;
          for ( int j = 0; j < lte[i1].isize( ); j++ )
          {    if ( paths_index[ lte[i1][j] ].size( ) > 0
                    || paths_index[ inv[ lte[i1][j] ] ].size( ) > 0 )
               {    has_reads = True;    }    }
          if ( maximal && has_reads ) 
               cout << "left tree-end: " << printSeq( lte[i1] ) << endl;
          else dlte[i1] = True;    }
     EraseIf( lte, dlte );

     // Find links between tree ends.

     cout << "\nlinks between tree ends:\n";
     int count = 0;
     const int max_dist = 500;
     for ( int i1 = 0; i1 < rte.isize( ); i1++ )
     for ( int i2 = 0; i2 < lte.isize( ); i2++ )
     {    const vec<int> &x = rte[i1], &y = lte[i2];
          vec<int64_t> p1, p2;
          for ( int j = 0; j < x.isize( ); j++ )
          {    int e = x[j];
               for ( int l = 0; l < (int) paths_index[e].size( ); l++ )
               {    int64_t id = paths_index[e][l];
                    const ReadPath& q = paths[id];
                    int start = q.getOffset( );
                    for ( int m = 0; m < (int) q.size( ); m++ )
                    {    if ( q[m] == e )
                         {    if ( hb.Bases(e) - start - bases[id].isize( ) 
                                   <= max_dist )
                              {    p1.push_back(id);    }    }
                         start -= hb.Kmers(q[m]);    }    }    }
          for ( int j = 0; j < y.isize( ); j++ )
          {    int e = inv[ y[j] ];
               for ( int l = 0; l < (int) paths_index[e].size( ); l++ )
               {    int64_t id = paths_index[e][l];
                    const ReadPath& q = paths[id];
                    int start = q.getOffset( );
                    for ( int m = 0; m < (int) q.size( ); m++ )
                    {    if ( q[m] == e )
                         {    if ( hb.Bases(e) - start - bases[id].isize( ) 
                                   <= max_dist )
                              {    int64_t idx = ( id % 2 == 0 ? id+1 : id-1 );
                                   p2.push_back(idx);    }    }
                         start -= hb.Kmers(q[m]);    }    }    }
          UniqueSort(p1), UniqueSort(p2);
          vec<int64_t> p = Intersection( p1, p2 ), pp;
          for ( int j = 0; j < p.isize( ); j++ )
               pp.push_back( p[j]/2 );
          UniqueSort(pp);
          if ( pp.size( ) == 0 || rte[i1] == lte[i2] ) continue;
          cout << "[" << ++count << "] " << printSeq(rte[i1]) << " ==[" << pp.size( )
               << "]==> " << printSeq(lte[i2]) << " (" << pp[0];
          if ( pp.size( ) > 1 ) cout << "...";
          cout << ")" << endl;    }

     // Build read stacks associated to tree ends.

     for ( int i = 0; i < rte.isize( ); i++ )
     {    const vec<int>& z = rte[i];
          int e = z[0];
          int w = to_right[e];
          vec<int> ws = {w}, es = {e}, wpos = { hb.Kmers(e) };
          for ( int i = 0; i < ws.isize( ); i++ )
          {    int x = ws[i];
               for ( int j = 0; j < hb.From(x).isize( ); j++ )
               {    int y = hb.From(x)[j], f = hb.IFrom(x,j);
                    ws.push_back(y), es.push_back(f);    
                    wpos.push_back( wpos[i] + hb.Kmers(f) );    }    }

          // Compute locs = {(id,pos,fw)}, where pos is the start of read id, 
          // relative to the last base on edge e.

          vec< triple<int64_t,int,Bool> > locs; // (id,pos,fw)
          for ( int j = 0; j < es.isize( ); j++ )
          {    int f = es[j];
               int fpos = wpos[j];
               for ( int l = 0; l < (int) paths_index[f].size( ); l++ )
               {    int64_t id = paths_index[f][l];
                    const ReadPath& p = paths[id];
                    int offset = p.getOffset( );
                    for ( int m = 0; m < (int) p.size( ); m++ )
                    {    if ( p[m] == f )
                         {    int start = offset + fpos - hb.Kmers(f) - hb.Bases(e);
                              if ( start + 2 * bases[id].isize( ) >= 0 )
                                   locs.push( id, start, True );    }
                         offset -= hb.Kmers( p[m] );     }    }
               int rf = inv[f];
               for ( int l = 0; l < (int) paths_index[rf].size( ); l++ )
               {    int64_t id = paths_index[rf][l];
                    const ReadPath& p = paths[id];
                    int offset = hb.Bases(f) - p.getOffset( ) - bases[id].isize( );
                    for ( int m = 0; m < (int) p.size( ); m++ )
                    {    if ( p[m] == rf )
                         {    int start = offset + fpos - hb.Kmers(f) - hb.Bases(e);
                              if ( start + 2 * bases[id].isize( ) >= 0 )
                                   locs.push( id, start, False );    }
                         offset += hb.Kmers( p[m] );    }    }    }

          // Now build the stack.

          /*
          cout << "\nlocs for rte " << i << ":\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          for ( int j = 0; j < locs.isize( ); j++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          {    cout << locs[j].first << " @ " << locs[j].second // XXXXXXXXXXXXXXXXX
                    << " " << ( locs[j].third ? "fw" : "rc" ) << endl;    } // XXXXX
          */
          // Could happen legitimately:
          ForceAssert( !locs.empty( ) );
          int min_start = locs[0].second;
          int max_stop = locs[0].second + bases[ locs[0].first ].isize( );
          for ( int j = 1; j < locs.isize( ); j++ )
          {    min_start = Min( min_start, locs[j].second );
               max_stop = Max( max_stop, 
                    locs[j].second + bases[ locs[j].first ].isize( ) );    }
          int n = locs.size( ), k = max_stop - min_start;
          readstack stack( n, k );
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int64_t id = locs[i].first;
               qualvector q;
               qualsp[id].unpack(&q);
               int pos = locs[i].second - min_start;
               Bool fw = locs[i].third;
               stack.SetLen( i, bases[id].size( ) );
               stack.SetId( i, id );
               stack.SetRc2( i, !fw );
               if (fw)
               {    stack.SetOffset( i, pos );
                    for ( int j = 0; j < bases[id].isize( ); j++ )
                    {    int p = pos + j;
                         stack.SetBase( i, p, bases[id][j] );
                         stack.SetQual( i, p, q[j] );    }    }
               else
               {    stack.SetOffset( i, pos );
                    basevector br = bases[id];
                    br.ReverseComplement( );
                    qualvector qr = q;
                    qr.ReverseMe( );
                    for ( int j = 0; j < br.isize( ); j++ )
                    {    int p = pos + j;
                         stack.SetBase( i, p, br[j] );
                         stack.SetQual( i, p, qr[j] );    }    }    }

          cout << "\nstack for right tree-end " << i << endl;
          stack.Print(cout);    }

     // Summarize.

     cout << "\n" << hb.E( ) << " edges" << endl;
     cout << unsat << " unsatisfied links of " << bases.size( ) << " (" 
          << PERCENT_RATIO( 2, unsat, (int) bases.size( ) ) << ")" << endl << endl;

     Scram(0);    }
