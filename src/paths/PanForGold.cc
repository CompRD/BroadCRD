/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// PanForGold.  Look for short-fragment read pairs for which there does
/// not exist a closure (short enough path from one read to the other), but for
/// which there would exist such a path if kmers from some read were reclaimed.

#include <map>

#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmerRecord.h"
#include "MainTools.h"
#include "STLExtensions.h"
#include "graph/Digraph.h"
#include "kmer_freq/KAdjGraph.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "math/Functions.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"

template<int K> void ReportGaps( const vecbasevector& genome, 
     const vecbasevector& unibases, const digraph& A )
{    vec<unsigned int> sizes;
     for ( size_t i = 0; i < unibases.size( ); i++ )
          sizes.push_back( unibases[i].size( ) );
     Sort(sizes);
     cout << "\nN50 unibase = " << N50(sizes) << endl;
     vecbitvector cov;
     Mimic( genome, cov );
     vec< kmer_record< K+1, 2 > > G, U;
     basevector b(K+1), brc(K+1);
     kmer_record< K+1, 2 > rec;
     for ( size_t t = 0; t < genome.size( ); t++ )
     {    for ( int j = 0; j <= genome[t].isize( ) - (K+1); j++ )
          {    b.SetToSubOf( genome[t], j, K+1 );
               rec.Set( b, t, j );
               G.push_back(rec);    }    }
     for ( size_t u = 0; u < unibases.size( ); u++ )
     {    for ( int j = 0; j <= unibases[u].isize( ) - (K+1); j++ )
          {    b.SetToSubOf( unibases[u], j, K+1 );
               brc.ReverseComplement(b);
               rec.Set( b, 0, 0 );
               U.push_back(rec);
               rec.Set( brc, 0, 0 );
               U.push_back(rec);    }
          for ( int j = 0; j < A.From(u).isize( ); j++ )
          {    int v = A.From(u)[j];
               for ( int l = 0; l < K; l++ )
                    b.Set( l, unibases[u][ unibases[u].isize( ) - K + l ] );
               b.Set( K, unibases[v][K-1] );
               brc.ReverseComplement(b);
               rec.Set( b, 0, 0 );
               U.push_back(rec);
               rec.Set( brc, 0, 0 );
               U.push_back(rec);    }     }
     Sort(G), Sort(U);
     longlong j = 0;
     for ( longlong i = 0; i < (longlong) G.size( ); i++ )
     {    nextj:
          int c = U[j].CmpKmers( G[i] );
          if ( c > 0 ) continue;
          if ( c == 0 ) cov[ G[i].GetId( ) ].Set( G[i].GetPos( ), True );
          else 
          {    ++j;
               if ( j == (longlong) U.size( ) ) break;
               goto nextj;    }    }
     int gaps = 0;
     cout << "\nGaps:\n";
     for ( size_t t = 0; t < genome.size( ); t++ )
     {    for ( unsigned int i = 0; i < genome[t].size( ); i++ )
          {    if ( cov[t][i] ) continue;
               unsigned int j;
               for ( j = i + 1; j < genome[t].size( ); j++ )
                    if ( cov[t][j] ) break;
               cout << t << "." << i << "-" << j << " [" << j-i << "]\n";
               ++gaps;
               i = j - 1;    }    }
     cout << "\n" << gaps << " gaps" << endl << endl;    }

template<int K> void MapPairs( const vec<read_pairing>& pairs, 
     const vec<Bool>& pairs_to_map, const vecbasevector& reads, 
     const vecbasevector& unibases, const vec<int>& to_rc, 
     const HyperBasevector& hb,
     vec<int>& U1, vec<int>& U2, vec<int>& P1, vec<int>& P2 )
{    
     // Initialize answers.

     U1.resize_and_set( pairs.size( ), -1 ), U2.resize_and_set( pairs.size( ), -1 );
     P1.resize_and_set( pairs.size( ), -1 ), P2.resize_and_set( pairs.size( ), -1 );

     // Set up to find the start point of each read in the unibases.

     int nreads = reads.size( );
     vec<int> start_u(nreads, -1), start_p(nreads, -1);

     // Determine which unibases can follow which.  (This could be read off hb.)

     int nuni = unibases.size( );
     vec< vec<int> > nexts(nuni);
     {    vec< kmer_record<K-1,1> > unistarts, unistops;
          unistarts.reserve(nuni), unistops.reserve(nuni);
          for ( size_t i = 0; i < unibases.size( ); i++ )
          {    static basevector b;
               b.SetToSubOf( unibases[i], 0, K-1 );
               static kmer_record<K-1,1> rec;
               rec.Set( b, i, 0 );
               unistarts.push_back(rec);
               int n = unibases[i].size( );
               b.SetToSubOf( unibases[i], n - (K-1), K-1 );
               rec.Set( b, i, n - (K-1) );
               unistops.push_back(rec);    }
          Sort(unistarts), Sort(unistops);
          int ulast = 0;
          for ( int i = 0; i < unistops.isize( ); i++ )
          {    int j;
               Bool eq = False;
               for ( j = ulast; j < unistarts.isize( ); j++ )
               {    if ( unistarts[j] > unistops[i] ) break;
                    if ( unistarts[j].EqualKmers( unistops[i] ) )
                    {    eq = True;
                         break;    }    }
               if ( j == unistarts.isize( ) ) break;
               if (eq)
               {    for ( int z = j; z < unistarts.isize( ); z++ )
                    {    if ( !unistarts[z].EqualKmers( unistops[i] ) ) break;
                         nexts[ unistops[i].GetId( ) ].
                              push_back( unistarts[z].GetId( ) );    }    }
               ulast = j;    }    }

     // Form kmer records associated with first kmers in the reads.

     vec< kmer_record<K,1> > readstarts, unistarts;
     readstarts.reserve(nreads);
     for ( size_t j = 0; j < reads.size( ); j++ )
     {    if ( reads[j].isize( ) >= K )
          {    static basevector b;
               b.SetToSubOf( reads[j], 0, K );
               static kmer_record<K,1> rec;
               rec.Set( b, j, 0 );
               readstarts.push_back(rec);    }    }

     // Form kmer records associated with kmers in the unibases.

     int nrecs = 0;
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    if ( unibases[i].isize( ) >= K )
               nrecs += ( unibases[i].isize( ) - K + 1 );    }
     unistarts.reserve(nrecs);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    for ( int j = 0; j <= unibases[i].isize( ) - K; j++ )
          {    static basevector b;
               b.SetToSubOf( unibases[i], j, K );
               static kmer_record<K,1> rec;
               rec.Set( b, i, j );
               unistarts.push_back(rec);    }    }

     // Match them up.

     Sort(readstarts), Sort(unistarts);
     int ulast = 0;
     for ( int i = 0; i < readstarts.isize( ); i++ )
     {    int j;
          Bool eq = False;
          for ( j = ulast; j < unistarts.isize( ); j++ )
          {    if ( unistarts[j] > readstarts[i] ) break;
               if ( unistarts[j].EqualKmers( readstarts[i] ) )
               {    eq = True;
                    break;    }    }
          if ( j == unistarts.isize( ) ) break;
          if (eq)
          {    int id = readstarts[i].GetId( );
               start_u[id] = unistarts[j].GetId( );
               start_p[id] = unistarts[j].GetPos( );    }
          ulast = j;    }

     // Create index mapping edges to vertices.  Set up integer graph.

     vec<int> L;
     vec<int> to_left, to_right;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          L.push_back( hb.EdgeObject(i).isize( ) - (K-1) );
     digraphE<int> G( hb, L );
     G.ToLeft(to_left), G.ToRight(to_right);

     // For each short-fragment read pair, map the two reads as far as possible           // along unipaths, then deduce an implied separation between two unipaths 
     // u1 and u2.

     const int max_short_fragment = 1000;
     for ( int pi = 0; pi < pairs.isize( ); pi++ )
     {    if ( pairs[pi].sep > max_short_fragment ) continue;
          int id1 = pairs[pi].id1, id2 = pairs[pi].id2;
          vec<int> U(2, -1), P(2, -1);
          for ( int side = 0; side < 2; side++ )
          {    int j = ( side == 0 ? id1 : id2 );
               int N = reads[j].size( );
               int u = start_u[j], p = start_p[j];
               U[side] = u;
               P[side] = -p;

               // Check for case where first kmer of read can't be placed.

               if ( u < 0 ) break; 

               int up = to_rc[u];
               int nu = unibases[u].size( );
               ++p;
               for ( int v = K; v < reads[j].isize( ); v++ )
               {    char next = reads[j][v];
                    if ( p + K <= unibases[u].isize( ) )
                    {    char n = unibases[u][p+K-1];
                         if ( n != next )
                         {
                              // We've aligned the last alignable base of the read.
                              // Define the read base that would map to the first
                              // base of unipath u.

                              U[side] = u;
                              P[side] = v - (p+K-1);
                              break;    }

                         ++p;    }
                    else
                    {    int z, ns = nexts[u].size( );
                         for ( z = 0; z < ns; z++ )
                         {    if ( next == unibases[ nexts[u][z] ][K-1] )
                              {    u = nexts[u][z];
                                   up = to_rc[u];
                                   nu = unibases[u].size( );
                                   U[side] = u;
                                   P[side] = v;
                                   p = 1;
                                   break;    }    }
                         if ( z == ns ) 
                         {    
                              // We've aligned the last alignable base of the read.
                              // Define the read base that would map to the first
                              // base of unipath u.

                              U[side] = u;
                              P[side] = v - (p+K-1);
                              break;    }    }    }    }

          // Record answers.

          U1[pi] = U[0], U2[pi] = U[1], P1[pi] = P[0], P2[pi] = P[1];    }    }

template<int K> void Pan( const vec<read_pairing>& pairs, 
     const vecbasevector& reads, vecbasevector& unibases, digraph& A,
     const vec<int>& to_rc, const HyperBasevector& hb,
     const vec<look_align>& aligns, const vec< vec<int> >& aligns_index )
{    
     // Map pairs.
     
     vec<int> U1, U2, P1, P2;
     vec<Bool> pairs_to_map( pairs.size( ), True );
     MapPairs<K>( pairs, pairs_to_map, reads, unibases, to_rc, hb, U1, U2, P1, P2 );

     // Create index mapping edges to vertices.  Set up integer graph.

     vec<int> L;
     vec<int> to_left, to_right;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          L.push_back( hb.EdgeObject(i).isize( ) - (K-1) );
     digraphE<int> G( hb, L );
     G.ToLeft(to_left), G.ToRight(to_right);

     // For each short-fragment read pair, map the two reads as far as possible           // along unipaths, then deduce an implied separation between two unipaths 
     // u1 and u2.

     vec< kmer_record<K,1> > allgapkmers;
     const int max_short_fragment = 1000;
     vec<Bool> interesting( pairs.size( ), False );
     for ( int pi = 0; pi < pairs.isize( ); pi++ )
     {    if ( pairs[pi].sep > max_short_fragment ) continue;
          int id1 = pairs[pi].id1, id2 = pairs[pi].id2;
          vec<int> U(2), P(2);
          U[0] = U1[pi], U[1] = U2[pi], P[0] = P1[pi], P[1] = P2[pi];
          
          // Reject pairs that appear to be duplicate reads.

          if ( U[0] == U[1] && P[0] == P[1] ) continue;

          // Compute predicted separation in K-mers between the end of u1 and the 
          // beginning of u2.

          if ( U[0] < 0 || U[1] < 0 ) continue;
          int u1 = U[0], u2 = to_rc[ U[1] ];
          int nu1 = unibases[u1].size( ), nu2 = unibases[u2].size( );
          int sep = pairs[pi].sep - ( nu1 - reads[id1].isize( ) + P[0] )
               - ( nu2 - reads[id2].isize( ) + P[1] ) - (K-1);

          // Find paths from u1 to u2.

          int v = to_right[u1], w = to_left[u2];

          int dev = pairs[pi].sd;
          const int max_devs = 6, max_devs_high = 8;
          Bool good = False, exploded = False;

          if ( u1 == u2 )
          {    int SEP = sep + unibases[u1].isize( );
               double err = Abs( double(SEP)/double(dev) );
               if ( err < (double) max_devs_high ) good = True;
               // not fully implemented yet; have to allow jump case
	  }
          else if ( v == w )
          {    double err = Abs( double(sep)/double(dev) );
               if ( err < (double) max_devs_high ) good = True;
	  }
          else if ( u1 != u2 )
          {    static vec< vec<int> > paths;
               const int max_partials = 10000;
               exploded = !G.AllPathsLengthRangeAlt( v, w, 0, 
                    sep + max_devs * pairs[pi].sd, to_right, paths, 1, 0, True, 
                    True, max_partials );
               if ( paths.nonempty( ) && !exploded ) good = True;    
	  }
          if ( good || exploded ) continue;

          // Build left side of pair.

          const int max_side = 20;
          int max_extend = ( sep + max_devs * pairs[pi].sd ) / 2;
          vec< triple<int,int,Bool> > vplus;  // (edge, dist[stop,v], processed)
          vplus.push( u1, 0, 0 >= max_extend );
          Bool side_too_big = False;
          while(1)
          {    Bool changed = False;
               for ( int j = 0; j < vplus.isize( ); j++ )
               {    if ( vplus[j].third ) continue;
                    int x = to_right[ vplus[j].first ];
                    for ( int l = 0; l < hb.From(x).isize( ); l++ )
                    {    int e = hb.EdgeObjectIndexByIndexFrom( x, l );
                         int d = vplus[j].second + L[e];
                         Bool found = False;
                         for ( int m = 0; m < vplus.isize( ); m++ )
                         {    if ( e != vplus[m].first ) continue;
                              found = True;
                              if ( d >= vplus[m].second ) continue;
                              changed = True;
                              vplus[m].second = d;
                              vplus[m].third = (d >= max_extend);    }
                         if ( !found )
                         {    vplus.push( e, d, (d >= max_extend) );
                              if ( vplus.isize( ) > max_side )
                              {    side_too_big = True;
                                   break;    }
                              changed = True;    }
                         if (side_too_big) break;    }
                    if (side_too_big) break;    }
               if ( !changed || side_too_big ) break;    }
          if (side_too_big) continue;

          // Build right side of pair.

          vec< triple<int,int,Bool> > wminus;  // (edge, dist[start,w], processed)
          wminus.push( u2, 0, 0 >= max_extend );
          while(1)
          {    Bool changed = False;
               for ( int j = 0; j < wminus.isize( ); j++ )
               {    if ( wminus[j].third ) continue;
                    int x = to_left[ wminus[j].first ];
                    for ( int l = 0; l < hb.To(x).isize( ); l++ )
                    {    int e = hb.EdgeObjectIndexByIndexTo( x, l );
                         int d = wminus[j].second + L[e];
                         Bool found = False;
                         for ( int m = 0; m < wminus.isize( ); m++ )
                         {    if ( e != wminus[m].first ) continue;
                              found = True;
                              if ( d >= wminus[m].second ) continue;
                              changed = True;
                              wminus[m].second = d;
                              wminus[m].third = (d >= max_extend);    }
                         if ( !found )
                         {    wminus.push( e, d, (d >= max_extend) );
                              if ( wminus.isize( ) > max_side )
                              {    side_too_big = True;
                                   break;    }
                              changed = True;    }
                         if (side_too_big) break;    }
                    if (side_too_big) break;    }
               if ( !changed || side_too_big ) break;    }
          if (side_too_big) continue;
          interesting[pi] = True;

          // Find kmer records in the two sides.

          const int tailpart = 30;
          vec< kmer_record<K,1> > leftkmers, rightkmers;
          for ( int j = 0; j < vplus.isize( ); j++ )
          {    const basevector& u = unibases[ vplus[j].first ];
               static basevector b, brc;
               static kmer_record<K,1> rec;
               int stop = vplus[j].second;
               if ( stop == 0 )
               {    for ( int m = Max( 0, u.isize( ) - tailpart - K ); 
                         m <= u.isize( ) - K; m++ )
                    {    b.SetToSubOf( u, m, K );
                         brc.ReverseComplement(b);
                         if ( brc < b || b == brc ) rec.Set( brc, 0, 0 );
                         else rec.Set( b, 0, 0 );
                         leftkmers.push_back(rec);    }    }
               else
               {    int start = stop - ( u.isize( ) - (K-1) );
                    for ( int m = 0; m <= u.isize( ) - K; m++ )
                    {    if ( start + m > max_extend ) break;
                         b.SetToSubOf( u, m, K );
                         brc.ReverseComplement(b);
                         if ( brc < b || b == brc ) rec.Set( brc, 0, 0 );
                         else rec.Set( b, 0, 0 );
                         leftkmers.push_back(rec);    }    }    }
          for ( int j = 0; j < wminus.isize( ); j++ )
          {    const basevector& u = unibases[ wminus[j].first ];
               static basevector b, brc;
               static kmer_record<K,1> rec;
               int start = wminus[j].second;
               if ( start == 0 )
               {    for ( int m = 0; m <= Min( tailpart - 1, u.isize( ) - K ); m++ )
                    {    b.SetToSubOf( u, m, K );
                         brc.ReverseComplement(b);
                         if ( brc < b || b == brc ) rec.Set( brc, 0, 0 );
                         else rec.Set( b, 0, 0 );
                         rightkmers.push_back(rec);    }    }
               else
               {    int stop = start + ( u.isize( ) - (K-1) );
                    for ( int m = u.isize( ) - K; m >= 0; m-- )
                    {    if ( stop + m > max_extend ) break;
                         b.SetToSubOf( u, m, K );
                         brc.ReverseComplement(b);
                         if ( brc < b || b == brc ) rec.Set( brc, 0, 0 );
                         else rec.Set( b, 0, 0 );
                         rightkmers.push_back(rec);    }    }    }
          allgapkmers.append(leftkmers), allgapkmers.append(rightkmers);

          // Find unibases common to both sides.
     }
     // Find reads containing gap kmers.

     UniqueSort(allgapkmers);
     cout << Date( ) << ": sorting allgapkmers" << endl;
     PRINT( allgapkmers.size( ) );
     cout << Date( ) << ": going through reads" << endl;
     // TODO: potentially dangerous truncation of index by matches
     vec<int> matches;
     for ( size_t id = 0; id < reads.size( ); id++ )
     {    if ( id % 1000000 == 0 ) DPRINT2( id, reads.size( ) );
          for ( unsigned int j = 0; j <= reads[id].size( ) - K; j++ )
          {    static basevector b, brc;
               static kmer_record<K,1> rec;
               b.SetToSubOf( reads[id], j, K );
               brc.ReverseComplement(b);
               if ( brc < b || b == brc ) rec.Set( brc, 0, 0 );
               else rec.Set( b, 0, 0 );
               if ( BinMember( allgapkmers, rec ) )
               {    matches.push_back(id);
                    break;    }    }    }
     cout << "\n";
     PRINT( matches.size( ) );

     // Rebuild unibases.

     vecbasevector unibases_plus(unibases);
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int i1 = 0; i1 < hb.To(v).isize( ); i1++ )
          {    int e1 = hb.EdgeObjectIndexByIndexTo( v, i1 );
               for ( int i2 = 0; i2 < hb.From(v).isize( ); i2++ )
               {    int e2 = hb.EdgeObjectIndexByIndexFrom( v, i2 );
                    static basevector b1, b2, b;
                    b1 = hb.EdgeObject(e1), b2 = hb.EdgeObject(e2);
                    b1.resize( b1.size( ) - (K-1) );
                    b = Cat( b1, b2 );
                    unibases_plus.push_back_reserve(b);    }    }    }
     vecbasevector unibases_plus0(unibases_plus);
     for ( int i = 0; i < matches.isize( ); i++ )
          unibases_plus.push_back_reserve( reads[ matches[i] ] );
     vecKmerPath paths, pathsrc, unipaths2;
     vec<tagged_rpint> pathsdb, unipathsdb2;
     cout << Date( ) << ": calling ReadsToPaths" << endl;
     ReadsToPathsCoreY( unibases_plus, K, paths, pathsrc, pathsdb );
     cout << Date( ) << ": unipathing" << endl;
     Unipath( paths, pathsrc, pathsdb, unipaths2, unipathsdb2 );
     digraph A2;
     BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths2,
          unipathsdb2, A2 );
     HyperKmerPath h2;
     BuildUnipathAdjacencyHyperKmerPath( K, A2, unipaths2, h2 );
     KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, unibases_plus );
     vecbasevector unibases2;
     for ( size_t i = 0; i < unipaths2.size( ); i++ )
          unibases2.push_back_reserve( kbb.Seq( unipaths2[i] ) );
     cout << Date( ) << ": building hb2" << endl;
     HyperBasevector hb2( h2, kbb );

     // Remap the interesting pairs.

     vec<int> to_rc2;
     UnibaseInvolution( unibases2, to_rc2 );
     cout << Date( ) << ": remapping pairs" << endl; 
     MapPairs<K>( 
          pairs, interesting, reads, unibases2, to_rc2, hb2, U1, U2, P1, P2 );

     // Create index mapping edges to vertices.  Set up integer graph.

     vec<int> L2, to_left2, to_right2;
     for ( int i = 0; i < hb2.EdgeObjectCount( ); i++ )
          L2.push_back( hb2.EdgeObject(i).isize( ) - (K-1) );
     digraphE<int> G2( hb2, L2 );
     G2.ToLeft(to_left2), G2.ToRight(to_right2);

     // For each short-fragment read pair, map the two reads as far as possible           // along unipaths, then deduce an implied separation between two unipaths 
     // u1 and u2.

     cout << "\nNOW TRYING TO CLOSE UNCLOSED PAIRS\n" << endl;
     vecbasevector extras;
     for ( int pi = 0; pi < pairs.isize( ); pi++ )
     {    if ( !interesting[pi] ) continue;
          int id1 = pairs[pi].id1, id2 = pairs[pi].id2;
          vec<int> U(2), P(2);
          U[0] = U1[pi], U[1] = U2[pi], P[0] = P1[pi], P[1] = P2[pi];

          // Reject pairs that appear to be duplicate reads.

          if ( U[0] == U[1] && P[0] == P[1] ) continue;

          // Compute predicted separation in K-mers between the end of u1 and the 
          // beginning of u2.

          if ( U[0] < 0 || U[1] < 0 ) continue;
          int u1 = U[0], u2 = to_rc2[ U[1] ];
          int nu1 = unibases2[u1].size( ), nu2 = unibases2[u2].size( );
          int sep = pairs[pi].sep - ( nu1 - reads[id1].isize( ) + P[0] )
               - ( nu2 - reads[id2].isize( ) + P[1] ) - (K-1);

          // Find paths from u1 to u2.

          int v = to_right2[u1], w = to_left2[u2];
          int dev = pairs[pi].sd;
          const int max_devs = 6, max_devs_high = 8;
          Bool good = False, exploded = False;
          //cout << "\n";
          //PRINT2( pi, pairs.size( ) );
          //cout << BaseAlpha(u1) << "[" << nu1 << "] --- (" << sep << " +/- " 
          //     << dev << ") ---> " << BaseAlpha(u2) << "[" << nu2 << "]" << endl;
          //cout << align_report[pi] << endl;
          if ( u1 == u2 )
          {    int SEP = sep + unibases2[u1].isize( );
               double err = Abs( double(SEP)/double(dev) );
               if ( err < (double) max_devs_high ) 
               {    good = True;
                    cout << "GOOD!\n";    }
               //cout << "same edge, sep = " << SEP << ", off by "
               //     << ToString( err, 1 ) << " devs" << endl;
               // not fully implemented yet; have to allow jump case
                      }
          else if ( v == w )
          {    double err = Abs( double(sep)/double(dev) );
               if ( err < (double) max_devs_high ) 
               {    good = True;
                    cout << "GOOD!\n";    }
               //cout << "edges adjacent, sep = " << sep << ", off by "
               //     << ToString( err, 1 ) << " devs" << endl;
	  }
          else if ( u1 != u2 )
          {    static vec< vec<int> > paths;
               const int max_partials = 10000, max_paths = 20;
               exploded = !G2.AllPathsLengthRangeAlt( v, w, 0, 
                    sep + max_devs * pairs[pi].sd, to_right2, paths, max_paths, 
                    0, True, True, max_partials );
               if ( paths.nonempty( ) && !exploded ) 
               {    good = True;    
		 //cout << "GOOD!\n";
                    if ( paths.isize( ) < max_paths )
		      {    //cout << "Saving extras!\n";
                         for ( int x = 0; x < paths.isize( ); x++ )
                         {    basevector extra = unibases2[u1];
                              extra.resize( extra.size( ) - (K-1) );
                              for ( int j = 0; j < paths[x].isize( ); j++ )
                              {    extra = Cat( extra, unibases2[ paths[x][j] ] );
                                   extra.resize( extra.size( ) - (K-1) );    }
                              extra = Cat( extra, unibases2[u2] );
                              extras.push_back_reserve(extra);    }    }    }
	  }
     }
     PRINT( extras.size( ) );

     // Build new unipaths, adding in extras.

     unibases_plus0.Append(extras);
     vecKmerPath unipaths3;
     vec<tagged_rpint> unipathsdb3;
     cout << Date( ) << ": calling ReadsToPaths" << endl;
     ReadsToPathsCoreY( unibases_plus0, K, paths, pathsrc, pathsdb );
     cout << Date( ) << ": unipathing" << endl;
     Unipath( paths, pathsrc, pathsdb, unipaths3, unipathsdb3 );
     digraph A3;
     BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths3,
          unipathsdb3, A3 );
     HyperKmerPath h3;
     BuildUnipathAdjacencyHyperKmerPath( K, A3, unipaths3, h3 );
     KmerBaseBroker kbb3( K, paths, pathsrc, pathsdb, unibases_plus0 );
     vecbasevector unibases3;
     for ( size_t i = 0; i < unipaths3.size( ); i++ )
          unibases3.push_back_reserve( kbb3.Seq( unipaths3[i] ) );
     cout << Date( ) << ": building hb3" << endl;
     HyperBasevector hb3( h3, kbb3 );
     unibases = unibases3;
     A = A3;    }

int main( int argc, char *argv[] )
{
     RunTime( ); BeginCommandArguments;
     CommandArgument_String(PRE);   
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(READS, "reads_orig");
     CommandArgument_String_OrDefault(UNIBASES_IN, "trusted_unibases.pan2");
     CommandArgument_String_OrDefault(UNIGRAPH_IN, "trusted_unipathGraph.pan2");
     CommandArgument_String_OrDefault(UNIBASES_OUT, 
          "trusted_unibases.pan3");
     CommandArgument_String_OrDefault(UNIGRAPH_OUT, 
          "trusted_unipathGraph.pan3");
     CommandArgument_Bool_OrDefault(USE_GENOME, True);
     CommandArgument_Bool_OrDefault(WRITE, True);
     EndCommandArguments;

     cout << Date( ) << ": Loading reads" << endl;
     // Define directories.

     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;

     // Load the reads.

     vecbasevector reads( run_dir + "/" + READS + ".fastb" );
     cout << Date( ) << ": reads loaded" << endl;

     // Load genome and compute its size.

     vecbasevector genome;

     // Load the unibases and their graph.

     vecbasevector unibases;
     String KS = ToString(K);
     unibases.ReadAll( run_dir + "/" + READS + "." + UNIBASES_IN + ".k" 
          + KS + ".fastb" );
     digraph unigraph;
     BinaryReader::readFile( run_dir + "/" + READS + "." + UNIGRAPH_IN + ".k" + KS, &unigraph );

     // Build the HyperKmerPath.

     digraph& A = unigraph;
     HyperBasevector hb;
     BuildUnibaseAdjacencyHyperBasevector( K, A, unibases, hb );

     // Find perfect alignments, build index.
          
     vec<look_align> aligns;
     vec< vec<int> > aligns_index( unibases.size( ) );
     if (USE_GENOME)
     {    ForceAssertLe( 12, K );
          PerfectLookup( 12, unibases, data_dir + "/genome.lookup", 
               aligns, FW_OR_RC );
          for ( int i = 0; i < aligns.isize( ); i++ )
               aligns_index[ aligns[i].query_id ].push_back(i);    }

     // Report gaps.

     if (USE_GENOME)
     {    if ( K == 20 ) ReportGaps<20>( genome, unibases, A );    }

     // For each read pair, ...

     vec<read_pairing> pairs;
     ReadPairsFile( run_dir + "/" + READS + ".pairto", pairs );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc );
     if ( K == 20 )
     {    Pan<20>( pairs, reads, unibases, A, to_rc, hb, aligns, 
               aligns_index );    }
     else
     {    cout << "K = " << K << " not implemented" << endl;
          cout << "Abort." << endl;
          exit(1);    }

     // Write out new unibases and graph.

     if (WRITE)
     {    cout << Date( ) << ": Writing unibases and unibase adjacency graph files" 
               << endl;
          unibases.WriteAll( run_dir + "/" + READS + "." + UNIBASES_OUT + ".k" 
               + KS + ".fastb" );
	  BinaryWriter writer( ( run_dir + "/" + READS + "." + UNIGRAPH_OUT + ".k"  + KS ).c_str( ) );
	  writer.write( unigraph );
	  writer.close( );
     }

     // Report gaps.

     if (USE_GENOME)
     {    if ( K == 20 ) ReportGaps<20>( genome, unibases, A );    }
     cout << Date( ) << ": DONE!" << endl;    }
