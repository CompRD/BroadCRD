///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "math/Functions.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/PlaceReads1.h"
#include "paths/long/ReadPath.h"
#include "paths/long/SupportedHyperBasevector.h"

namespace {
    struct CompareReadPlaceByQsum1 {
        bool operator() (const read_place& a, const read_place& b) 
        { return a.Qsum() < b.Qsum(); }
    };
}

typedef priority_queue< 
     read_place, std::vector<read_place>, CompareReadPlaceByQsum1 > pqsum;

void ExtendReadPlaces1( vec<read_place>& places_part, pqsum& candidates,
     const basevector& b, const qualvector& q, const HyperBasevector& hb, 
     const vec<int>& to_right, const int max_partials,
     const int min_qual = 3, const double prox = 0 )
{
     // First score partials for first K bases.  Pick the best.

     if ( places_part.empty( ) ) return;
     int K = hb.K( );
     vec<int> qsum0( places_part.size( ), 0 );
     vec<int> ids( places_part.size( ), vec<int>::IDENTITY );
     for ( int i = 0; i < places_part.isize( ); i++ ) 
     {    const read_place& p = places_part[i];
          int e = p.E(0), s = p.P( );
          for ( int j = 0; j < K; j++ )
	      if ( b[j] != hb.EdgeObject(e)[s+j] ) qsum0[i] += q[j];    }

     // Pick best placement based on qsum

     SortSync( qsum0, ids );
     if ( places_part.size( ) > 1 && qsum0[0] == qsum0[1] ) return;
     const int max_qsum_K = 50;
     const double max_qsum_ratio = 0.3;
     if ( qsum0[0] > max_qsum_K ) return;
     if ( places_part.size( ) > 1 && double(qsum0[0])/double(qsum0[1]) 
          > max_qsum_ratio )
     {    return;    }
     read_place p = places_part[ ids[0] ];

     // Extend to end of edge or read

     int e = p.E(0), s = p.P( );
     for ( int j = K; j < b.isize( ); j++ ) 
     {    if ( s+j == hb.EdgeObject(e).isize( ) ) break;   //  reached end of edge
          if ( b[j] != hb.EdgeObject(e)[s+j] ) qsum0[0] += q[j];    // add to qsum
               }
     p.SetQsum( qsum0[0] ); 
   
     // Stop if read is contained in a single edge
     
     int bpos = hb.EdgeObject(e).isize( ) - s;
     if ( bpos >= b.isize( ) ) 
     {    candidates.push(p);
          return;    }
     
     // Extend along graph

     while(1) 
     {    int v = to_right[ p.E( p.N( ) - 1 ) ];  // get next edge

          if ( hb.From(v).empty( ) )  
          {    // No more edges, stop
	       candidates.push(p);
	       return;   }
       
          // What is the shortest edge, or remaining read bases if smaller

          int min_ext = b.isize( ) - bpos;
          for ( int j = 0; j < hb.From(v).isize( ); j++ ) 
          {    int e = hb.EdgeObjectIndexByIndexFrom( v, j );
	       min_ext = Min( min_ext, hb.EdgeLengthKmers(e) );    }

          // Extend along the shortest edge (and same distance along other edges)

          vec<int> qsum( hb.From(v).isize( ), 0 );
          vec<int> ids( hb.From(v).isize( ), vec<int>::IDENTITY );
          for ( int j = 0; j < hb.From(v).isize( ); j++ ) 
          {    int e = hb.EdgeObjectIndexByIndexFrom( v, j );
	       for ( int l = 0; l < min_ext; l++ ) 
               {    if ( b[bpos+l] != hb.EdgeObject(e)[K-1+l] ) 
                         qsum[j] += q[bpos+l];    }    }

          // Pick best edge based on qsum

          SortSync( qsum, ids );
          if ( qsum.size( ) > 1 && qsum[0] == qsum[1] ) 
          {    candidates.push(p);
	       return;    }
          p.SetQsum( p.Qsum( ) + qsum[0] );

          // Extend along remainder of the current edge

          int e = hb.EdgeObjectIndexByIndexFrom( v, ids[0] );
          for ( int l = min_ext; l < b.isize( ) - bpos; l++ ) 
          {    if ( K-1+l == hb.EdgeObject(e).isize( ) ) break;
	       if ( b[bpos+l] != hb.EdgeObject(e)[K-1+l] ) 
 	            p.SetQsum( p.Qsum( ) + q[bpos+l] );    }
          p.AddEdge(e);

          // Stop if we have no more of the read to align.

          bpos += hb.EdgeLengthKmers(e);
          if ( bpos >= b.isize( ) ) 
          {    candidates.push(p);
	       return;    }    }    }

void PlaceReads1( const HyperBasevector& hb, const vecbasevector& bases,
     const vecqualvector& quals, ReadPathVec& paths )
{
     // double clock1 = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Define heuristics.

     const int DIVINE_MAX_LOCS = 1000;
     const int max_partials = 100000;

     // Build data structures.

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Find alignments seeded on 60-mers.

     paths.resize(0);
     paths.resize( bases.size( ) );
     vec<Bool> found( bases.size( ), False );
     vec< vec<read_place> > partials( bases.size( ) );

     const int K0 = 60;
     ForceAssertGe( hb.K( ), K0 );
     vecbasevector all;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          all.push_back( hb.EdgeObject(i) );

     for ( int i = 0; i < (int) bases.size( ); i++ )
     {    basevector b = bases[i];
          b.resize(K0);
          all.push_back(b);    }
     // cout << TimeSince(clock1) << " used in initial setup" << endl; // XXXXXXXXXXX

     // double clock2 = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     vec< triple<kmer<K0>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );
     // cout << TimeSince(clock2) << " used making lookup" << endl; // XXXXXXXXXXXXXX

     // double clock3 = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int kpass = 0; kpass < 5; kpass++ )
     {    int K = 60;
     
          // Note that these numbers must all be divisible by four.

          //if ( kpass == 0 ) K = 60;
          if ( kpass == 1 ) K = 40;
          if ( kpass == 2 ) K = 20;
          if ( kpass == 3 ) K = 16;
          if ( kpass == 4 ) K = 12;

          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               {    if ( K == 60 )
                    {    if ( kmers_plus[j].first != kmers_plus[i].first ) 
                              break;    }
                    else
                    {    Bool diff = False;
                         for ( int l = 0; l < K/4; l++ )
                         {    if ( kmers_plus[j].first.Bytes( )[l]
                                   != kmers_plus[i].first.Bytes( )[l] )
                              {    diff = True;
                                   break;    }    }
                         if (diff) break;    }    }
               for ( int64_t k1 = i; k1 < j; k1++ )
               {    if ( kmers_plus[k1].second < hb.EdgeObjectCount( ) ) continue;
                    int id = kmers_plus[k1].second - hb.EdgeObjectCount( );
                    if ( found[id] ) continue;
                    for ( int64_t k2 = i; k2 < j; k2++ )
                    {    if ( kmers_plus[k2].second >= hb.EdgeObjectCount( ) ) 
                                   continue;
                         int e = kmers_plus[k2].second, p2 = kmers_plus[k2].third;
                         if ( hb.To( to_left[e] ).nonempty( ) && p2 < hb.K( ) - K ) 
                              continue;
                         vec<int> v(1);
                         v[0] = e;
                         partials[id].push( v, p2, True, 0 );    }    }
               i = j - 1;    }
          for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          {    if ( partials[id].nonempty( ) )
               {    if ( !found[id] ) 
		    found[id] = True;    }    }    }
     // cout << TimeSince(clock3) << " used finding partials" << endl; // XXXXXXXXXXX

     // double clock4 = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     #pragma omp parallel for schedule (dynamic, 1000)
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
     {    const basevector& b = bases[id];
          const qualvector& q = quals[id];
          if ( partials[id].isize( ) > DIVINE_MAX_LOCS ) continue;
          pqsum candidates;
          ExtendReadPlaces1( partials[id], candidates, b, q, hb, 
               to_right, max_partials, 3, 0 );

          /*
          PLACES[id].reserve( candidates.size( ) );
          while ( candidates.size( ) > 0 ) 
	  {   PLACES[id].push_back( candidates.top( ) );
              candidates.pop( );   }    }    }
          */

          if ( candidates.size( ) > 0 )
          {    read_place p = candidates.top( );
               IntVec x;
               for ( int j = 0; j < p.N( ); j++ )
                    x.push_back( p.E(j) );
               (IntVec&) paths[id] = x;
               paths[id].setOffset( p.P( ) );    }    }    

     // cout << TimeSince(clock4) << " used extending" << endl; // XXXXXXXXXXXXXXXXXX

}
