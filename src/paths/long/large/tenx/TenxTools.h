///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef TENX_TOOLS_H
#define TENX_TOOLS_H

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"

void FindTenxNhood( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<vec<vec<vec<int>>>>& lines, const vec<int>& lens,
     const vec< pair<int,int> >& places, const vec<int>& lhits,
     const vec<double>& bc_frac,
     const vec<int64_t>& all, const vec<uint32_t>& bcs, const int nbc,
     const int nbc_total,
     const vec<String>& genome_names, const vec< vec< pair<int,int> > >& aligns );

template<int K> void Match( const int start, const HyperBasevectorX& hb,
     const vecbasevector& tigs, const vec<int>& inv, vecbasevector& bases, 
     const int Kbig, const int max_freq, const int max_frag, 
     vec< pair<int,int> >& places )
{
     cout << Date( ) << ": building kmers" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     double kclock = WallClockTime( );
     vec<int64_t> xstarts;
     xstarts.push_back(0);
     for ( size_t e = 0; e < tigs.size( ); e++ )
     {    const basevector& u = tigs[e];
          int v = hb.ToLeft(e), w = hb.ToRight(e), N = u.isize( );
          int cstart = 0, cstop = N - K;
          if ( hb.To(v).size( ) > 0 && hb.From(v).size( ) > 1 ) cstart = Kbig - K;
          if ( hb.From(w).size( ) > 0 && hb.To(w).size( ) > 1 ) cstop = N - Kbig;
          int n = ( ( u.size( ) > 0 && cstart <= cstop ) ? cstop - cstart + 1 : 0 );
          xstarts.push_back( xstarts.back( ) + n );    }
     int64_t goods = 0;
     vec<int64_t> to_pos( bases.size( ), -1 );
     for ( size_t i = 0; i < bases.size( ); i++ )
     {    if ( bases[i].isize( ) >= K + start ) 
          {    to_pos[i] = goods;
               goods++;    }    }
     kmers_plus.resize( xstarts.back( ) + goods );
     #pragma omp parallel for
     for ( size_t e = 0; e < tigs.size( ); e++ )
     {    const basevector& u = tigs[e];
          if ( u.size( ) == 0 ) continue;
          int v = hb.ToLeft(e), w = hb.ToRight(e), N = u.isize( );
          int cstart = 0, cstop = N - K;
          if ( hb.To(v).size( ) > 0 && hb.From(v).size( ) > 1 ) cstart = Kbig - K;
          if ( hb.From(w).size( ) > 0 && hb.To(w).size( ) > 1 ) cstop = N - Kbig;
          if ( cstart > cstop ) continue;
          kmer<K> x;
          for ( int j = cstart; j <= cstop; j++ )
          {    int64_t r = xstarts[e] + j - cstart;
               x.SetToSubOf( u, j ); 
               kmers_plus[r].first = x;
               kmers_plus[r].second = e; 
               kmers_plus[r].third = j;    }    }
     #pragma omp parallel for
     for ( size_t i = 0; i < bases.size( ); i++ )
     {    if ( bases[i].isize( ) < K + start ) continue;
          kmer<K> x;
          x.SetToSubOf( bases[i], start );
          int64_t r = xstarts.back( ) + to_pos[i];
          kmers_plus[r].first = x;
          kmers_plus[r].second = (int64_t) tigs.size( ) + i; 
          kmers_plus[r].third = start;    }
     ParallelSort(kmers_plus);
     cout << TimeSince(kclock) << " used building kmers" << endl;

     // Find matches.

     cout << Date( ) << ": finding matches" << endl;
     const int64_t batches = 1000;
     vec<int64_t> bstart(batches+1);
     for ( int64_t i = 0; i <= batches; i++ )
          bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
     for ( int64_t i = 1; i < batches; i++ )
     {    int64_t& s = bstart[i];
          while( s > bstart[i-1] && kmers_plus[s].first == kmers_plus[s-1].first )
          {    s--;    }    }
     double mclock = WallClockTime( );
     vec< vec< triple<int,int,int> > > resultsx(batches);
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t bi = 0; bi < batches; bi++ )
     {    vec<int64_t> ok;
          for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j, m;
               for ( j = i + 1; j < bstart[bi+1]; j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( m = i; m < j; m++ )
                    if ( kmers_plus[m].second >= (int) tigs.size( ) ) break;
               if ( 1 <= m - i && m - i <= max_freq )
               {    for ( int64_t k2 = m; k2 < j; k2++ )
                    {    int rid = kmers_plus[k2].second - (int) tigs.size( );
                         int ridp = ( rid % 2 == 0 ? rid + 1 : rid - 1 );
                         if ( m - i == 1 )
                         {    for ( int64_t k1 = i; k1 < m; k1++ )
                              {    int tid = kmers_plus[k1].second;
                                   int offset = kmers_plus[k1].third 
                                        - kmers_plus[k2].third;
                                   resultsx[bi].push( rid, tid, offset );    }    }
                         else if ( places[ridp].first >= 0 )
                         {    ok.clear( );
                              int tidp = places[ridp].first;
                              int offsetp = places[ridp].second;
                              for ( int64_t k1 = i; k1 < m; k1++ )
                              {    int tid = kmers_plus[k1].second;
                                   int offset = kmers_plus[k1].third 
                                        - kmers_plus[k2].third;
                                   if ( tid != inv[tidp] )
                                   {    if ( tigs[tid].isize( ) - Kbig + K
                                             >= offset + max_frag )
                                        {    continue;    }
                                        if ( tigs[tidp].isize( ) - Kbig + K
                                             >= offsetp + max_frag )
                                        {    continue;    }    }
                                   ok.push_back(k1);    }
                              if ( ok.solo( ) )
                              {    int64_t k1 = ok[0];
                                   int tid = kmers_plus[k1].second;
                                   int offset = kmers_plus[k1].third 
                                        - kmers_plus[k2].third;
                                   resultsx[bi].push( 
                                        rid, tid, offset );    }    }    }    }
               i = j - 1;    }    }
     cout << TimeSince(mclock) << " used matching kmers" << endl;

     // Collate results.

     cout << Date( ) << ": combining results" << endl;
     vec<int64_t> starts( batches+1, 0 );
     for ( int j = 0; j < batches; j++ )
          starts[j+1] = starts[j] + (int64_t) resultsx[j].size( );
     vec< triple<int,int,int> > results( starts.back( ) );
     #pragma omp parallel for
     for ( int j = 0; j < batches; j++ )
     {    memcpy( &results[ starts[j] ], &resultsx[j][0],
               resultsx[j].size( ) * sizeof( triple<int,int,int> ) );    }

     // Sort.  This seems unnecessary.

     // cout << Date( ) << ": sorting" << endl;
     // ParallelSort(results);

     // Attempt to push some alignments onto bubble.

     cout << Date( ) << ": pushing bubbles" << endl;
     #pragma omp parallel for
     for ( int i = 0; i < results.isize( ); i++ )
     {    int rid = results[i].first;
          int &e = results[i].second, &pos = results[i].third;
          int w = hb.ToRight(e);
          if ( hb.From(w).size( ) != 2 || hb.To(w).size( ) != 1 ) continue;
          int f1 = hb.IFrom( w, 0 ), f2 = hb.IFrom( w, 1 );
          int o = pos - hb.Kmers(e);
          const int min_match = 10;
          Bool match1 = True, match2 = True;
          for ( int j = 0; j < min_match; j++ )
          {    int fpos = hb.K( ) - 1 + j;
               int rpos = hb.K( ) - 1 - o + j;
               if ( rpos < 0 || rpos >= bases[rid].isize( ) )
               {    match1 = match2 = False;
                    break;    }
               if ( fpos < 0 || fpos >= hb.Bases(f1) ) match1 = False;
               if ( fpos < 0 || fpos >= hb.Bases(f2) ) match2 = False;
               if ( bases[rid][rpos] != hb.O(f1)[fpos] ) match1 = False;
               if ( bases[rid][rpos] != hb.O(f2)[fpos] ) match2 = False;
               if ( !match1 && !match2 ) break;    }
          if (match1)
          {    pos -= hb.Kmers(e);
               e = f1;    }
          if (match2)
          {    pos -= hb.Kmers(e);
               e = f2;    }    }

     // Mark places.
 
     cout << Date( ) << ": marking" << endl;
     for ( int i = 0; i < results.isize( ); i++ )
     {    places[ results[i].first ] 
               = make_pair( results[i].second, results[i].third );    
          bases[ results[i].first ].resize(0);    }    }

#endif
