///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "paths/BigMapTools.h"
#include "paths/Jumpster.h"

template<int K> void MakeKmerLookupJ( const vecbasevector& tigs,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    kmer<K> x, xrc;
               int64_t r = starts[i] + j;
               x.SetToSubOf( u, j ); 
               kmers_plus[r].first = x;
               kmers_plus[r].second = i; 
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);    }

void Jumpster( const int K, const vecbasevector& unibases, 
     const vec< vec<int> >& nexts, const vecbasevector& jumps, 
     const vecqualvector& quals, vec<basevector>& bridges, const Bool verbose, 
     const vecbasevector& genome2, const int LG, 
     const VecIntPairVec& Glocs )

{
     // Heuristics.

     const int L = 40;
     const int flank = 5;
     const int min_good = 5;
     const int max_next = 2;
     const int max_cycles = 1000;

     // Form table of kmers in the jumps.

     if (verbose) cout << Date( ) << ": building kmer table" << endl;
     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookupJ( jumps, kmers_plus );
     cout << Date( ) << ": building kmers" << endl;
     vec< kmer<L> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;

     // Go through the unibases.

     if (verbose) cout << Date( ) << ": traversing unibases" << endl;
     for ( size_t u = 0; u < unibases.size( ); u++ )
     {    
          // Consider only dead ends.

          if ( nexts[u].nonempty( ) ) continue;
          if (verbose)
          {    cout << "\nu = " << u << ", nkmers = " 
                    << unibases[u].isize( ) - K + 1 << endl;    }

          // Define start.

          basevector start;
          start.SetToSubOf( unibases[u], unibases[u].isize( ) - K, L );
          kmer<L> cur(start);
          basevector full(start);
          Bool enable_repeat_readthrough = False;
          Bool entering_repeat = False;

          // Walk right.

          int count(1);
          while(1)
          {    if (verbose)
               {    cout << "\n" << count << ". looking forward from " 
                         << cur.ToString( ) << endl;    }
               basevector b;
               cur.GetBasevector(b);
               vec<basevector> exts0;
               vec<qualvector> exts0q;
               int64_t low = LowerBound(kmers, cur), high = UpperBound(kmers, cur);
               for ( int64_t j = low; j < high; j++ )
               {    int id = kmers_plus[j].second, pos = kmers_plus[j].third;
                    const basevector& r = jumps[id];
                    basevector f;
                    if ( pos + L + flank > r.isize( ) ) continue;
                    f.SetToSubOf( r, pos+L, r.isize( ) - (pos+L) );
                    exts0.push_back(f);    
                    qualvector q;
                    q.SetToSubOf( quals[id], pos+L, r.isize( ) - (pos+L) );
                    exts0q.push_back(q);    }
               kmer<L> cur_rc(cur);
               cur_rc.ReverseComplement( );
               low = LowerBound( kmers, cur_rc ), high = UpperBound( kmers, cur_rc );
               for ( int64_t j = low; j < high; j++ )
               {    int id = kmers_plus[j].second, pos = kmers_plus[j].third;
                    const basevector& r = jumps[id];
                    basevector f;
                    if ( pos < flank ) continue;
                    f.SetToSubOf( r, 0, pos );
                    f.ReverseComplement( );
                    exts0.push_back(f);    
                    qualvector q;
                    q.SetToSubOf( quals[id], 0, pos );
                    q.ReverseMe( );
                    exts0q.push_back(q);    }

               // At this point 'exts0' are the read extensions beyond the L-mer, of
               // size at least 'flank'.

               int wcount = 0;
               while(1)
               {    if ( ++wcount > 1 && !enable_repeat_readthrough ) break;
                    if ( exts0.empty( ) ) break;

                    // Get just the flank.

                    vec<basevector> exts(exts0);
                    vec<qualvector> extsq(exts0q);
                    for ( int i = 0; i < exts.isize( ); i++ )
                    {    exts[i].resize(flank);
                         extsq[i].resize(flank);    }

                    // Condense to uniques.

                    vec<basevector> extsx;
                    vec< vec<qualvector> > extsxq;
                    vec<int> extsx_mult;
                    SortSync( exts, extsq );
                    for ( int i = 0; i < exts.isize( ); i++ )
                    {    int j = exts.NextDiff(i);
                         extsx.push_back( exts[i] );
                         extsx_mult.push_back(j-i);
                         vec<qualvector> qs;
                         for ( int k = i; k < j; k++ )
                              qs.push_back( extsq[k] );
                         extsxq.push_back(qs);
                         i = j - 1;    }

                    // Vote.  This normalizes counts using quality scores.

                    ReverseSortSync( extsx_mult, extsx, extsxq );
                    vec<int> topq;
                    for ( int j = 0; j < extsxq[0].isize( ); j++ )
                         topq.push_back( extsxq[0][j].back( ) );
                    double mq = Mean(topq);
                    vec<double> extsx_dmult;
                    for ( int l = 0; l < extsxq.isize( ); l++ )
                    {    vec<int> topq;
                         for ( int j = 0; j < extsxq[l].isize( ); j++ )
                              topq.push_back( extsxq[l][j].back( ) );
                         double mlq = Mean(topq);
                         extsx_dmult.push_back( extsx_mult[l] * mlq/mq );    }
                    ReverseSortSync( extsx_dmult, extsx_mult, extsx, extsxq );
                    if (verbose)
                    {    for ( int j = 0; j < extsx_mult.isize( ); j++ )
                         {    cout << "[" << wcount << "] see " 
                                   << setiosflags(ios::fixed) << setprecision(1)
                                   << extsx_dmult[j] << resetiosflags(ios::fixed)
                                   << "<" << extsx_mult[j] 
                                   << "> x " << extsx[j] << " (";
                                   for ( int l = 0; l < extsxq[j].isize( ); l++ )
                                   {    if ( l > 0 ) cout << " ";
                                        cout << int( extsxq[j][l].back( ) );    }
                                   cout << ")" << endl;    }    }

                    // Do we lack enough coverage?

                    if ( extsx.empty( ) ) break;
                    if ( extsx_mult[0] < min_good ) break;

                    // Does it look like we're at a repeat?  

                    Bool repeat = False;
                    if ( extsx.size( ) > 1 )
                    {    if ( extsx_dmult[1] > max_next )
                              repeat = True;
                         if ( extsx_dmult[0] < min_good * extsx_dmult[1] ) 
                         {    repeat = True;    }    }

                    // If we're at a repeat, enable repeat readthrough and back up.

                    if (repeat)
                    {    if ( full.isize( ) > L && !enable_repeat_readthrough )
                         {    enable_repeat_readthrough = True;
                              entering_repeat = True;
                              basevector b2(L);
                              b2.Set( 0, full[ full.isize( ) - L - 1 ] );
                              for ( int j = 0; j < b.isize( ) - 1; j++ )
                                   b2.Set( j+1, b[j] );
                              b = b2;
                              cur.Set(b);
                              full.resize( full.isize( ) - 1 );
                              count--;    }
                         break;    }

                    // Success!  Save a base.

                    basevector bnew(L);
                    for ( int j = 1; j < L; j++ )
                         bnew.Set( j-1, b[j] );
                    bnew.Set( L-1, extsx[0][0] );
                    cur.Set(bnew);
                    cur.GetBasevector(b);
                    full.push_back( extsx[0][0] );

                    // Far enough?

                    if ( ++count == max_cycles ) break;    

                    // Advance on this readstack.

                    vec<Bool> to_delete( exts0.size( ), False );
                    for ( int j = 0; j < exts0.isize( ); j++ )
                    {    if ( !exts0[j].ToString( ).
                              Contains( extsx[0].ToString( ), 0 ) ) 
                         {    to_delete[j] = True;    }
                         else if ( exts0[j].isize( ) < flank + 1 ) 
                         {    to_delete[j] = True;    }
                         else
                         {    exts0[j].SetToSubOf( 
                                   exts0[j], 1, exts0[j].isize( ) - 1 ); 
                              exts0q[j].SetToSubOf( exts0q[j], 1, 
                                   (int) exts0q[j].size( ) - 1 );    }    }
                    EraseIf( exts0, to_delete );
                    EraseIf( exts0q, to_delete );    }

               // Far enough?

               if ( enable_repeat_readthrough && !entering_repeat && wcount == 2 )
                    break;
               if ( enable_repeat_readthrough && !entering_repeat )
                    enable_repeat_readthrough = False;
               if ( ( wcount == 1 && !entering_repeat ) || count == max_cycles ) 
                    break;
               entering_repeat = False;    }

          // Save bridge.

          bridges.push_back(full);
          if (verbose)
          {    cout << "\n";
               full.Print( cout, "bridge" );
               cout << "\n";    }    
          if ( genome2.size( ) > 0 )
          {    if ( full.isize( ) < K ) cout << "Tiny, ignoring\n";
               else
               {    vec<placementy> places
                         = FindGenomicPlacementsY( 0, full, LG, genome2, Glocs );    
                    if ( places.nonempty( ) ) cout << "Has perfect placement.\n";
                    else
                    {    cout << "IMPERFECT!\n";    }    }    }    }    }
