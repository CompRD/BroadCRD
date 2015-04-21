///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// HQHBV.  Given a set of reads (bases and quality scores), truncate the reads
// at the first base having quality < 20, and then build a HyperBasevector from
// the remaining bases.
//
// Tested on 'standard' human data set.  Run time exclusive of i/o.
// crd26 run time 2.06 hours, peak mem 484.1 GB: see HQHBV.log.
//
// Notes.  
// - None of the byte5 stuff is used.
// - ParallelSort uses 2x the memory of the object.  Could reduce memory usage
//   by sorting in a small number of chunks and then merging the chunks on the fly.
// - passes_to_use is probably not computed correctly.
// - K must be odd, but the HyperBasevector is constructed using K-1.
//
// Useful HEAD values (add /frag_reads_orig):
// - tmp.xxx
// - /wga/scr4/human_data/NA12878_H01UJ/test-sets-2/reads.chr15,16,17,18,19,20,21,22
// - /wga/scr4/human_data/NA12878_H01UJ/test-sets-2/reads.chr15
// - /wga/scr4/tsharpe/qgraph/whole.genome.

#define _GLIBCXX_PARALLEL

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/KmerCount.h"
#include "random/Random.h"
#include "random/Shuffle.h"
#include "system/SortInPlace.h"

// A byte5 is a five-byte integer.

class byte5 {

     public:

     byte5( ) { }

     byte5( int64_t val )
     {    unsigned char* p = (unsigned char*) &val;
          for ( int j = 0; j < 5; j++ )
               x_[j] = p[j];    }
     
     int64_t Val( ) const
     {    int64_t y;
          unsigned char* p = (unsigned char*) &y;
          for ( int j = 0; j < 5; j++ )
               p[j] = x_[j];
          if ( x_[4] == 0 || x_[4] == 255 )
          {    for ( int j = 5; j < 8; j++ )
                    p[j] = x_[4];    }
          else
          {    for ( int j = 5; j < 8; j++ )
                    p[j] = 0;    }
          return y;    }

     friend Bool operator==( const byte5& x, const byte5& y )
     {    for ( int j = 0; j < 5; j++ )
               if ( x.x_[j] != y.x_[j] ) return False;
          return True;    }

     friend Bool operator<( const byte5& x, const byte5& y )
     {    return x.Val( ) < y.Val( );    }

     friend Bool operator<( const byte5& x, const int64_t& y )
     {    return x.Val( ) < y;    }

     friend Bool operator>=( const byte5& x, const byte5& y )
     {    return x.Val( ) >= y.Val( );    }

     friend byte5 operator-( const byte5& x, const byte5& y )
     {    return x.Val( ) - y.Val( );    }

     friend byte5 operator-( const byte5& x, const int64_t& y )
     {    return x.Val( ) - y;    }

     friend byte5 operator-( const int64_t& x, const byte5& y )
     {    return x - y.Val( );    }

     private:

     unsigned char x_[5];

};

// A vec4 is a vector that can hold up to four T-integers and is optimized to hold
// <= 1 integers.  Not thread-safe.  Given the value of over_floor, this should
// be used with T of size 5 bytes or larger.  Don't try to stick integers larger
// than over_floor in this.  Don't try to put more than four elements in a vec4.

template<class T> class vec4 {

     private:

     static vec< quad<T,T,T,T> > overflow;
     const int64_t over_floor = 100l * 1000l * 1000l * 1000l;

     public:

     T val;

     vec4( ) { val = -1; }
     Bool empty( ) const { return val < 0; }
     Bool nonempty( ) const { return val >= 0; }
     Bool solo( ) const { return val >= 0 && val < over_floor; }

     int size( ) const
     {    if ( val < 0 ) return 0;
          else if ( val < over_floor ) return 1;
          else if ( overflow[val-over_floor].third < 0 ) return 2;
          else if ( overflow[val-over_floor].fourth < 0 ) return 3;
          else return 4;    }

     void push_back( const T& x )
     {    if ( val < 0 ) val = x;
          else if ( val < over_floor )
          {    overflow.push( val, x, -1, -1 );
               val = over_floor + overflow.size( ) - 1;    }
          else if ( overflow[val-over_floor].third < 0 ) 
               overflow[val-over_floor].third = x;
          else overflow[val-over_floor].fourth = x;    }

     T operator[ ]( const int i ) // unchecked!
     {    if ( i == 0 )
          {    if ( val < over_floor ) return val;
               else return overflow[val-over_floor].first;    }
          if ( i == 1 ) return overflow[val-over_floor].second;
          if ( i == 2 ) return overflow[val-over_floor].third;
          return overflow[val-over_floor].fourth;    };

     void destroy( ) { Destroy(overflow); }

};

// Separately implemented vec4_byte5 since I couldn't get vec4<byte5> to work.

class vec4_byte5 {

     private:

     static vec< quad<byte5,byte5,byte5,byte5> > overflow;
     const int64_t over_floor = 100l * 1000l * 1000l * 1000l;

     public:

     byte5 val;

     vec4_byte5( ) { val = -1; }
     Bool empty( ) const { return val.Val( ) < 0; }
     Bool nonempty( ) const { return val.Val( ) >= 0; }
     Bool solo( ) const { return val.Val( ) >= 0 && val.Val( ) < over_floor; }

     int size( ) const
     {    if ( val.Val( ) < 0 ) return 0;
          else if ( val.Val( ) < over_floor ) return 1;
          else if ( overflow[val.Val()-over_floor].third < 0 ) return 2;
          else if ( overflow[val.Val()-over_floor].fourth < 0 ) return 3;
          else return 4;    }

     void push_back( const byte5& x )
     {    if ( val.Val( ) < 0 ) val = x;
          else if ( val.Val( ) < over_floor )
          {    overflow.push( val, x, -1, -1 );
               val = over_floor + overflow.size( ) - 1;    }
          else if ( overflow[val.Val()-over_floor].third < 0 ) 
               overflow[val.Val()-over_floor].third = x;
          else overflow[val.Val()-over_floor].fourth = x;    }

     byte5 operator[ ]( const int i ) // unchecked!
     {    if ( i == 0 )
          {    if ( val.Val( ) < over_floor ) return val;
               else return overflow[val.Val()-over_floor].first;    }
          if ( i == 1 ) return overflow[val.Val()-over_floor].second;
          if ( i == 2 ) return overflow[val.Val()-over_floor].third;
          return overflow[val.Val()-over_floor].fourth;    };

     void destroy( ) { Destroy(overflow); }

};

template<> vec< quad<int64_t,int64_t,int64_t,int64_t> > vec4<int64_t>::overflow
     = vec< quad<int64_t,int64_t,int64_t,int64_t> >( );

template<> vec< quad<byte5,byte5,byte5,byte5> > vec4<byte5>::overflow
     = vec< quad<byte5,byte5,byte5,byte5> >( );

vec< quad<byte5,byte5,byte5,byte5> > vec4_byte5::overflow
     = vec< quad<byte5,byte5,byte5,byte5> >( );

// KmerClass.  Consider the K-mer s starting at position p on basevector b.
// Consider the front and back d-mers x and y of s.  Return the minimum of x and y'.

int KmerClass( const basevector& b, const int p, const int K, const int d )
{    int x = 0, y = 0;
     for ( int i = 0; i < d; i++ )
     {    x = x << 2;
          x += b[p+i];
          y = y << 2;
          y += 3 - b[p+K-1-i];    }
     return Min( x, y );    }

const int K = 61;
vec< kmer<K> > bkmers;

template<int K> unsigned char GetBase( const kmer<K>& x, const int i )
{    int d = i % 4;
     unsigned char c = x.Bytes( )[i/4];
     return ( c >> (2*d) ) & 3;    }

Bool cmp_bpairs( const int64_t i1, const int64_t i2 )
{    int d1 = ( i1 >= 0 ? 0 : 1 );
     int d2 = ( i2 >= 0 ? 0 : 1 );
     int64_t j1 = i1, j2 = i2;
     if ( j1 < 0 ) j1 = -i1-1;
     if ( j2 < 0 ) j2 = -i2-1;
     for ( int j = 0; j < K-1; j++ )
     {    unsigned char c1 = GetBase( bkmers[j1], j + d1 );
          unsigned char c2 = GetBase( bkmers[j2], j + d2 );
          if ( c1 < c2 ) return True;
          if ( c1 > c2 ) return False;    }
     return False;    }

int main(int argc, char *argv[])
{    double clock = WallClockTime( );
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD, "looks for IN_HEAD.{fastb,qualb}, "
          "in ParseStringSet format");
     CommandArgument_String_OrDefault_Doc(OUT_HEAD, "", 
          "generates OUT_HEAD.{hbv,fastb}");
     CommandArgument_Int_OrDefault_Doc(MINQ, 7, "minimum quality score");
     CommandArgument_Int_OrDefault_Doc(MIN_FREQ, 3, "minimum kmer frequency");
     CommandArgument_Int_OrDefault_Doc(QUAL_KEY_SIZE, 6, "quality score key size");
     CommandArgument_Double_OrDefault_Doc(QUAL_DIST_CUTOFF, 0.25,
          "cutoff for quality score distribution");
     CommandArgument_Bool_OrDefault(PRETRIMMED, False);
     CommandArgument_Bool_OrDefault_Doc(CLEAN, False, "remove hanging ends");
     CommandArgument_Int_OrDefault_Doc(NPASSES, 8, "number of passes");
     CommandArgument_Bool_OrDefault_Doc(TAIL_FILTER, False, "trim after last good kmer");
     EndCommandArguments;

     // Heuristics affecting results (see also command-line args).

     // const int K = 61;              // kmer size

     // Heuristics affecting computational performance.

     const int kks = 5;             // kmer key size
     const int sample = 500000;     // number of reads for kmer frequency estimate
     const double mem_fudge = 1.2;  // safety factor for memory usage
     const int nbstarts = 100;      // passes for bpairs

     // Load and trim reads.

     cout << Date( ) << ": loading reads" << endl;
     vec<String> heads;
     ParseStringSet( IN_HEAD, heads );
     cout << Date( ) << ": found " << heads.size( ) << " heads" << endl;
     vecbasevector bases;
     for ( int i = 0; i < heads.isize( ); i++ )
          bases.ReadAll( heads[i] + ".fastb", True );

     cout << Date( ) << ": total reads = " << ToStringAddCommas( bases.size( ) ) 
          << endl;
     if ( !PRETRIMMED )
     {    
          // Define distributions

          /*
          vec< vec<int> >  qdist( IPow( 4, QUAL_KEY_SIZE ) );
          for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          for ( int p = 0; p <= bases[id].isize( ) - QUAL_KEY_SIZE; p++ )
          {    int v = 0;
               for ( int j = 0; j < QUAL_KEY_SIZE; j++ )
               {    v *= 4;
                    v += bases[id][p+j];    }
               qdist[v].push_back( quals[id][ p + QUAL_KEY_SIZE - 1 ] );    }
          */
          
          // Define cutoffs.

          /*
          const int cutoff_default = 0;
          vec<int> cutoff( qdist.size( ), cutoff_default );
          for ( int64_t i = 0; i < qdist.jsize( ); i++ )
          {    if ( qdist[i].empty( ) ) continue;
               Sort( qdist[i] );
               int mid = int( round( qdist[i].size( ) * QUAL_DIST_CUTOFF ) );
               cutoff[i] = qdist[i][mid];    }
          */

          int64_t start = 0;
          for ( int i = 0; i < heads.isize( ); i++ )
          {    cout << Date( ) << ": loading quals " << i+1 << endl;
               vecqualvector quals( heads[i] + ".qualb" );
               cout << Date( ) << ": trimming reads" << endl;
               #pragma omp parallel for
               for ( int64_t id = 0; id < (int64_t) quals.size( ); id++ )
               /*
               for ( int p = QUAL_KEY_SIZE - 1; p < bases[id].isize( ); p++ )
               {    int v = 0;
                    for ( int j = 0; j < QUAL_KEY_SIZE; j++ )
                    {    v *= 4;
                         v += bases[id][ p + j - QUAL_KEY_SIZE + 1 ];    }
                    if ( quals[id][p] < Min( 15, Max( 10, cutoff[v] ) ) )
                    {    bases[id].resize(p);
                         break;    }    }    }
               */
               if ( !TAIL_FILTER )
               {    for ( int p = 0; p < bases[start+id].isize( ); p++ )
                    {    if ( quals[id][p] < MINQ )
                         {    bases[start+id].resize(p);
                              break;    }    }    }
               else
               {    int goods = 0;
                    for ( int p = bases[start+id].isize( ) - 1; p >= 0; p-- )
                    {    if ( quals[id][p] < MINQ ) goods = 0;
                         else
                         {    goods++;
                              if ( goods == K )
                              {    bases[start+id].resize( p + K );
                                   break;    }    }    }
                    if ( goods < K ) bases[start+id].resize(0);    }
               start += quals.size( );    }    }
     int64_t total = 0;
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          total += bases[id].size( );
     cout << "mean trimmed read length = " << total / bases.size( ) << endl;

     // Estimate distribution of kmer classes.

     cout << Date( ) << " estimating kmer class distribution" << endl;
     vec<uint64_t> rids;
     Shuffle64( bases.size( ), rids );
     int nclasses = IPow( 4, kks );
     vec<int> freq( nclasses, 0 );
     for ( uint64_t i = 0; i < Min( bases.size( ), (uint64_t) sample ); i++ )
     {    uint64_t id =  rids[i];
          for ( int p = 0; p <= bases[id].isize( ) - K; p++ )
               freq[ KmerClass( bases[id], p, K, kks ) ]++;    }
     PRINT2( Min(freq), Max(freq) );

     // Compute number of passes.  (FOR NOW, ADVISORY!)

     cout << Date( ) << ": computing number of passes" << endl;
     int64_t mem = MemAvailable( );
     cout << "see " << ToStringAddCommas(mem) << " available bytes" << endl;
     int64_t total_kmers = 0;
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          if ( bases[id].isize( ) >= K ) total_kmers += bases[id].isize( ) - K + 1;
     int64_t total_kmer_bytes = total_kmers * (K+3)/4;
     cout << "see " << ToStringAddCommas(total_kmer_bytes) << " bytes in kmers"
          << endl;
     int passes_to_use = int( ceil( ( total_kmer_bytes * mem_fudge ) / mem ) );
     PRINT(passes_to_use);

     // Define kmer bins.  Trying to split into relatively even piles.

     cout << Date( ) << ": defining kmer bins" << endl;
     ForceAssertLe( NPASSES, nclasses );
     vec<int> kbin( nclasses, -1 ), fids( nclasses, vec<int>::IDENTITY );
     vec<int64_t> inbin( NPASSES, 0 );
     ReverseSortSync( freq, fids );
     for ( int j = 0; j < NPASSES; j++ )
     {    kbin[ fids[j] ] = j;
          inbin[j] = freq[j];    }
     for ( int j = NPASSES; j < nclasses; j++ )
     {    int m = 1000000000, mi = 0;
          for ( int l = 0; l < NPASSES; l++ )
          {    if ( inbin[l] < m ) 
               {    m = inbin[l];
                    mi = l;    }    }
          kbin[ fids[j] ] = mi;
          inbin[mi] += freq[j];    }
     PRINT2( Min(inbin), Max(inbin) );
     cout << "max to min bin ratio = " 
          <<  double(Max(inbin)) / double(Min(inbin)) << endl;

     // Go through passes.

     cout << Date( ) << ": starting main passes" << endl;
     vec<int64_t> kfreq( 100, 0 );
     int64_t kbig = 0;
     vec< kmer<K> > passing;
     for ( int pass = 0; pass < NPASSES; pass++ )
     {    cout << "starting pass " << pass+1 << " of " << NPASSES << endl;

          // Index kmers in the reads.

          vec< kmer<K> > kmers;
          vec<int> counts( bases.size( ), 0 );
          #pragma omp parallel for
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    const basevector& u = bases[i];
               int count = 0;
               for ( int j = 0; j <= u.isize( ) - K; j++ )
               {    int c = kbin[ KmerClass( u, j, K, kks ) ];
                    if ( c == pass ) counts[i]++;    }    }
          vec<int64_t> starts( bases.size( ) + 1 );
          starts[0] = 0;
          for ( size_t i = 0; i < bases.size( ); i++ )
               starts[i+1] = starts[i] + counts[i];
          kmers.resize( starts.back( ) );

          // Find the kmers.

          #pragma omp parallel for
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    const basevector& u = bases[i];
               int count = 0;
               kmer<K> x, xrc;
               for ( int j = 0; j <= u.isize( ) - K; j++ )
               {    int c = kbin[ KmerClass( u, j, K, kks ) ];
                    if ( c != pass ) continue;
                    int64_t r = starts[i] + count++;
                    x.SetToSubOf( u, j );
                    xrc = x;
                    xrc.ReverseComplement( );   
                    if ( x <= xrc ) kmers[r] = x;
                    else kmers[r] = xrc;    }    }
          ParallelSort(kmers);

          // Find the passing kmers.

          int groups = omp_get_max_threads( );
          vec<int64_t> kstarts(groups+1);
          for ( int i = 0; i <= groups; i++ )
               kstarts[i] = ( i * kmers.size( ) ) / groups;
          for ( int i = groups - 1; i >= 0; i-- )
          {    while( kstarts[i] > 0 
                    && kmers[ kstarts[i] ] == kmers[ kstarts[i] - 1 ] )
               {    kstarts[i]--;    }    }
          vec< vec< kmer<K> > > this_passing(groups);
          #pragma omp parallel for
          for ( int p = 0; p < groups; p++ )
          {    vec<int64_t> this_kfreq( 100, 0 );
               int64_t this_kbig = 0;
               for ( int64_t i = kstarts[p]; i < kstarts[p+1]; i++ )
               {    int64_t j = kmers.NextDiff(i);
                    int64_t n = j - i;
                    if ( n < kfreq.isize( ) ) this_kfreq[n]++;
                    else this_kbig++;
                    if ( n >= MIN_FREQ ) this_passing[p].push_back( kmers[i] );
                    i = j - 1;    }    }

          // Copy passing kmers.

          vec<int64_t> pstarts(groups+1);
          pstarts[0] = passing.size( );
          for ( int j = 1; j <= groups; j++ )
               pstarts[j] = pstarts[j-1] + this_passing[j-1].jsize( );
          passing.resize( pstarts.back( ) );
          #pragma omp parallel for
          for ( int j = 0; j < groups; j++ )
          {    if ( this_passing[j].nonempty( ) )
               {    SafeMemcpy( &passing[pstarts[j]], &this_passing[j][0],
                         this_passing[j].size( ) * sizeof( kmer<K> ) );
                         }    }    }
     bases.destroy( ); // done with bases
     // for ( int i = 1; i < kfreq.isize( ); i++ )
     //      PRINT2( i, kfreq[i] );
     // PRINT(kbig);
     cout << ToStringAddCommas( passing.size( ) ) << " passing kmers" << endl;
     cout << "peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;

     // For now, sort.  Should be possible to eliminate, but doesn't take much time.

     cout << Date( ) << ": sorting passing" << endl;
     ParallelSort(passing);

     // Create bkmers.  

     cout << Date( ) << ": creating bkmers" << endl;
     bkmers.resize( 2 * passing.size( ) );
     #pragma omp parallel for
     for ( int64_t i = 0; i < passing.jsize( ); i++ )
     {    bkmers[ 2*i ] = passing[i];
          bkmers[ 2*i + 1 ] = passing[i];
          bkmers[ 2*i + 1 ].ReverseComplement( );    }
     cout << "peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;
     Destroy(passing); // done with this

     // Set up to find backs and nexts.

     cout << Date( ) << ": pushing" << endl;
     vec<int64_t> bpairs( bkmers.size( ) * 2 );
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) bkmers.size( ); i++ )
     {    bpairs[ 2*i ] = i;
          bpairs[ 2*i + 1 ] = -i-1;    }
     cout << Date( ) << ": sorting to find nexts" << endl;
     ParallelSort( bpairs, cmp_bpairs );
     cout << "peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;

     // Find backs and nexts, which define the kmers that come before and after
     // a given kmer.  In effect uses constant length vectors of size four.
     // First find starts.

     cout << Date( ) << ": finding nexts and backs" << endl;
     vec< vec4<int64_t> > nexts( bkmers.size( ) ), backs( bkmers.size( ) );
     vec<int64_t> bstarts(nbstarts+1);
     for ( int i = 0; i <= nbstarts; i++ )
          bstarts[i] = ( i * bpairs.jsize( ) ) / nbstarts;
     cout << Date( ) << ": defining bstarts" << endl;
     for ( int p = nbstarts - 1; p >= 1; p-- )
     {    int64_t i = bstarts[p];
          kmer<K-1> x1, x2;
          basevector b;
          int64_t i1 = bpairs[i];
          if ( i1 >= 0 ) 
          {    bkmers[i1].GetBasevector(b);
               x1.SetToSubOf( b, 0 );    }
          else
          {    bkmers[-i1-1].GetBasevector(b);
               x1.SetToSubOf( b, 1 );    }
          for ( int64_t j = i - 1; j >= 0; j-- )
          {    int64_t i2 = bpairs[j];
               if ( i2 >= 0 ) 
               {    bkmers[i2].GetBasevector(b);
                    x2.SetToSubOf( b, 0 );    }
               else
               {    bkmers[-i2-1].GetBasevector(b);
                    x2.SetToSubOf( b, 1 );    }
               if ( x2 != x1 ) break;
               bstarts[p]--;    }    }
     cout << Date( ) << ": defining cstarts" << endl;
     vec< vec<int64_t> > cstarts(nbstarts);
     #pragma omp parallel for
     for ( int p = 0; p < nbstarts; p++ )
     {    for ( int64_t i = bstarts[p]; i < bstarts[p+1]; i++ )
          {    cstarts[p].push_back(i);
               kmer<K-1> x1, x2;
               basevector b;
               int64_t i1 = bpairs[i];
               if ( i1 >= 0 ) 
               {    bkmers[i1].GetBasevector(b);
                    x1.SetToSubOf( b, 0 );    }
               else
               {    bkmers[-i1-1].GetBasevector(b);
                    x1.SetToSubOf( b, 1 );    }
               int64_t j;
               for ( j = i + 1; j < bpairs.jsize( ); j++ )
               {    int64_t i2 = bpairs[j];
                    if ( i2 >= 0 ) 
                    {    bkmers[i2].GetBasevector(b);
                         x2.SetToSubOf( b, 0 );    }
                    else
                    {    bkmers[-i2-1].GetBasevector(b);
                         x2.SetToSubOf( b, 1 );    }
                    if ( x2 != x1 ) break;    }
               i = j - 1;    }    }
     cout << "peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;

     // Now go ahead and build the nexts and backs.  Not parallelized very well.

     cout << Date( ) << ": and now building nexts and backs" << endl;
     for ( int p = 0; p < nbstarts; p++ )
     {    vec< pair<int64_t,int64_t> > N, B;
          for ( int64_t j = 0; j < cstarts[p].jsize( ); j++ )
          {    int64_t a = cstarts[p][j], b;
               if ( j < cstarts[p].jsize( ) - 1 ) b = cstarts[p][j+1];
               else if ( p < nbstarts - 1 ) b = cstarts[p+1][0];
               else b = bpairs.jsize( );
               for ( int64_t k1 = a; k1 < b; k1++ )
               {    if ( bpairs[k1] < 0 ) continue;
                    int64_t n1 = bpairs[k1];
                    for ( int64_t k2 = a; k2 < b; k2++ )
                    {    if ( bpairs[k2] >= 0 ) continue;
                         int64_t n2 = -bpairs[k2] - 1;
                         N.push(n1,n2), B.push(n2,n1);    }    }    }
          ParallelSort(N), ParallelSort(B);
          #pragma omp parallel for
          for ( int64_t j = 0; j < N.jsize( ); j++ )
          {    int64_t n1 = N[j].first, n2 = N[j].second;
               if ( j > 0 && n1 == N[j-1].first ) continue;
               nexts[n1].push_back(n2);    }
          for ( int64_t j = 0; j < N.jsize( ); j++ )
          {    int64_t n1 = N[j].first, n2 = N[j].second;
               if ( j > 0 && n1 == N[j-1].first ) nexts[n1].push_back(n2);    }
          #pragma omp parallel for
          for ( int64_t j = 0; j < B.jsize( ); j++ )
          {    int64_t b1 = B[j].first, b2 = B[j].second;
               if ( j > 0 && b1 == B[j-1].first ) continue;
               backs[b1].push_back(b2);    }
          for ( int64_t j = 0; j < B.jsize( ); j++ )
          {    int64_t b1 = B[j].first, b2 = B[j].second;
               if ( j > 0 && b1 == B[j-1].first ) backs[b1].push_back(b2);    }    }
     cout << "peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;
     Destroy(bpairs), Destroy(cstarts); // done with these

     // Find chains of kmers.  These will define the unibases.

     cout << Date( ) << ": making chains 1" << endl;
     vec< vec<int64_t> > chains;
     vec<Bool> used( bkmers.size( ), False );
     for ( int64_t j = 0; j < (int64_t) bkmers.size( ); j++ )
     {    if ( backs[j].empty( ) && nexts[j].empty( ) )
          {    used[j] = True;
               vec<int64_t> x;
               x.push_back(j);
               chains.push_back(x);    }    }
     PRINT2( bkmers.size( ), chains.size( ) );
     cout << Date( ) << ": making chains 2" << endl;
     #pragma omp parallel for
     for ( int64_t j = 0; j < (int64_t) bkmers.size( ); j++ )
     {    if ( used[j] ) continue;
          if ( !backs[j].solo( )
               || ( backs[j].nonempty( ) && nexts[ backs[j][0] ].size( ) >= 2 ) )
          {    vec<int64_t> x;
               x.push_back(j);
               int64_t s = j;
               ForceAssert( !used[s] );
               used[s] = True;
               while( nexts[s].solo( ) && backs[ nexts[s][0] ].size( ) <= 1 )
               {    s = nexts[s][0];
                    x.push_back(s);
                    ForceAssert( !used[s] );
                    used[s] = True;    }
               #pragma omp critical
               {    chains.push_back(x);    }    }    }
     PRINT2( bkmers.size( ), chains.size( ) );
     Destroy(backs);

     // Find circles.

     cout << Date( ) << ": creating circles" << endl;
     for ( int64_t j = 0; j < bkmers.jsize( ); j++ )
     {    if ( used[j] ) continue;

          // Make circle.

          vec<int64_t> x = {j};
          int64_t s = j;
          used[s] = True;
          while(1)
          {    Assert( nexts[s].size( ) == 1 );
               s = nexts[s][0];
               if ( s == x[0] ) break;
               Assert( !used[s] );
               used[s] = True;
               x.push_back(s);    }

          // Canonicalize.

          vec< vec<int64_t> > c;
          c.push_back(x);
          int64_t m = Min(x);
          int64_t zrc;
          if ( x[0] % 2 == 0 ) zrc = x[0] + 1;
          else zrc = x[0] - 1;
          if ( !Member( x, zrc ) )
          {    vec<int64_t> y;
               for ( int i = x.isize( ) - 1; i >= 0; i-- )
               {    if ( x[i] % 2 == 0 ) zrc = x[i] + 1;
                    else zrc = x[i] - 1;
                    int64_t s = zrc;
                    used[s] = True;
                    y.push_back(s);    }
               c.push_back(y);
               m = Min( m, Min(y) );
               if ( m < Min(x) ) swap( c[0], c[1] );    }
          int p = Position( c[0], m );
          x.clear( );
          for ( int l = p; l < c[0].isize( ); l++ )
               x.push_back( c[0][p] );
          for ( int l = 0; l < p; l++ )
               x.push_back( c[0][p] );

          // Save.

          chains.push_back(x);    
          if ( c.size( ) == 2 )
          {    vec<int64_t> y;
               for ( int i = x.isize( ) - 1; i >= 0; i-- )
               {    if ( x[i] % 2 == 0 ) zrc = x[i] + 1;
                    else zrc = x[i] - 1;
                    y.push_back(zrc);    }
               chains.push_back(y);     }    }
     PRINT2( bkmers.size( ), chains.size( ) );

     // Sort chains.  

     cout << Date( ) << ": sorting chains" << endl;
     sortInPlaceParallel( chains.begin( ), chains.end( ) );
     cout << "peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;

     // Make unibases.

     cout << Date( ) << ": making unibases" << endl;
     vec<basevector> unibases( chains.size( ) );
     #pragma omp parallel for
     for ( int64_t i = 0; i < chains.jsize( ); i++ )
     {    vec<int64_t> x = chains[i];
          x.ReverseMe( ); // bug patch??????????????????????????????????????????????
          bkmers[ x[0] ].GetBasevector( unibases[i] );
          for ( int j = 1; j < x.isize( ); j++ )
          {    basevector b;
               bkmers[ x[j] ].GetBasevector(b);
               unibases[i].push_back( b.back( ) );    }    }
     Destroy(bkmers);
     vec<int> len;
     for ( int64_t i = 0; i < unibases.jsize( ); i++ )
          len.push_back( unibases[i].isize( ) - K + 1 );
     Sort(len);
     cout << ToStringAddCommas( unibases.size( ) ) << " unibases of median size " 
          << Median(len) << " and N50 size " << N50(len) << endl;

     // Define graph structure on unibases.

     cout << Date( ) << ": defining con" << endl;
     vec< triple<int64_t,int64_t,int> > con;
     int64_t nuni = chains.size( );
     for ( int64_t i = 0; i < nuni; i++ )
     {    con.push( chains[i][0], i, 0 );
          int64_t s = chains[i].back( );
          for ( int j = 0; j < nexts[s].size( ); j++ )
               con.push( nexts[s][j], i, 1 );    }
     Destroy(chains), Destroy(nexts);
     // vec4<int64_t> a_vec4;
     // a_vec4.destroy( ); // done with vec4
     ParallelSort(con);
     cout << Date( ) << ": building graph" << endl;
     equiv_rel e( 2 * nuni );
     for ( int64_t i = 0; i < con.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < con.jsize( ); j++ )
               if ( con[j].first != con[i].first ) break;
          for ( int64_t k1 = i; k1 < j; k1++ )
          for ( int64_t k2 = i; k2 < j; k2++ )
          {    if ( con[k1].third == 1 && con[k2].third == 0 )
               {    int64_t u1 = con[k1].second, u2 = con[k2].second;
                    e.Join( 2*u1 + 1, 2*u2 );    }    }
          i = j - 1;    }
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     Sort(reps);
     int64_t nreps = reps.size( );
     vec< vec<int> > from(nreps), to(nreps); 
     vec< vec<int> > from_edge_obj(nreps), to_edge_obj(nreps);
     for ( int64_t u = 0; u < nuni; u++ )
     {    int64_t v1 = BinPosition( reps, e.ClassId( 2*u ) ); 
          int64_t v2 = BinPosition( reps, e.ClassId( 2*u + 1 ) );
          swap( v1, v2 ); // bug patch??????????????????????????????????????????????
          from[v1].push_back(v2), to[v2].push_back(v1);
          from_edge_obj[v1].push_back(u), to_edge_obj[v2].push_back(u);    }
     for ( int64_t v = 0; v < nreps; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     HyperBasevector hb(K-1);
     digraphE<basevector> G( from, to, unibases, to_edge_obj, from_edge_obj );
     ( digraphE<basevector>& )(hb) = G;
     PRINT( hb.EdgeObjectCount( ) );

     // Remove hanging ends.

     if (CLEAN)
     {    cout << Date( ) << ": removing hanging ends" << endl;
          vec<kmer_count> kc( hb.EdgeObjectCount( ) );
          for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
               kc[e].n = hb.EdgeObject(e).isize( ) - K + 1;
          const int max_paths = 100;
          const double junk_ratio = 10.0;
          const int max_del = 1000;
          digraphE<kmer_count> hb_kc( hb, kc );
          RemoveHangingEnds3(hb_kc, &kmer_count::N, max_del, junk_ratio, max_paths);
          vec<int> e_to_delete;
          vec<Bool> used;
          hb_kc.Used(used);
          for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
               if ( !used[e] ) e_to_delete.push_back(e);
          hb.DeleteEdges(e_to_delete);
          hb.RemoveUnneededVertices( );
          hb.RemoveDeadEdgeObjects( );
          vec<int> len;
          for ( int64_t i = 0; i < hb.EdgeObjectCount( ); i++ )
               len.push_back( hb.EdgeObject(i).isize( ) - K + 1 );
          Sort(len);
          cout << ToStringAddCommas( hb.EdgeObjectCount( ) )
               << " edges of N50 size " << N50(len) << endl;    }

     // Write output.

     if ( OUT_HEAD != "" ) 
     {    BinaryWriter::writeFile( OUT_HEAD + ".hbv", hb );
          vecbasevector edges( hb.EdgeObjectCount( ) );
          for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
               edges[e] = hb.EdgeObject(e);
          edges.WriteAll( OUT_HEAD + ".fastb" );    }
     cout << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;    }
