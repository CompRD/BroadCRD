/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Take the reads from Snorgle PHASE=1 and walk through "unambiguous" 20-mers
// in them.

// 2.2 hours for 121,245,220 reads (crd28)

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ReadStack.h"

void BuildStack( const vecbasevector& bases, const vecqualvector& quals,
     const vec< triple<int,int,Bool> >& track, readstack& stack );
void AlignToRef( const basevector& bseq );
void PrintTrack( const vec< triple<int,int,Bool> >& track );
template<int K> Bool AnalyzeNext( const vec< triple<kmer<K>,int,int> >& kmers_plus,
     const int64_t low, const int64_t high, const vecbasevector& bases,
     const vecqualvector& quals, const Bool forward, const int max_alt,
     const int min_mul, ostringstream& out, vec<int>& ids );
template<int K> Bool Advance( const vec< triple<kmer<K>,int,int> >& kmers_plus,
     kmer<K>& x, const vecbasevector& bases, const vecqualvector& quals,
     Bool& forward, vec< triple<int,int,Bool> >& track, ostringstream& out,
     int& its, const int max_alt, const int min_mul, const int max_len, 
     String& seq );

int main( )
{    RunTime( );
     double clock = WallClockTime( );

     // Directory.

     String dir = "/wga/scr4/jaffe/GapToy/51400.newchem";

     // Heuristics.

     const int K = 20;
     const int min_mul = 5;
     const int min_len = 100;
     const int max_len = 1000;
     const int max_alt = 100;
     const int batches = 250;

     // Verbosity.

     Bool print_stack = False;
     Bool print_track = False;
     Bool qlt = True;

     // Special.

     String special =
"CTCCCATTCCCTTCCCCCCTTACCTTTTACCCCCCGCAACTCAGGCTGAACGCCCCCTTAGAGTCTCCTAAGGCCGTCGC"
"CACCAGGAGCCTTTCCCGCCCCCGCGATCTCAGTTTCTATTCCCTCCCCACTCTGCAGGAGCCCCCCAGCTCCCCAGTCT"
"TTACCCACCCTGTCTGGCCCCACCCCCCCTCCCCGCGTTGCTAGGTCCTCCCTCGGCCCTTCCCCCGCCCCCACAAACTG"
"CCCTTCGCAACCCCTTCCCCCCAGTTACCTTTCAGGTCGCCCACTTCCTCAGCCCCGTCCCCCTTTCCCATTGCTCCTCC"
"AGTCCTCCCTTCTCCTGTTGCCCACCCCCAACCTCAGGCATTTACAGGAAGACATTTTGGAGCTACGATCTCTCCTTAGA";
     vecbasevector sp;
     sp.push_back( basevector(special) );
     vec< triple<kmer<K>,int,int> > skmers_plus;
     MakeKmerLookup( sp, skmers_plus );
     vec<kmer<K>> skmers( skmers_plus.size( ) );
     for ( int i = 0; i < skmers_plus.isize( ); i++ )
          skmers[i] = skmers_plus[i].first;

     // Load reads.

     cout << Date( ) << ": loading reads" << endl;
     vec<int> ids;
     BinaryReader::readFile( dir + "/hard_ids", &ids );
     cout << ToStringAddCommas( ids.size( ) ) << " reads" << endl;
     // ids.resize(100000000);
     cout << ToStringAddCommas( ids.size( ) ) << " reads" << endl;

     ids.clear( );

     /*
     ids = {10547430,95492598,95492599,102634588,124133603,141722529,167063884,
          179811585,233633547,282077178,305114558,310184903,325838503,351033980,
          352920328,352920329,365433788,396956579,445269463,586006580,604117940,
          609087827,719627347,856102018};
     */

     ids = { // from RegionIds X=X:70473273-70473820
10547430,37903883,80472101,95492598,102634588,113318980,124133602,125336546,
127286527,141722528,153889594,159637946,160736225,167063884,171845182,179811584,
211370960,219550282,233633546,244117272,244129804,273117088,282077178,305114558,
310184902,316524086,324020363,325838502,342334935,351033980,352920328,362945590,
364614294,365433788,396956578,410048646,445269462,453646710,466701666,493721150,
502222978,524599583,586006580,587355756,587444914,603295767,604117940,609087826,
661594348,696017959,717086889,718146556,719627346,722345201,729219413,759076072,
802214063,827573894,832136544,856102018,866855746,871395155,879096877,879280984
          };

     int ni = ids.size( );
     for ( int j = 0; j < ni; j++ )
     {    int pid = ids[j]/2;
          ids.push_back( 2*pid, 2*pid+1 );    }
     UniqueSort(ids);

     vecbasevector bases;
     bases.Read( dir + "/data/frag_reads_orig.fastb", ids );
     int64_t N = bases.size( );
     cout << "mem usage = " << MemUsageGBString( ) << endl;

     // Make lookup.

     cout << Date( ) << ": making kmer lookup" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup( bases, kmers_plus );
     cout << "mem usage = " << MemUsageGBString( ) << endl;

     // Read and unpack qual scores.

     cout << Date( ) << ": reading quals" << endl;
     VecPQVec pquals;
     pquals.Read( dir + "/data/frag_reads_orig.qualp", ids );
     cout << "mem usage = " << MemUsageGBString( ) << endl;
     cout << Date( ) << ": unpacking quals" << endl;
     vecqualvector quals( pquals.size( ) );
     for ( int64_t id = 0; id < N; id++ )
          pquals[id].unpack( &quals[id] );
     Destroy(pquals);
     cout << "mem usage = " << MemUsageGBString( ) << endl;

     // Search.

     cout << Date( ) << ": searching" << endl;
     vec<int64_t> bstart(batches+1);
     for ( int64_t i = 0; i <= batches; i++ )
          bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t i = 1; i < batches; i++ )
     {    int64_t& s = bstart[i];
          while( s > 0 && kmers_plus[s].first == kmers_plus[s-1].first )
          {    s--;    }    }
     int scount = 0;
     // #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t bi = 0; bi < batches; bi++ )
     for ( int64_t li = bstart[bi]; li < bstart[bi+1]; li++ )
     {    int64_t lj;
          for ( lj = li + 1; lj < bstart[bi+1]; lj++ )
               if ( kmers_plus[lj].first != kmers_plus[li].first ) break;
          if ( lj - li == 1 ) continue;

          // There are two passes.  On the first pass we follow the kmer 
          // forward.  On the second pass we follow it backward.

          for ( int zpass = 1; zpass <= 2; zpass++ )
          {
               Bool forward = ( zpass == 1 );

               Bool fw1 = False, fw2 = False;
               for ( int64_t m = li; m < lj; m++ )
               {    int pos = kmers_plus[m].third;
                    if ( pos == 0 ) fw1 = True;
                    if ( pos == -1 ) fw2 = True;    }

               if ( ( !fw1 && forward ) || ( !fw2 && !forward ) ) continue;

               kmer<K> x = kmers_plus[li].first;
               String seq = x.ToString( );
               if ( !forward ) seq.ReverseComplement( );

               // Track stack.

               vec< triple<int,int,Bool> > track; // (pos,id,fw)

               ostringstream out;
               out << "\n=========================================================="
                    << "==========" << endl;

               int its = 0;
               while(1)
               {    if ( !Advance( kmers_plus, x, bases, quals, forward, track, 
                         out, its, max_alt, min_mul, max_len, seq ) )
                    {    break;    }    }
               basevector bseq(seq);
               bseq.Print( out, "seq" );
               Bool precious = False;
               for ( int j = 0; j <= bseq.isize( ) - K; j++ )
               {    basevector b( bseq, j, K );
                    kmer<K> x(b);
                    if ( BinMember( skmers, x ) ) precious = True;    }
               if ( seq.isize( ) >= min_len ) 
               {    
                    #pragma omp critical
                    {    if ( ++scount % 20000 == 0 || precious )
                         {    cout << out.str( );    

                              // Print track.

                              Sort(track);
                              if (print_track) PrintTrack(track);
     
                              // Build and print stack.

                              readstack stack;
                              BuildStack( bases, quals, track, stack );
                              if (print_stack)
                              {    cout << "stack:\n";
                                   stack.Print(cout);    }

                              if (precious) 
                              {    cout << "PRECIOUS!" << endl;    
                                   {    Ofstream( sout, "moo.fasta" );
                                        bseq.Print( sout, "seq" );    }
                                   if (qlt) AlignToRef(bseq);    }    
                                        }    }    }    }

          li = lj - 1;    }

     // Done.

     cout << "\n" << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem usage = " << PeakMemUsageGBString( ) << endl;
     Scram(0);    }

template<int K> Bool Advance( const vec< triple<kmer<K>,int,int> >& kmers_plus,
     kmer<K>& x, const vecbasevector& bases, const vecqualvector& quals,
     Bool& forward, vec< triple<int,int,Bool> >& track, ostringstream& out,
     int& its, const int max_alt, const int min_mul, const int max_len, String& seq )
{
     // Bound x.

     int64_t low = LowerBound1( kmers_plus, x ), high = UpperBound1( kmers_plus, x );

     // Append to track.

     for ( int64_t m = low; m < high; m++ )
     {    int id = kmers_plus[m].second, pos = kmers_plus[m].third;
          Bool fw = ( pos >= 0 );
          if ( !fw ) pos = -pos - 1;
          int n = bases[id].size( );
          triple<int,int,Bool> t;
          if ( (forward && fw) || (!forward && !fw) )
               t = make_triple( -pos + its, id, True );
          else t = make_triple( K - n + pos + its, id, False );
          if ( !Member( track, t ) ) track.push_back(t);    }

     // Proceed.

     String v = x.ToString( );
     if ( !forward ) v.ReverseComplement( );
     out << "\n" << ++its << ".  " << v << " [" << high-low << "]\n";
     vec<int> ids;
     if ( !AnalyzeNext( kmers_plus, low, high, bases, quals, 
          forward, max_alt, min_mul, out, ids ) )
     {    return False;    }
     String s = x.ToString( );
     if (forward) s = s.substr( 1, K-1 ) + as_base( ids[0] );
     else s = as_base( ids[0] ) + s.substr( 0, K-1 );
     seq.push_back( as_base( forward ? ids[0] : 3 - ids[0] ) );
     if ( s == x.ToString( ) ) return False;
     basevector b(s);
     x = b;
     kmer<K> xrc = x;
     xrc.ReverseComplement( );
     if ( xrc < x ) 
     {    x = xrc;
          forward = !forward;    }
     if ( x == kmers_plus[low].first ) return False; // ?????????????????????
     if ( its == max_len ) return False;
     return True;    }

void BuildStack( const vecbasevector& bases, const vecqualvector& quals,
     const vec< triple<int,int,Bool> >& track, readstack& stack )
{    int low = 1000000000, high = 0;
     for ( int i = 0; i < track.isize( ); i++ )
     {    low = Min( low, track[i].first );
          high = Max( high, track[i].first
               + bases[ track[i].second ].isize( ) );    }
     int n = track.size( ), k = high - low;
     stack.Initialize( n, k );
     for ( int j = 0; j < n; j++ )
     {    int id = track[j].second, offset = track[j].first;
          basevector b = bases[id];
          qualvector q = quals[id];
          Bool fw = track[j].third;
          if ( !fw )
          {    b.ReverseComplement( );
               q.ReverseMe( );    }
          for ( int p = 0; p < b.isize( ); p++ )
          {    stack.SetBase( j, p+offset-low, b[p] );
               stack.SetQual( j, p+offset-low, q[p] );    }
          stack.SetOffset( j, offset );
          stack.SetLen( j, b.size( ) );
          stack.SetId( j, id );
          stack.SetRc2( j, !fw );    }    }

void AlignToRef( const basevector& bseq )
{    SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=moo.fasta "
          "L=/wga/scr4/bigrefs/human19/genome.lookup TARGETS_TO_PROCESS=22 "
          "QUIET=True NH=True VISUAL=True SMITH_WAT=True" );    }

void PrintTrack( const vec< triple<int,int,Bool> >& track )
{    cout << "\ntrack:\n";
     for ( int i = 0; i < track.isize( ); i++ )
     {    cout << "[" << i+1 << "] " << track[i].second << " @ " << track[i].first 
               << " " << ( track[i].third ? "fw" : "rc" ) << endl;    }
     vec<int> ids;
     for ( int i = 0; i < track.isize( ); i++ )
          ids.push_back( track[i].second );
     Sort(ids);
     for ( int i = 1; i < ids.isize( ); i++ )
     {    if ( ids[i] == ids[i-1] )
          {    cout << "Warning: " << ids[i] << " appears more than once."
                    << endl;    }    }    }

template<int K> Bool AnalyzeNext( const vec< triple<kmer<K>,int,int> >& kmers_plus,
     const int64_t low, const int64_t high, const vecbasevector& bases,
     const vecqualvector& quals, const Bool forward, const int max_alt,
     const int min_mul, ostringstream& out, vec<int>& ids )
{
     vec<int> count(4,0), qcount(4,0);
     ids = vec<int>( 4, vec<int>::IDENTITY );
     for ( int64_t i = low; i < high; i++ )
     {    int id = kmers_plus[i].second, pos = kmers_plus[i].third;
          Bool fw = ( pos >= 0 );
          if ( !fw ) pos = -pos - 1;
          int np = ( forward == fw ? pos+K : pos-1 );
          if ( np == -1 || np == bases[id].isize( ) ) continue;
          char y = ( fw ? bases[id][np] : 3 - bases[id][np] );
          int q = quals[id][np];
          count[y]++;
          qcount[y] += q;    }
     ReverseSortSync( qcount, count, ids );
     for ( int j = 0; j < 4; j++ )
     {    if ( count[j] == 0 ) continue;
          if (forward) out << as_base( ids[j] ); 
          else out << as_base( 3 - ids[j] );
          out << " n=" << count[j] << " q=" << qcount[j] << endl;    }
     if ( qcount[0] == 0 || qcount[1] >= max_alt ) return False;
     if ( qcount[0] < min_mul * qcount[1] ) return False;
     if ( Sum(count) == 1 ) return False;
     return True;    }
